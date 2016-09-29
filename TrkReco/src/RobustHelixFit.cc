//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.cc,v 1.12 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
//
// the following has to come before other BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "TrkReco/inc/RobustHelixFit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
// root
#include "TH1F.h"
// C++
#include <vector>
#include <string>
#include <math.h>
using CLHEP::Hep3Vector;
using namespace std;
namespace mu2e 
{
  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<HelixHit,HelixHit,bool> {
    bool operator()(HelixHit const& p1, HelixHit const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  
  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("UseFlag",vector<string>{"Outlier","OtherBackground"})),
    _mindelta(pset.get<double>("minDelta",5000.0)),
    _minnhit(pset.get<unsigned>("minNHit",10)),
    _lambda0(pset.get<double>("lambda0",1.0)),
    _lstep(pset.get<double>("lstep",0.2)),
    _minlambda(pset.get<double>("minlambda",0.01)),
    _maxniter(pset.get<unsigned>("maxniter",50)),
    _minzsep(pset.get<double>("minzsep",100.0)),
    _maxzsep(pset.get<double>("maxzsep",700.0)),
    _mindphi(pset.get<double>("mindphi",0.25)),
    _maxdphi(pset.get<double>("maxdphi",2.5)),
    _mindist(pset.get<double>("mindist",50.0)),
    _maxdist(pset.get<double>("maxdist",500.0)),
    _rmin(pset.get<double>("minR",200.0)),
    _rmax(pset.get<double>("maxR",500.0)),
    _tdmin(pset.get<double>("minAbsTanDip",0.3)),
    _tdmax(pset.get<double>("maxAbsTanDip",2.0)),
    _filterxy(pset.get<bool>("filterxy",true)),
    _filterz(pset.get<bool>("filterz",true)),
    _stereoinit(pset.get<bool>("stereoinit",false)),
    _stereofit(pset.get<bool>("stereofit",false)),
    _targetinit(pset.get<bool>("targetinit",true)),
    _targetinter(pset.get<bool>("targetintersect",true)),
    _targetradius(pset.get<double>("targetradius",75.0)), // target radius: include some buffer
    _trackerradius(pset.get<double>("trackerradius",700.0)), // tracker out radius; include some buffer
    _rwind(pset.get<double>("RadiusWindow",10.0)), // window for calling a point to be 'on' the helix
    _rout(pset.get<double>("OutlierRadius",40.0)), 
    _pout(pset.get<double>("OutlierDPhi",2.0)), 
    _helicity(pset.get<int>("Helicity",Helicity::unknown))
  {
    if(_helicity._value == Helicity::unknown){
      throw cet::exception("RECO")<<"mu2e::RobustHelix: Invalid Helicity specified"<< std::endl;
    }
    if(_stereofit)_useflag = StrawHitFlag(StrawHitFlag::stereo);
    if(_helicity._value > 0){
      _smin = 1.0/(_rmax*_tdmax);
      _smax = 1.0/(_rmin*_tdmin);
    } else {
      _smax = -1.0/(_rmax*_tdmax);
      _smin = -1.0/(_rmin*_tdmin);
    }
  }

  RobustHelixFit::~RobustHelixFit()
  {}

  void
  RobustHelixFit::findHelix(HelixSeed& myseed) {
    RobustHelix& myhel = myseed._helix;
// Convert the hit positions to what we need for helix finding
    HelixHitCollection& hhits = myseed._hhits;
// filter by geometry
    if(_filterxy)filterDist(hhits);
    if(hitCount(hhits) >= _minnhit){
      myseed._status.merge(TrkFitFlag::hitsOK);
	// solve for the circle parameters
      if(initCircle(hhits,myhel)){
	if(findXY(hhits,myhel)){
	  myseed._status.merge(TrkFitFlag::circleOK);
	  if(findZ(hhits,myhel)) {
	    // set the success
	    myseed._status.merge(TrkFitFlag::helixOK);
	  }
	}
      }
    }
  }

  bool RobustHelixFit::findXY(HelixHitCollection& hhits,RobustHelix& myhel) {
    double rmed, age;
    Hep3Vector center = myhel.center();
    bool changed(true);
    unsigned niter(0);
    while(niter < _maxniter && changed && hitCount(hhits) > _minnhit){
      findCenterAGE(hhits,center,rmed,age);
      if(_filterxy)
	filterXY(hhits,center,rmed,changed);
      else
	changed = false;
      niter++;
    }
    myhel.center() = center;
    myhel.radius() = rmed;
    return true;
  }

  bool
  RobustHelixFit::findCenterAGE(HelixHitCollection const& hhit,Hep3Vector& center, double& rmed, double& age) {
// this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
// initialize step
    double lambda = _lambda0;
// find median and AGE for the initial center
    findAGE(hhit,center,rmed,age);
// loop while step is large
    unsigned niter(0);
    Hep3Vector descent(1.0,0.0,0.0);
    while(lambda*descent.mag() > _minlambda && niter < _maxniter){
// fill the sums for computing the descent vector
      AGESums sums;
      fillSums(hhit,center,rmed,sums);
// descent vector cases: if the inner vs outer difference is significant (compared to the median), damp using the median sums,
// otherwise not.  These expressions take care of the undiferentiable condition on the boundary. 
      double dx(sums._sco-sums._sci);
      double dy(sums._sso-sums._ssi);
      if(fabs(dx) < sums._scc)
        dx += (sums._sco < sums._sci) ? -sums._scc : sums._scc;
      if(fabs(dy) < sums._ssc)
        dy += (sums._sso < sums._ssi) ? -sums._ssc : sums._ssc;
      descent = Hep3Vector(dx,dy,0.0);
// compute error function, decreasing lambda until this is better than the previous
      double agenew;
      Hep3Vector cnew = center + lambda*descent;
      findAGE(hhit,cnew,rmed,agenew);
// if we've improved, increase the step size and iterate
      if(agenew < age){
        lambda *= (1.0+_lstep);
      } else {
// if we haven't improved, keep reducing the step till we do
        unsigned miter(0);
        while(agenew > age && miter < _maxniter && lambda*descent.mag() > _minlambda){
          lambda *= (1.0-_lstep);
          cnew = center + lambda*descent;
          findAGE(hhit,cnew,rmed,agenew);
          ++miter;
        }
// if this fails, reverse the descent drection and try again
        if(agenew > age){
	  descent *= -1.0;
	  lambda *= (1.0 +_lstep);
          cnew = center + lambda*descent;
          findAGE(hhit,cnew,rmed,agenew);
        }
      }
// prepare for next iteration
      if(agenew < age){
        center = cnew;
        age = agenew;
      } else {
	static const double minage(0.1);
	if(_debug > 0 && agenew-age>minage)
	  std::cout << "iteration did not improve AGE!!! lambda = " 
	  << lambda  << " age = " << age << " agenew = " << agenew << std::endl;
	break;
      }
// if we're constraining to intersect the target, adjust the center and radius if necessary
      if(_targetinter) {
	double rperigee = center.perp()-rmed;
	if(fabs(rperigee) > _targetradius){
// adjust both center position and radius till they touch the target, holding phi constant.  This keeps the circle near the hits.  Sign matters!
	  double dr;
	  if(rperigee > 0)
	    dr = 0.5*(rperigee - _targetradius);  // change is 1/2 the difference
	  else
	    dr = 0.5*(rperigee + _targetradius);
	  rmed += dr;
	  // direction radially outwards from origin to the center
	  Hep3Vector pdir = Hep3Vector(center.x(),center.y(),0.0).unit();
	  // the center moves opposite the radius
	  center -= dr*pdir;
	}
      }
      ++niter;
    }
// check for convergence
    if(_debug > 0 && niter > _maxniter ){
      std::cout << "AGE didn't converge!!! " << std::endl;
//      return false;
    }
    return true;
  }
  
  bool
  RobustHelixFit::findZ(HelixHitCollection& hhits,RobustHelix& myhel) {
    using namespace boost::accumulators;
// sort points by z
    std::sort(hhits.begin(),hhits.end(),zcomp());
// find phi information
    std::vector<FZ> finfo;
    finfo.reserve(hhits.size());
    for(auto hhit : hhits) {
      if(use(hhit) && (stereo(hhit) || (!_stereoinit)) ){
	FZ fz;
	phiInfo(myhel.center(),hhit,fz._phi);
	fz._z = hhit._pos.z();
	finfo.push_back(fz);
      }
    }
// make initial estimate of dfdz using 'nearby' pairs
    accumulator_set<double, stats<tag::median(with_p_square_quantile)> > accf;
    for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
      for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	double dz = finfo[jphi]._z - finfo[iphi]._z;
	if(dz > _minzsep && dz < _maxzsep){
	  double dphi = deltaPhi(finfo[iphi]._phi._val,finfo[jphi]._phi._val);
	  if(fabs(dphi) > _mindphi && fabs(dphi) < _maxdphi){
	    double slope = dphi/dz;
	    if(slope > _smin && slope < _smax){ 
	      accf(slope);
	    }
	  }
	}
      }
    }
    if(boost::accumulators::extract::count(accf) > _minnhit){
      double  dfdz = extract_result<tag::median>(accf);
      //      double dfdztest = extract_result<tag::mean>(acctest);
      // if the sign of dfdz disagrees, abort
      if( dfdz * _helicity._value < 0.0 || dfdz > _smax || dfdz < _smin) return false;
// find phi at z intercept.  Use a histogram technique since phi looping
// hasn't been resolved yet
      TH1F hphi("hphi","phi value",50,-1.1*CLHEP::pi,1.1*CLHEP::pi);
      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
      	double phiex = finfo[iphi]._z*dfdz;
	double dphi = deltaPhi(phiex,finfo[iphi]._phi._val);
	hphi.Fill(dphi);
	hphi.Fill(dphi-CLHEP::twopi);
	hphi.Fill(dphi+CLHEP::twopi);
      }
      double fz0 = hphi.GetBinCenter(hphi.GetMaximumBin());
      if(fz0>CLHEP::pi)fz0 -= CLHEP::twopi;
      if(fz0<-CLHEP::pi)fz0 += CLHEP::twopi;
// reset for full search
      if(_stereoinit){
	finfo.clear();
	for(auto hhit : hhits) {
	  if(use(hhit)){
	    FZ fz;
	    phiInfo(myhel.center(),hhit,fz._phi);
	    fz._z = hhit._pos.z();
	    finfo.push_back(fz);
	  }
	}
      }
  // iterate over slope and ambiguity resolution
      bool changed(true);
      unsigned niter(0);
      while(changed && niter < _maxniter){
  // resolve phi over the full z range using the initial slope and intercept
	changed = false;
	++niter;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  double phiex = fz0 + finfo[iphi]._z*dfdz;
	  int nloop = (int)rint((phiex-finfo[iphi]._phi._val)/CLHEP::twopi);
	  finfo[iphi]._phi._val += nloop*CLHEP::twopi;
	  changed |= nloop != 0;
	}
	// make a long-range estimate of slope
	accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accf2;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	    double dz = finfo[jphi]._z -finfo[iphi]._z;
	    double dphi = finfo[jphi]._phi._val-finfo[iphi]._phi._val;
	    if(dz > _minzsep && fabs(dphi) > _mindphi ){
	      double dphiex = dz*dfdz;
	      if(!_filterz || fabs(dphi-dphiex) < _pout){
		double slope = dphi/dz;
		if(slope > _smin && slope < _smax){ 
		  accf2(slope);
		}
	      }
	    }
	  }
	}
	dfdz = extract_result<tag::median>(accf2);
	accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acci2;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  double phiex = fz0+finfo[iphi]._z*dfdz;
	  if(!_filterz || fabs(finfo[iphi]._phi._val-phiex) < _pout ){
	    acci2(finfo[iphi]._phi._val - finfo[iphi]._z*dfdz);
	  }
	}
	fz0 = fmod(extract_result<tag::median>(acci2),CLHEP::twopi);
	if(fz0>CLHEP::pi)fz0 -= CLHEP::twopi;
	if(fz0<-CLHEP::pi)fz0 += CLHEP::twopi;
      }
      // inverse convention: should solve directly for lambda, FIXME!
      myhel.lambda() = 1.0/dfdz;
      myhel.fz0() = fz0;
      // fix the phi for the hit points
      for(auto hhit : hhits) {
	double phiex = myhel.fz0() + hhit._pos.z()/myhel.lambda();
	FZ fz;
	phiInfo(myhel.center(),hhit,fz._phi);
	int nloop = (int)rint((phiex - fz._phi._val)/CLHEP::twopi);
	hhit._phi = fz._phi._val + nloop*CLHEP::twopi;
	if(_filterz && fabs(hhit._phi-phiex)> _pout) setOutlier(hhit);
      }
      return true;
    } else
      return false;
  }


  bool
  RobustHelixFit::initCircle(HelixHitCollection const& hhits,RobustHelix& myhel) {
    bool retval(false);
    static const double mind2 = _mindist*_mindist;
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx, accy, accr;
    // form all triples, and compute the circle center for unaligned hits.  I can aford to be choosy
    unsigned ntriple(0);
    std::vector<CLHEP::Hep3Vector> pos;
    pos.reserve(hhits.size());
    for(auto hhit : hhits) {
      if(use(hhit) && (stereo(hhit) || (!_stereoinit)) ){
	pos.push_back(hhit._pos);
      }
    }
    // loop over all triples
    size_t np = pos.size();
    for(size_t ip=0;ip<np;++ip){
      // pre-compute some values
      double ri2 = pos[ip].perp2();
      for(size_t jp=ip+1;jp<np; ++jp){
	if(pos[ip].perpPart().diff2(pos[jp].perpPart()) > mind2){
	  double rj2 = pos[jp].perp2();
	  size_t mink = jp+1;
	  for(size_t kp=mink;kp<np; ++kp){
	    if(pos[ip].perpPart().diff2(pos[kp].perpPart()) > mind2 &&
		pos[jp].perpPart().diff2(pos[kp].perpPart()) > mind2){
	      // this effectively measures the slope difference
	      double delta = (pos[kp].x() - pos[jp].x())*(pos[jp].y() - pos[ip].y()) - 
		(pos[jp].x() - pos[ip].x())*(pos[kp].y() - pos[jp].y());
	      if(fabs(delta) > _mindelta){
		double rk2 = pos[kp].perp2();
		// find circle center for this triple
		double cx = 0.5* (
		    (pos[kp].y() - pos[jp].y())*ri2 + 
		    (pos[ip].y() - pos[kp].y())*rj2 + 
		    (pos[jp].y() - pos[ip].y())*rk2 ) / delta;
		double cy = -0.5* (
		    (pos[kp].x() - pos[jp].x())*ri2 + 
		    (pos[ip].x() - pos[kp].x())*rj2 + 
		    (pos[jp].x() - pos[ip].x())*rk2 ) / delta;
		double rho = sqrt(pow(pos[ip].x()-cx,2)+pow(pos[ip].y()-cy,2));
		double rc = sqrt(cx*cx + cy*cy);
		double rmin = fabs(rc-rho);
		double rmax = rc+rho;
		// test circle parameters for this triple: should be inside the tracker,
		// optionally consistent with the target
		if(rho > _rmin && rho< _rmax && rmax < _trackerradius
		  && ( (!_targetinit) || rmin < _targetradius) ) {
		  // accumulate 
		  ++ntriple;
		  accx(cx);
		  accy(cy);
		  accr(rho);
		}
	      }
	    }
	  }
	}
      }
    }
    if(ntriple > _minnhit){
      retval = true;
      double centx = extract_result<tag::median>(accx);
      double centy = extract_result<tag::median>(accy);
      double rho = extract_result<tag::median>(accr);
      myhel.center() = CLHEP::Hep3Vector(centx,centy,0.0);
      myhel.radius() = rho;
    }
    return retval;
  }

  void RobustHelixFit::findAGE(HelixHitCollection const& hhits, Hep3Vector const& center,double& rmed, double& age) {
    using namespace boost::accumulators;
    // protection against empty data
    if(hhits.size() == 0)return;
    // fill radial information for all points, given this center
    std::vector<VALERR> radii;
    for(auto hhit : hhits) {
      if(use(hhit)){
	// find radial information for this point
	VALERR rad;
	radInfo(center,hhit,rad);
	radii.push_back(rad);
	// compute the normalization too
      }
    }
    // find the median radius
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accr;
    for(unsigned irad=0;irad<radii.size();++irad){
      accr(radii[irad]._val); 
    }
    rmed = extract_result<tag::median>(accr);
    // now compute the AGE (Absolute Geometric Error)
    age = 0.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      age += fabs(radii[irad]._val-rmed);
    }
  }

  void
  RobustHelixFit::fillSums(HelixHitCollection const& hhits, Hep3Vector const& center,double rmed,AGESums& sums) {
    // initialize sums
    sums.clear();
    // compute the transverse sums
    for(auto hhit : hhits) {
      if(use(hhit)){
	// find radial information for this point
	VALERR rad;
	radInfo(center,hhit,rad);
	// now x,y projections
	double pcos = (hhit._pos.x()-center.x())/rad._val;
	double psin = (hhit._pos.y()-center.y())/rad._val;
	// 3 conditions: either the radius is inside the median, outside the median, or 'on' the median.  We define 'on'
	// in terms of a window
	if(fabs(rmed - rad._val) < _rwind  ){
	  sums._scc += fabs(pcos);
	  sums._ssc += fabs(psin);
	  ++sums._nc;
	} else if (rad._val > rmed  ) {
	  sums._sco += pcos;
	  sums._sso += psin;
	  ++sums._no;
	} else {
	  sums._sci += pcos;
	  sums._ssi += psin;        
	  ++sums._ni;
	}
      }  
    }
  }

  void
  RobustHelixFit::filterDist(HelixHitCollection& hhits) {
    using namespace boost::accumulators;
    // first, resolve phi.   Use the average X and Y to define the initial
    // phi value, to avoid looping issues
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
    for(auto hhit : hhits ) {
      if(use(hhit)){
	accx(hhit._pos.x());
	accy(hhit._pos.y());
      }
    }
    double mx = extract_result<tag::median>(accx);
    double my = extract_result<tag::median>(accy);
    double mphi = atan2(my,mx);
    for(auto hhit : hhits) {
      if(use(hhit)){
	double dphi = hhit._phi - mphi;
	if(fabs(dphi) > CLHEP::pi){
	  if(dphi > 0)
	    hhit._phi -= CLHEP::twopi;
	  else
	    hhit._phi += CLHEP::twopi;
	}
	setResolvedPhi(hhit);
      }
    }
    // now cut
    CLHEP::Hep3Vector mh(mx,my,0.0);
    for(auto hhit : hhits) {
      double dist = sqrt(hhit._pos.perpPart().diff2(mh));
      if(dist > _maxdist)setOutlier(hhit);
    }
  }

  void
  RobustHelixFit::filterXY(HelixHitCollection& hhits, Hep3Vector const& center,double rmed,bool& changed) {
    static StrawHitFlag other(StrawHitFlag::other);
    changed = false;
    for(auto hhit : hhits) {
      VALERR rad;
      radInfo(center,hhit,rad);
      bool olduse = use(hhit);
  // update the flag each iteration, until we've converged
      if(fabs(rad._val -rmed) > _rout)
	hhit._flag.merge(other);
      else
	hhit._flag.clear(other);
      changed |= olduse != use(hhit);
    }
  }

  double
  RobustHelixFit::deltaPhi(double phi1, double phi2){
    double dphi = fmod(phi2-phi1,CLHEP::twopi);
    if(dphi>CLHEP::pi)dphi -= CLHEP::twopi;
    if(dphi<-CLHEP::pi)dphi += CLHEP::twopi;
    return dphi;
  }

  bool
  RobustHelixFit::use(HelixHit const& hhit) const {
     return (!hhit._flag.hasAnyProperty(_dontuseflag))
      && (hhit._flag.hasAllProperties(_useflag) || _useflag.empty());
  }

  bool
  RobustHelixFit::stereo(HelixHit const& hhit) const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return hhit._flag.hasAllProperties(stereo);
  }

  void RobustHelixFit::setResolvedPhi(HelixHit& hhit) const {
    static StrawHitFlag resphi(StrawHitFlag::resolvedphi);
    hhit._flag.merge(resphi);
  }

  void RobustHelixFit::setOutlier(HelixHit& hhit) const {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    hhit._flag.merge(outlier);
  }

 void
  RobustHelixFit::radInfo(Hep3Vector const& center,HelixHit const& hhit, VALERR& rad) const {
    double rvec = Hep3Vector(hhit._pos - center).perp();
    rad._val = rvec;
    rad._err = hhit._wres; // should take the projection FIXME!
    if(_debug > 3)std::cout << "radInfo : r = " << rad._val << " rerr = " << rad._err  << std::endl;
  }

  void
  RobustHelixFit::phiInfo(Hep3Vector const& center, HelixHit const& hhit, VALERR& phi) const {
    CLHEP::Hep3Vector rad = hhit._pos - center;
    double phi0 = rad.phi();
    phi._val = phi0;
    phi._err = hhit._wres/rad.mag(); // should take the projection FIXME!
    if(_debug > 3)std::cout << "phiInfo : phi = " << phi._val << " ferr = " << phi._err << std::endl;
  }

  unsigned RobustHelixFit::hitCount(HelixHitCollection const& hhits) const {
    unsigned retval(0);
    for(auto hhit : hhits)
      if(use(hhit))++retval;
    return retval;
  }
 
}

