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
#include <cmath>
using CLHEP::Hep3Vector;
using namespace std;
using namespace boost::accumulators;
namespace mu2e 
{
  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<HelixHit,HelixHit,bool> {
    bool operator()(HelixHit const& p1, HelixHit const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("UseFlag",vector<string>{"Outlier"})),
    _mindelta(pset.get<double>("minDelta",5000.0)),
    _minnhit(pset.get<unsigned>("minNHit",10)),
    _maxphisep(pset.get<double>("MaxPhiHitSeparation",1.8)),
    _lambda0(pset.get<double>("lambda0",1.0)),
    _lstep(pset.get<double>("lstep",0.2)),
    _minlambda(pset.get<double>("minlambda",0.01)),
    _nphibins(pset.get<unsigned>("NPhiHistBins",30)),
    _minnphi(pset.get<unsigned>("MinNPhi",3)),
    _maxniter(pset.get<unsigned>("maxniter",50)),
    _minzsep(pset.get<double>("minzsep",100.0)),
    _maxzsep(pset.get<double>("maxzsep",700.0)),
    _mindphi(pset.get<double>("mindphi",0.25)),
    _maxdphi(pset.get<double>("maxdphi",2.5)),
    _mindist(pset.get<double>("mindist",50.0)),
    _maxdist(pset.get<double>("maxdist",500.0)),
    _rmin(pset.get<double>("minR",200.0)),
    _rmax(pset.get<double>("maxR",500.0)),
    _lmin(pset.get<double>("minAbsLambda",50.0)),
    _lmax(pset.get<double>("maxAbsLambda",500.0)),
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
    if(_helicity._value < 0){
      _lmin *= -1.0;
      _lmax *= -1.0;
    }
  }

  RobustHelixFit::~RobustHelixFit()
  {}

  void RobustHelixFit::fitHelix(HelixSeed& hseed) {
    // reset the fit status flags, in case this is called iteratively
    hseed._status.clear(TrkFitFlag::circleOK);
    hseed._status.clear(TrkFitFlag::phizOK);
    hseed._status.clear(TrkFitFlag::helixOK);
    // gross outlier removal on init
    if(!hseed._status.hasAllProperties(TrkFitFlag::initOK))
      filterSector(hseed._hhits);
    // count what's left
    if(hitCount(hseed._hhits) >= _minnhit){
      hseed._status.merge(TrkFitFlag::hitsOK);
      // solve for the circle parameters
      fitXY(hseed);
      if(hseed._status.hasAnyProperty(TrkFitFlag::circleOK)){
	// solve for the longitudinal parameters
	fitFZ(hseed);
	if(hseed._status.hasAnyProperty(TrkFitFlag::phizOK)){
	  // final test
	  if (goodHelix(hseed._helix)){
	    hseed._status.merge(TrkFitFlag::helixOK);
	  }
	}
      }
    }
  }

  void RobustHelixFit::fitXY(HelixSeed& hseed) {
    HelixHitCollection& hhits = hseed._hhits;
    RobustHelix& rhel = hseed._helix;
    // if needed, initialize
    if(!hseed._status.hasAllProperties(TrkFitFlag::initOK)){
      if(!initXY(hhits,rhel))
	return;
    }
// this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
// setup working variables
    double age;
    Hep3Vector center = rhel.center();
    double rmed = rhel.radius();
// initialize step
    double lambda = _lambda0;
// find median and AGE for the initial center
    findAGE(hhits,center,rmed,age);
// loop while step is large
    unsigned niter(0);
    Hep3Vector descent(1.0,0.0,0.0);
    while(lambda*descent.mag() > _minlambda && niter < _maxniter){
// fill the sums for computing the descent vector
      AGESums sums;
      fillSums(hhits,center,rmed,sums);
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
      findAGE(hhits,cnew,rmed,agenew);
// if we've improved, increase the step size and iterate
      if(agenew < age){
        lambda *= (1.0+_lstep);
      } else {
// if we haven't improved, keep reducing the step till we do
        unsigned miter(0);
        while(agenew > age && miter < _maxniter && lambda*descent.mag() > _minlambda){
          lambda *= (1.0-_lstep);
          cnew = center + lambda*descent;
          findAGE(hhits,cnew,rmed,agenew);
          ++miter;
        }
// if this fails, reverse the descent drection and try again
        if(agenew > age){
	  descent *= -1.0;
	  lambda *= (1.0 +_lstep);
          cnew = center + lambda*descent;
          findAGE(hhits,cnew,rmed,agenew);
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
	forceTargetInter(center,rmed);
      }
      ++niter;
    }
// check for convergence
    if(_debug > 0 && niter > _maxniter ){
      std::cout << "AGE didn't converge!!! " << std::endl;
    }
    // update parameters
    rhel.center() = center;
    rhel.radius() = rmed;
    // update flag
    if(goodCircle(rhel))
      hseed._status.merge(TrkFitFlag::circleOK);
    if(niter < _maxniter)
      hseed._status.merge(TrkFitFlag::circleConverged);
  }

  void RobustHelixFit::forceTargetInter(Hep3Vector& center, double& radius) {
  // find the smallest radius
    double rperigee = center.perp()-radius;
    // if this is outside the target radius, adjust the parameters
    if(fabs(rperigee) > _targetradius){
      // adjust both center position and radius till they touch the target, holding phi constant.  This keeps the circle near the hits.  Sign matters!
      double dr;
      if(rperigee > 0)
	dr = 0.5*(rperigee - _targetradius);  // change is 1/2 the difference
      else
	dr = 0.5*(rperigee + _targetradius);
      radius += dr;
      // direction radially outwards from origin to the center
      Hep3Vector pdir = Hep3Vector(center.x(),center.y(),0.0).unit();
      // the center moves opposite the radius
      center -= dr*pdir;
    }
  }

  bool RobustHelixFit::initFZ(HelixHitCollection& hhits,RobustHelix& rhel) {
    bool retval(false);
    // initialize z parameters
    rhel.lambda() = 1.0e12; //infinite slope for now
    rhel.fz0() = 0.0;
    static TrkFitFlag circleOK(TrkFitFlag::circleOK);
    static TrkFitFlag helixOK(TrkFitFlag::helixOK);
    // sort points by z
    std::sort(hhits.begin(),hhits.end(),zcomp());
    // initialize phi for this center (without resolving loop);
    for(auto& hhit : hhits){
      initPhi(hhit,rhel);
    }
    // make initial estimate of dfdz using 'nearby' pairs.  This insures they are
    // on the same loop
    accumulator_set<double, stats<tag::median(with_p_square_quantile)> > accf;
    for(auto ihit=hhits.begin(); ihit != std::prev(hhits.end()); ++ihit) {
      if(use(*ihit) && ( (!_stereoinit) || stereo(*ihit))) {
	for(auto jhit = std::next(ihit); jhit != hhits.end(); ++jhit) {
	  if(use(*jhit) && ( (!_stereoinit) || stereo(*jhit))) {
	    double dz = jhit->_pos.z() - ihit->_pos.z();
	    if(dz > _minzsep && dz < _maxzsep){
	      double dphi = fabs(deltaPhi(ihit->_phi,jhit->_phi));
	      if(dphi > _mindphi && dphi < _maxdphi){
		// compute the slope for this pair of hits
		double lambda = dz/dphi;
		if(lambda > _lmin && lambda < _lmax){ 
		  accf(lambda);
		} // lambda in range
	      } // good phi separation
	    } // good z separation
	  } // good 2nd hit 
	} // loop over 2nd hit
      } // good first hit
    } // loop over first hit
    // extract the lambda
    if(boost::accumulators::extract::count(accf) > _minnhit){
      double lambda = extract_result<tag::median>(accf);
      // check the range 
      if( lambda < _lmax && lambda > _lmin) {
	// update helix
	rhel.lambda() = lambda;
	// find phi at z intercept.  Use a histogram technique since phi looping
	// hasn't been resolved yet
	TH1F hphi("hphi","phi value",_nphibins,-1.1*CLHEP::pi,1.1*CLHEP::pi);
	for(auto const& hhit : hhits) {
	  if(use(hhit) && ( (!_stereoinit) || stereo(hhit))) {
	    double phiex = rhel.circleAzimuth(hhit._pos.z());
	    double dphi = deltaPhi(phiex,hhit._phi);
	    hphi.Fill(dphi);
	    // make sure to cover the edges
	    hphi.Fill(dphi-CLHEP::twopi);
	    hphi.Fill(dphi+CLHEP::twopi);
	  }
	}
	// extract the intercept as the maximum bin's center
	unsigned nmax = (unsigned)(ceil(hphi.GetMaximum()));
	double fz0 = hphi.GetBinCenter(hphi.GetMaximumBin());
	if(nmax > _minnphi) {
	  // choose the intercept to have |fz0| < pi.  This is purely a convention
	  rhel.fz0() = deltaPhi(0.0,fz0);
	  // update the hits to resolved phi looping with these parameters
	  for(auto& hhit : hhits) {
	    resolvePhi(hhit,rhel);
	  }
	  retval = true;
	} // good # for fz0
      } // good slope
    } // good # of hits
    return retval;
  }

  void RobustHelixFit::fitFZ(HelixSeed& hseed) {
    HelixHitCollection& hhits = hseed._hhits;
    RobustHelix& rhel = hseed._helix;
    // if required, initialize
    if(!hseed._status.hasAllProperties(TrkFitFlag::initOK)){
      if(initFZ(hhits,rhel))
	hseed._status.merge(TrkFitFlag::initOK);
      else
	return;
    }
  // iterate over lambda and loop resolution
    bool changed(true);
    unsigned niter(0);
    while(changed && niter < _maxniter){
      changed = false;
      ++niter;
      // make a long-range estimate of slope
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accf2;
      for(auto ihit=hhits.begin(); ihit != std::prev(hhits.end()); ++ihit) {
	if(use(*ihit)) {
	  for(auto jhit = std::next(ihit); jhit != hhits.end(); ++jhit) {
	    if(use(*jhit)) {
	      double dz = jhit->_pos.z() - ihit->_pos.z();
	      double dphi = jhit->_phi-ihit->_phi; // includes loop determination
	      if(dz > _minzsep && fabs(dphi) > _mindphi){
		double lambda = dz/dphi;
		if(lambda > _lmin && lambda < _lmax){ 
		  accf2(lambda);
		} // good slope
	      } // minimum hit separation
	    } // good 2nd hit 
	  } // loop on 2nd hit
	} // good 1st hit
      } // loop on 1st hit
      // extract slope
      rhel.lambda() = extract_result<tag::median>(accf2);
//  test parameters for return value
      // now extract intercept.  Here we solve for the difference WRT the previous value
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > acci2;
      for(auto const& hhit : hhits) {
	if(use(hhit)) {
	  double phiex = rhel.circleAzimuth(hhit._pos.z());
	  double dphi = deltaPhi(phiex,hhit._phi);
	  acci2(dphi);// accumulate the difference WRT the current intercept
	}
      }
      double dphi = extract_result<tag::median>(acci2);
      rhel.fz0() += dphi;
      // resolve the hit loops again
      for(auto& hhit : hhits)
	changed |= resolvePhi(hhit,rhel);
    }
    // set the flags
    if(goodFZ(rhel))
      hseed._status.merge(TrkFitFlag::phizOK);
    if(niter < _maxniter)
      hseed._status.merge(TrkFitFlag::phizConverged);
  }

  bool RobustHelixFit::initXY(HelixHitCollection const& hhits,RobustHelix& rhel) {
    bool retval(false);
    static const double mind2 = _mindist*_mindist;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx, accy, accr;
    // pick out a subset of hits.  I can aford to be choosy
    unsigned ntriple(0);
    std::vector<CLHEP::Hep3Vector> pos;
    pos.reserve(hhits.size());
    for(auto const& hhit : hhits) {
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
		double rho = sqrt(std::pow(pos[ip].x()-cx,(int)2)+std::pow(pos[ip].y()-cy,(int)2));
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
		} // radius meets requirements
	      } // slope is not the same
	    } // points 1 and 3 and 2 and 3 are separated
	  } // 3rd point
	} // points 1 and 2 are separated
      } // second point
    } // first point
    // median calculation needs a reasonable number of points to function
    if(ntriple > _minnhit){
      retval = true;
      double centx = extract_result<tag::median>(accx);
      double centy = extract_result<tag::median>(accy);
      double rho = extract_result<tag::median>(accr);
      rhel.center() = CLHEP::Hep3Vector(centx,centy,0.0);
      rhel.radius() = rho;
    }
    return retval;
  }

  void RobustHelixFit::findAGE(HelixHitCollection const& hhits, Hep3Vector const& center,double& rmed, double& age) {
  // fill radial information for all points, given this center
    std::vector<double> radii;
    for(auto const& hhit : hhits) {
      if(use(hhit)){
	// find radial information for this point
	double rad = Hep3Vector(hhit._pos - center).perp();
	radii.push_back(rad);
      }
    }
    if(radii.size() > _minnhit){
      // find the median radius
      accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accr;
      for(unsigned irad=0;irad<radii.size();++irad){
	accr(radii[irad]); 
      }
      rmed = extract_result<tag::median>(accr);
      // now compute the AGE (Absolute Geometric Error)
      age = 0.0;
      for(unsigned irad=0;irad<radii.size();++irad){
	age += fabs(radii[irad]-rmed);
      }
    }
  }

  void RobustHelixFit::fillSums(HelixHitCollection const& hhits, Hep3Vector const& center,double rmed,AGESums& sums) {
    // initialize sums
    sums.clear();
    // compute the transverse sums
    for(auto const& hhit : hhits) {
      if(use(hhit)){
	// find radial information for this point
	double rad = Hep3Vector(hhit._pos - center).perp();
	// now x,y projections
	double pcos = (hhit._pos.x()-center.x())/rad;
	double psin = (hhit._pos.y()-center.y())/rad;
	// 3 conditions: either the radius is inside the median, outside the median, or 'on' the median.  We define 'on'
	// in terms of a window
	if(fabs(rmed - rad) < _rwind  ){
	  sums._scc += fabs(pcos);
	  sums._ssc += fabs(psin);
	  ++sums._nc;
	} else if (rad > rmed  ) {
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

  double RobustHelixFit::deltaPhi(double phi1, double phi2){
    double dphi = fmod(phi2-phi1,CLHEP::twopi);
    if(dphi>CLHEP::pi)dphi -= CLHEP::twopi;
    if(dphi<-CLHEP::pi)dphi += CLHEP::twopi;
    return dphi;
  }

  bool RobustHelixFit::use(HelixHit const& hhit) const {
     return (!hhit._flag.hasAnyProperty(_dontuseflag))
      && (hhit._flag.hasAllProperties(_useflag) || _useflag.empty());
  }

  bool RobustHelixFit::stereo(HelixHit const& hhit) const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return hhit._flag.hasAllProperties(stereo);
  }

  void RobustHelixFit::setOutlier(HelixHit& hhit) const {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    hhit._flag.merge(outlier);
  }

  unsigned RobustHelixFit::hitCount(HelixHitCollection const& hhits) const {
    unsigned retval(0);
    for(auto hhit : hhits)
      if(use(hhit))++retval;
    return retval;
  }

  void RobustHelixFit::initPhi(HelixHit& hhit, RobustHelix const& rhel) const {
  // ray from the circle center to the point
    hhit._phi = Hep3Vector(hhit._pos - rhel.center()).phi();
  }

  void RobustHelixFit::filterSector(HelixHitCollection& hhits) {
    // Use the average X and Y to define the average phi of the hits
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
    for(auto const& hhit : hhits ) {
      if(use(hhit)){
	accx(hhit._pos.x());
	accy(hhit._pos.y());
      }
    }
    double mx = extract_result<tag::median>(accx);
    double my = extract_result<tag::median>(accy);
    double mphi = atan2(my,mx);
    // require the hits be within a maximum of this global azimuth
    for(auto& hhit : hhits) {
      if(use(hhit)){
	double dphi = fabs(deltaPhi(hhit.pos().phi(),mphi));
	if(dphi > _maxphisep)setOutlier(hhit);
      }
    }
  }

  bool RobustHelixFit::resolvePhi(HelixHit& hhit, RobustHelix const& rhel) const {
    // find phi expected
    double phiex = rhel.circleAzimuth(hhit._pos.z());
    int nloop = (int)rint((phiex-hhit._phi)/CLHEP::twopi);
    hhit._phi += nloop*CLHEP::twopi;
    // now update the flag
    static StrawHitFlag resphi(StrawHitFlag::resolvedphi);
    hhit._flag.merge(resphi);
    return nloop != 0;
  }

  bool RobustHelixFit::goodFZ(RobustHelix const& rhel) {
    return rhel.lambda() > _lmin && rhel.lambda() < _lmax;
  }

  bool RobustHelixFit::goodCircle(RobustHelix const& rhel) {
    return rhel.radius() > _rmin && rhel.radius() < _rmax;
  }

  bool RobustHelixFit::goodHelix(RobustHelix const& rhel) {
    return goodCircle(rhel) && goodFZ(rhel);
  }
}

