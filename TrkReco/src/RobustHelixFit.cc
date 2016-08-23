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
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
// root
#include "TH1F.h"
// C++
#include <vector>
#include <string>
using namespace CLHEP;
namespace mu2e 
{
// comparison functor for ordering points
  struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
    bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
  };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<XYZP,XYZP,bool> {
    bool operator()(XYZP const& p1, XYZP const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  
  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _mindelta(pset.get<double>("minDelta",5000.0)),
    _minnhit(pset.get<unsigned>("minNHit",10)),
    _lambda0(pset.get<double>("lambda0",1.0)),
    _lstep(pset.get<double>("lstep",0.2)),
    _minlambda(pset.get<double>("minlambda",0.01)),
    _maxniter(pset.get<unsigned>("maxniter",50)),
    _nsigma(pset.get<double>("nsigma",5.0)),
    _minzsep(pset.get<double>("minzsep",100.0)),
    _maxzsep(pset.get<double>("maxzsep",700.0)),
    _mindphi(pset.get<double>("mindphi",0.25)),
    _maxdphi(pset.get<double>("maxdphi",2.5)),
    _efac(pset.get<double>("ErrorFactor",1.0)),
    _mindist(pset.get<double>("mindist",50.0)),
    _maxdist(pset.get<double>("maxdist",500.0)),
    _rmin(pset.get<double>("minR",200.0)),
    _rmax(pset.get<double>("maxR",500.0)),
    _tdmin(pset.get<double>("minAbsTanDip",0.3)),
    _tdmax(pset.get<double>("maxAbsTanDip",2.0)),
    _force(pset.get<bool>("forceP",false)),
    _xyweights(pset.get<bool>("xyWeights",false)),
    _zweights(pset.get<bool>("zWeights",false)),
    _filterxy(pset.get<bool>("filterxy",true)),
    _filterz(pset.get<bool>("filterz",true)),
    _stereoinit(pset.get<bool>("stereoinit",false)),
    _stereofit(pset.get<bool>("stereofit",false)),
    _targetinit(pset.get<bool>("targetinit",true)),
    _targetinter(pset.get<bool>("targetintersect",true)),
    _targetradius(pset.get<double>("targetradius",75.0)),
    _trackerradius(pset.get<double>("trackerradius",700.0)),
    _helicity((Helicity::helicity)pset.get<int>("Helicity",0))
  {
    XYZP::_efac = _efac;
    std::vector<std::string> bitnames;
    // these should be fcl parameters, FIXME!!!
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    XYZP::_dontuseflag = StrawHitFlag(bitnames);
    if(_stereofit)XYZP::_useflag = StrawHitFlag(StrawHitFlag::stereo);
    XYZP::_debug = _debug;
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
  RobustHelixFit::findHelix(StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,HelixSeed& myseed) {
    RobustHelix& myhel = myseed._helix;
// Convert the hit positions to what we need for helix finding
    XYZPVector  xyzp;
    fillXYZP(shcol,shpcol,myseed._timeCluster._strawHitIdxs,xyzp);
// filter by geometry
    if(_filterxy)filterDist(xyzp);
    if(xyzp.size() >= _minnhit){
      myseed._status.merge(TrkFitFlag::hitsOK);
      // initialize the circle parameters if required.
      if(myseed._status.hasAllProperties(TrkFitFlag::initialized) ||
	  initCircle(xyzp,myhel)){
	myseed._status.merge(TrkFitFlag::initialized);
	// solve for the circle parameters
	if(findXY(xyzp,myhel)){
	  myseed._status.merge(TrkFitFlag::radiusOK);
	  if(findZ(xyzp,myhel)) {
	    // set the success
	    myseed._status.merge(TrkFitFlag::phizOK);
	    myseed._status.merge(TrkFitFlag::fitOK);
	  }
	}
      }
    }
  }

  bool
    RobustHelixFit::findXY(XYZPVector& xyzp,RobustHelix& myhel) {
      double rmed, age;
      Hep3Vector center = myhel.center();
      // then, refine that using weights
      bool changed(true);
    unsigned niter(0);
    while(niter < _maxniter && changed){
      findCenterAGE(xyzp,center,rmed,age,_xyweights);
      if(_filterxy)
	filterXY(xyzp,center,rmed,changed);
      else
	changed = false;
      niter++;
    }
    myhel.center() = center;
    myhel.radius() = rmed;
    return true;
  }

  bool
  RobustHelixFit::findCenterAGE(XYZPVector const& xyzp,Hep3Vector& center, double& rmed, double& age,bool useweights) {
// this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
// initialize step
    double lambda = _lambda0;
// find median and AGE for the initial center
    findAGE(xyzp,center,rmed,age,useweights);
// loop while step is large
    unsigned niter(0);
    Hep3Vector descent(1.0,0.0,0.0);
    while(lambda*descent.mag() > _minlambda && niter < _maxniter){
// fill the sums for computing the descent vector
      SUMS sums;
      fillSums(xyzp,center,rmed,sums,useweights);
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
      findAGE(xyzp,cnew,rmed,agenew,useweights);
// if we've improved, increase the step size and iterate
      if(agenew < age){
        lambda *= (1.0+_lstep);
      } else {
// if we haven't improved, keep reducing the step till we do
        unsigned miter(0);
        while(agenew > age && miter < _maxniter && lambda*descent.mag() > _minlambda){
          lambda *= (1.0-_lstep);
          cnew = center + lambda*descent;
          findAGE(xyzp,cnew,rmed,agenew,useweights);
          ++miter;
        }
// if this fails, reverse the descent drection and try again
        if(agenew > age){
	  descent *= -1.0;
	  lambda *= (1.0 +_lstep);
          cnew = center + lambda*descent;
          findAGE(xyzp,cnew,rmed,agenew,useweights);
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
  RobustHelixFit::findZ(XYZPVector& xyzp,RobustHelix& myhel) {
    using namespace boost::accumulators;
// sort points by z
    std::sort(xyzp.begin(),xyzp.end(),zcomp());
// find phi information
    std::vector<FZ> finfo;
    finfo.reserve(xyzp.size());
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use() && (xyzp[ixyzp].stereo() || (!_stereoinit)) ){
	FZ fz;
	xyzp[ixyzp].finfo(myhel.center(),fz._phi);
	fz._z = xyzp[ixyzp]._pos.z();
	finfo.push_back(fz);
      }
    }
// make initial estimate of dfdz using 'nearby' pairs
    accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf;
    for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
      for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	double dz = finfo[jphi]._z - finfo[iphi]._z;
	if(dz > _minzsep && dz < _maxzsep){
	  double dphi = deltaPhi(finfo[iphi]._phi._val,finfo[jphi]._phi._val);
	  if(fabs(dphi) > _mindphi && fabs(dphi) < _maxdphi){
	    double slope = dphi/dz;
	    if(slope > _smin && slope < _smax){ 
	      double wt = _zweights ? dz/(finfo[iphi]._phi._err+finfo[jphi]._phi._err) : 1.0;
	      accf(slope,weight=wt);
	    }
	  }
	}
      }
    }
    if(count(accf) > _minnhit){
      double  dfdz = extract_result<tag::weighted_median>(accf);
      //      double dfdztest = extract_result<tag::mean>(acctest);
      // if the sign of dfdz disagrees, abort
      if( dfdz * _helicity._value < 0.0)
	return false;
      // if requested, restrict the range
      if(_force)
	dfdz = std::max(std::min(dfdz,_smax),_smin);
      else
	if(dfdz > _smax || dfdz < _smin) return false;
// find phi at z intercept.  Use a histogram technique since phi looping
// hasn't been resolved yet
      TH1F hphi("hphi","phi value",50,-1.1*pi,1.1*pi);
      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
      	double phiex = finfo[iphi]._z*dfdz;
	double dphi = deltaPhi(phiex,finfo[iphi]._phi._val);
	hphi.Fill(dphi);
	hphi.Fill(dphi-twopi);
	hphi.Fill(dphi+twopi);
      }
      double fz0 = hphi.GetBinCenter(hphi.GetMaximumBin());
// refine this using the median.  Not sure this helps, needs to be tested
//      accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > acci;
//      for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
//	double wt = _zweight ? 1.0/finfo[iphi]._phi._err : 1.0;
//	double phiex = fz0 + finfo[iphi]._z*dfdz;
//	double dphi = deltaPhi(phiex,finfo[iphi]._phi._val);
//	acci(dphi, weight=wt);
//      }
//      fz0 += extract_result<tag::weighted_median>(acci);
      if(fz0>pi)fz0 -= twopi;
      if(fz0<-pi)fz0 += twopi;
// reset for full search
      if(_stereoinit){
	finfo.clear();
	for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	  if(xyzp[ixyzp].use()){
	    FZ fz;
	    xyzp[ixyzp].finfo(myhel.center(),fz._phi);
	    fz._z = xyzp[ixyzp]._pos.z();
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
	  int nloop = (int)rint((phiex-finfo[iphi]._phi._val)/twopi);
	  finfo[iphi]._phi._val += nloop*twopi;
	  changed |= nloop != 0;
	}
	// make a long-range estimate of slope
	accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf2;
//	accumulator_set<double, stats<tag::mean > > acctest2;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  for(unsigned jphi=iphi+1; jphi < finfo.size(); ++jphi){
	    double dz = finfo[jphi]._z -finfo[iphi]._z;
	    double dphi = finfo[jphi]._phi._val-finfo[iphi]._phi._val;
	    if(dz > _minzsep && fabs(dphi) > _mindphi ){
	      double dphiex = dz*dfdz;
	      double ferr = finfo[iphi]._phi._err+finfo[jphi]._phi._err;
	      if(!_filterz || fabs(dphi-dphiex) < _nsigma*ferr){
		double slope = dphi/dz;
		if(slope > _smin && slope < _smax){ 
// limit the weight so as not to count more than 1 loop
		  double wt = _zweights ? std::min(dz,_maxzsep)/ferr : 1.0;
		  accf2(slope,weight=wt);
//		  acctest2(slope);
		}
	      }
	    }
	  }
	}
	dfdz = extract_result<tag::weighted_median>(accf2);
//	dfdztest = extract_result<tag::mean>(acctest2);
	// find phi at z intercept
	accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > acci2;
//	accumulator_set<double, stats<tag::mean>> acci2test;
	for(unsigned iphi=0; iphi < finfo.size(); ++iphi){
	  double phiex = fz0+finfo[iphi]._z*dfdz;
	  if(!_filterz || fabs(finfo[iphi]._phi._val-phiex) < _nsigma*finfo[iphi]._phi._err){
	    double wt = _zweights ? 1.0/finfo[iphi]._phi._err : 1.0;
	    acci2(finfo[iphi]._phi._val - finfo[iphi]._z*dfdz, weight=wt);
//	    acci2test(finfo[iphi]._phi._val - finfo[iphi]._z*dfdz);
	  }
	}
//	double fz0test = fmod(extract_result<tag::mean>(acci2test),twopi);
	fz0 = fmod(extract_result<tag::weighted_median>(acci2),twopi);
	if(fz0>pi)fz0 -= twopi;
	if(fz0<-pi)fz0 += twopi;
      }
      // inverse convention: should solve directly for lambda, FIXME!
      myhel.lambda() = 1.0/dfdz;
      myhel.fz0() = fz0;
      // fix the phi for the hit points
      for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
	double phiex = myhel.fz0() + xyzp[ixyzp]._pos.z()/myhel.lambda();
	FZ fz;
	xyzp[ixyzp].finfo(myhel.center(),fz._phi);
	int nloop = (int)rint((phiex - fz._phi._val)/twopi);
	xyzp[ixyzp]._phi = fz._phi._val + nloop*twopi;
	if(_filterz && fabs(xyzp[ixyzp]._phi-phiex)> _nsigma*fz._phi._err) xyzp[ixyzp].setOutlier();
      }
      return true;
    } else
      return false;
  }


  bool
  RobustHelixFit::initCircle(XYZPVector const& xyzp,RobustHelix& myhel) {
    bool retval(false);
    static const double mind2 = _mindist*_mindist;
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx, accy, accr;
    // form all triples, and compute the circle center for unaligned hits.  I can aford to be choosy
    unsigned ntriple(0);
    unsigned nxyzp = xyzp.size();
    std::vector<CLHEP::Hep3Vector> pos;
    pos.reserve(xyzp.size());
    for(size_t ixyzp=0; ixyzp < nxyzp; ++ixyzp){
      if(xyzp[ixyzp].use() && (xyzp[ixyzp].stereo() || (!_stereoinit)) ){
	pos.push_back(xyzp[ixyzp]._pos);
      }
    }
    size_t np = pos.size();
    for(size_t ip=0;ip<np;++ip){
      // pre-compute some values
      double ri2 = pow(pos[ip].x(),2) + pow(pos[ip].y(),2);
      for(size_t jp=ip+1;jp<np; ++jp){
	if(pos[ip].perpPart().diff2(pos[jp].perpPart()) > mind2){
	  double rj2 = pow(pos[jp].x(),2) + pow(pos[jp].y(),2);
	  size_t mink = jp+1;
  // we can force the initialization to use the target position in every triple
	  for(size_t kp=mink;kp<np; ++kp){
	    if(pos[ip].perpPart().diff2(pos[kp].perpPart()) > mind2 &&
		pos[jp].perpPart().diff2(pos[kp].perpPart()) > mind2){
	      // this effectively measures the slope difference
	      double delta = (pos[kp].x() - pos[jp].x())*(pos[jp].y() - pos[ip].y()) - 
		(pos[jp].x() - pos[ip].x())*(pos[kp].y() - pos[jp].y());
	      if(fabs(delta) > _mindelta){
		double rk2 = pow(pos[kp].x(),2) + pow(pos[kp].y(),2);
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
		if(rho > _rmin && rho<_rmax && rmax < _trackerradius
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
      double centx = extract_result<tag::median>(accx);
      double centy = extract_result<tag::median>(accy);
      double rho = extract_result<tag::median>(accr);
      myhel.center() = CLHEP::Hep3Vector(centx,centy,0.0);
      myhel.radius() = rho;
      if(_force){
	myhel.radius() = std::max(std::min(myhel.radius(),_rmax),_rmin);
	retval = true;
      } else {
	retval = myhel.radius() >= _rmin && myhel.radius() <= _rmax;
      }
    }
    return retval;
  }

  void fillXYZP(StrawHitCollection const& shcol,
    StrawHitPositionCollection const& shpcol, std::vector<hitIndex> hits, XYZPVector& xyzp) {
    const Tracker& tracker = getTrackerOrThrow();
    // loop over straw hits, and store their positions
    for(auto istr : hits) { 
      StrawHit const& sh = shcol.at(istr._index);
      Straw const& straw= tracker.getStraw(sh.strawIndex());
      StrawHitPosition const& shp = shpcol.at(istr._index);
      XYZP pos(istr._index,sh,shp,straw);
      xyzp.push_back(pos);
    } 
  }

  void
    RobustHelixFit::findAGE(XYZPVector const& xyzp, Hep3Vector const& center,double& rmed, double& age,bool useweights) {
      using namespace boost::accumulators;
      // protection against empty data
      if(xyzp.size() == 0)return;
      // fill radial information for all points, given this center
      std::vector<VALERR> radii;
      unsigned nxyzp = xyzp.size();
      double wtot(0.0);
      for(unsigned ixyzp=0; ixyzp < nxyzp; ++ixyzp){
	if(xyzp[ixyzp].use()){
	// find radial information for this point
	VALERR rad;
	xyzp[ixyzp].rinfo(center,rad);
	radii.push_back(rad);
	// compute the normalization too
	wtot += useweights ? 1.0/rad._err : 1.0;
      }
    }
    // find the weighted median radius
    accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accr;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._err : 1.0;
      accr(radii[irad]._val, weight = wt); 
    }
    rmed = extract_result<tag::weighted_median>(accr);
    // if requested force radius into range
    if(_force)
      rmed = std::max(std::min(rmed,_rmax),_rmin);
    // now compute the AGE (Absolute Geometric Error)
    age = 0.0;
    for(unsigned irad=0;irad<radii.size();++irad){
      double wt = useweights ? 1.0/radii[irad]._err : 1.0;
      age += wt*fabs(radii[irad]._val-rmed);
    }
    // normalize
    age *= radii.size()/wtot;
  }

  void
  RobustHelixFit::fillSums(XYZPVector const& xyzp, Hep3Vector const& center,double rmed,SUMS& sums,bool useweights) {
    // initialize sums
    sums.clear();
    // compute the transverse sums
    double wtot(0.0);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
	// find radial information for this point
	VALERR rad;
	xyzp[ixyzp].rinfo(center,rad);
	double wt = useweights ? 1.0/rad._err : 1.0;
	wtot += wt;
	// now x,y projections
	double pcos = (xyzp[ixyzp]._pos.x()-center.x())/rad._val;
	double psin = (xyzp[ixyzp]._pos.y()-center.y())/rad._val;
	// 3 conditions: either the radius is inside the median, outside the median, or 'on' the median.  We define 'on'
	// in terms of the error
	if(fabs(rad._val -rmed) < rad._err ){
	  sums._scc += wt*fabs(pcos);
	  sums._ssc += wt*fabs(psin);
	  ++sums._nc;
	} else if (rad._val > rmed) {
	  sums._sco += wt*pcos;
	  sums._sso += wt*psin;
	  ++sums._no;
	} else {
	  sums._sci += wt*pcos;
	  sums._ssi += wt*psin;        
	  ++sums._ni;
	}
      }  
    }
    // normalize to unit weight
    unsigned nused = sums._nc + sums._no + sums._ni;
    sums._scc *= nused/wtot;
    sums._ssc *= nused/wtot;
    sums._sco *= nused/wtot;
    sums._sso *= nused/wtot;
    sums._sci *= nused/wtot;
    sums._ssi *= nused/wtot;
  }

  void
  RobustHelixFit::filterDist(XYZPVector& xyzp) {
    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(M_PI_2);
    // first, resolve phi.   Use the average X and Y to define the initial
    // phi value, to avoid looping issues
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accx;
    accumulator_set<double, stats<tag::median(with_p_square_quantile) > > accy;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
	accx(xyzp[ixyzp]._pos.x());
	accy(xyzp[ixyzp]._pos.y());
      }
    }
    double mx = extract_result<tag::median>(accx);
    double my = extract_result<tag::median>(accy);
    double mphi = atan2(my,mx);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      if(xyzp[ixyzp].use()){
	double dphi = xyzp[ixyzp]._phi - mphi;
	if(fabs(dphi) > pi){
	  if(dphi > 0)
	    xyzp[ixyzp]._phi -= twopi;
	  else
	    xyzp[ixyzp]._phi += twopi;
	}
      }
    }
    // now cut
    CLHEP::Hep3Vector mh(mx,my,0.0);
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      double dist = sqrt(xyzp[ixyzp]._pos.perpPart().diff2(mh));
      if(dist > _maxdist)xyzp[ixyzp].setOutlier();
    }
  }

  void
  RobustHelixFit::filterXY(XYZPVector& xyzp, Hep3Vector const& center,double rmed,bool& changed) {
    static StrawHitFlag other(StrawHitFlag::other);
    changed = false;
    for(unsigned ixyzp=0; ixyzp < xyzp.size(); ++ixyzp){
      VALERR rad;
      xyzp[ixyzp].rinfo(center,rad);
      bool olduse = xyzp[ixyzp].use();
  // update the flag each iteration, until we've converged
      if(fabs(rad._val -rmed) > _nsigma*rad._err)
	xyzp[ixyzp]._flag.merge(other);
      else
	xyzp[ixyzp]._flag.clear(other);
      changed |= olduse != xyzp[ixyzp].use();
    }
  }

  double
  RobustHelixFit::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }
}

