//
// Object to perform helix fit to straw hits
//
// $Id: HelixFit.cc,v 1.12 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $
// $Date: 2014/07/10 14:47:26 $
//
// mu2e
#include "TrkReco/inc/RobustHelixFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include "RecoDataProducts/inc/CaloCluster.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
// root
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"
//c++
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <cmath>

using namespace std;
using namespace boost::accumulators;
using namespace ROOT::Math::VectorUtil;

namespace {
  typedef ROOT::Math::XYVectorF  XYVec;
  // struct for weighted positions
  class XYWVec : public XYVec {
    public :
      XYWVec(XYZVec pos,float weight=1.0) : XYVec(pos.x(),pos.y()), _weight(weight) {}
      float weight() const { return _weight; }

    private :
      float _weight; // weight for this position
  };

}
namespace mu2e
{
 
 
  typedef std::pair<float,float> WVal;

  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _cinit(static_cast<CircleFit>(pset.get<int>("CircleInitType",median))),
    _cfit(static_cast<CircleFit>(pset.get<int>("CircleFitType",median))),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnhit(pset.get<unsigned>("minNHit",5)),
    _lambda0(pset.get<float>("lambda0",0.1)),
    _lstep(pset.get<float>("lstep",0.01)),
    _minlambda(pset.get<float>("minlambda",0.001)),
    _nphibins(pset.get<unsigned>("NPhiHistBins",25)),
    _phifactor(pset.get<float>("PhiHistRangeFactor",1.2)),
    _minnphi(pset.get<unsigned>("MinNPhi",5)),
    _maxniter(pset.get<unsigned>("maxniter",100)),
    _minzsep(pset.get<float>("minzsep",100.0)),
    _maxzsep(pset.get<float>("maxzsep",500.0)),
    _mindphi(pset.get<float>("mindphi",0.5)),
    _maxdphi(pset.get<float>("maxdphi",2.5)),
    _mindist(pset.get<float>("mindist",100.0)), // mm
    _maxdist(pset.get<float>("maxdist",500.0)), // mm
    _rmin(pset.get<float>("minR",150.0)), // mm
    _rmax(pset.get<float>("maxR",400.0)), // mm
    _rcmin(pset.get<float>("minCenterR",200.0)), // mm
    _rcmax(pset.get<float>("maxCenterR",400.0)), // mm
//    _mindelta(pset.get<float>("minDelta",500.0)),
    _lmin(pset.get<float>("minAbsLambda",100.0)),
    _lmax(pset.get<float>("maxAbsLambda",400.0)),
    _targetcon(pset.get<bool>("targeconsistent",true)),
    _targetinter(pset.get<bool>("targetintersect",false)),
    _tripler(pset.get<bool>("TripleRadius",false)),
    _errrwt(pset.get<bool>("HitErrorWeight",false)),
    _usecc(pset.get<bool>("UseCaloCluster",false)),
    _ccwt(pset.get<float>("CaloClusterWeight",10.0)), // Cluster weight in units of non-stereo hits
    _targetradius(pset.get<float>("targetradius",150.0)), // effective target radius (mm)
    _trackerradius(pset.get<float>("trackerradius",750.0)), // tracker out radius; include some buffer (mm)
    _rwind(pset.get<float>("RadiusWindow",10.0)), // window for calling a point to be 'on' the helix in the AGG fit (mm)
    _hphi("hphi","phi value",_nphibins,-_phifactor*CLHEP::pi,_phifactor*CLHEP::pi),
    _ntripleMin(pset.get<unsigned>("ntripleMin",10)),
    _ntripleMax(pset.get<unsigned>("ntripleMax",100000))
  {
    float minarea(pset.get<float>("minArea",2000.0));
    _minarea2 = minarea*minarea;
  }

  RobustHelixFit::~RobustHelixFit()
  {}


  void RobustHelixFit::fitHelix(HelixSeed& hseed) {
    hseed._status.clear(TrkFitFlag::helixOK);

    fitCircle(hseed);
    if (hseed._status.hasAnyProperty(TrkFitFlag::circleOK))
    {
      fitFZ(hseed);
      if (goodHelix(hseed._helix)) hseed._status.merge(TrkFitFlag::helixOK);
    }
  }

  bool RobustHelixFit::initCircle(HelixSeed& hseed) {
    bool retval(false);

    switch ( _cinit ) {
      case median : default :
	fitCircleMedian(hseed);
	retval = hseed._status.hasAllProperties(TrkFitFlag::circleOK);
	break;
    }
    return retval;
  }

  void RobustHelixFit::fitCircle(HelixSeed& hseed) {
    hseed._status.clear(TrkFitFlag::circleOK);

    // if required, initialize
    bool init(false);
    if (!hseed._status.hasAllProperties(TrkFitFlag::circleInit)) {
      init = true;
      if (initCircle(hseed))
	hseed._status.merge(TrkFitFlag::circleInit);
      else
	return;
    }
    // if we initialized and the initialization is the same as the fit type, we're done
    // make sure to refit (iterate) otherwise
    if(!init || _cfit != _cinit) {
      switch ( _cfit ) {
	case median : default :
	  fitCircleMedian(hseed);
	  break;
	case mean :
	  fitCircleMean(hseed);
	  break;
	case AGE :
	  fitCircleAGE(hseed);
	  break;
      }
    }
  }

  void RobustHelixFit::fitCircleMean(HelixSeed& hseed) {

  }

  void RobustHelixFit::fitCircleAGE(HelixSeed& hseed) {
    // this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
    // this algorithm needs extension to use the calorimeter cluster position FIXME!

    RobustHelix& rhel = hseed._helix;

    unsigned niter(0);
    float age;
    XYZVec center = rhel.center();
    float rmed = rhel.radius();
    // initialize step
    float lambda = _lambda0;
    // find median and AGE for the initial center
    findAGE(hseed,center,rmed,age);
    // loop while step is large
    XYZVec descent(1.0,0.0,0.0);
    while(lambda*sqrtf(descent.mag2()) > _minlambda && niter < _maxniter)
    {
      // fill the sums for computing the descent vector
      AGESums sums;
      fillSums(hseed,center,rmed,sums);
      // descent vector cases: if the inner vs outer difference is significant (compared to the median), damp using the median sums,
      // otherwise not.  These expressions take care of the undiferentiable condition on the boundary.
      float dx(sums._sco-sums._sci);
      float dy(sums._sso-sums._ssi);
      if(fabs(dx) < sums._scc)
	dx += (sums._sco < sums._sci) ? -sums._scc : sums._scc;
      if(fabs(dy) < sums._ssc)
	dy += (sums._sso < sums._ssi) ? -sums._ssc : sums._ssc;
      descent = XYZVec(dx,dy,0.0);
      // compute error function, decreasing lambda until this is better than the previous
      float agenew;
      XYZVec cnew = center + lambda*descent;
      findAGE(hseed,cnew,rmed,agenew);
      // if we've improved, increase the step size and iterate
      if(agenew < age){
	lambda *= (1.0+_lstep);
      } else {
	// if we haven't improved, keep reducing the step till we do
	unsigned miter(0);
	while(agenew > age && miter < _maxniter && lambda*sqrtf(descent.mag2()) > _minlambda){
	  lambda *= (1.0-_lstep);
	  cnew = center + lambda*descent;
	  findAGE(hseed,cnew,rmed,agenew);
	  ++miter;
	}
	// if this fails, reverse the descent drection and try again
	if(agenew > age){
	  descent *= -1.0;
	  lambda *= (1.0 +_lstep);
	  cnew = center + lambda*descent;
	  findAGE(hseed,cnew,rmed,agenew);
	}
      }
      // prepare for next iteration
      if(agenew < age){
	center = cnew;
	age = agenew;
      } else {
	static const float minage(0.1);
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
    rhel._rcent = sqrtf(center.perp2());
    rhel._fcent = center.phi();
    rhel._radius = rmed;
    // update flag
    if(goodCircle(rhel)) hseed._status.merge(TrkFitFlag::circleOK);
  }


  void RobustHelixFit::forceTargetInter(XYZVec& center, float& radius) {    
    float rperigee = sqrtf(center.perp2())-radius;    
    if (fabs(rperigee) > _targetradius)
    {
      // adjust both center position and radius till they touch the target, holding phi constant.  
      // This keeps the circle near the hits.  Sign matters!
      float dr;
      if(rperigee > 0)
	dr = 0.5*(rperigee - _targetradius);  // change is 1/2 the difference
      else
	dr = 0.5*(rperigee + _targetradius);
      radius += dr;

      // direction radially outwards from origin to the center, the center moves opposite the radius
      XYZVec pdir = XYZVec(center.x(),center.y(),0.0).unit();
      center -= dr*pdir;
    }
  }



  bool RobustHelixFit::initFZ(HelixSeed& hseed){
    bool retval(false);
    ComboHitCollection& hhits = hseed._hhits;
    RobustHelix& rhel         = hseed._helix;

    rhel._lambda = 1.0e12; //infinite slope for now
    rhel._fz0 = 0.0;
    static TrkFitFlag circleOK(TrkFitFlag::circleOK);
    static TrkFitFlag helixOK(TrkFitFlag::helixOK);

    std::sort(hhits.begin(),hhits.end(),[](const ComboHit& p1, const ComboHit& p2){return p1._pos.z() < p2._pos.z();});    
    for(auto& hhit : hhits)initPhi(hhit,rhel); 

    std::vector<const ComboHit*> validHhits;
    validHhits.reserve(hhits.size());
    for (const auto& hhit: hhits) 
      if (use(hhit) ) validHhits.push_back(&hhit);
    if (validHhits.empty()) return retval;

    // make initial estimate of dfdz using 'nearby' pairs.  This insures they are on the same loop
    accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accf;
    for (auto ihit=validHhits.begin(); ihit != std::prev(validHhits.end()); ++ihit) 
    {
      for(auto jhit = std::next(ihit); jhit != validHhits.end(); ++jhit)
      {
	float dz = (*jhit)->_pos.z() - (*ihit)->_pos.z();
	if (dz < _minzsep || dz > _maxzsep) continue;
	float dphi = (*jhit)->_hphi-(*ihit)->_hphi;
	if (dphi*int(rhel.helicity()._value) < 0 || fabs(dphi) < _mindphi || fabs(dphi) > _maxdphi) continue;

	float lambda = dz/dphi;
	if(goodLambda(rhel.helicity(),lambda)){
	  float wt = sqrtf((*ihit)->nStrawHits()*(*jhit)->nStrawHits());
	  //		float wt = (*ihit)->nStrawHits() + (*jhit)->nStrawHits();
	  accf(lambda, weight=wt);
	}
      }
    }

    if(boost::accumulators::extract::count(accf) < _minnhit) return retval;

    float lambda = extract_result<tag::weighted_median>(accf);

    if(!goodLambda( rhel.helicity(),lambda) ) return retval;
    rhel._lambda = lambda;

    // find phi at z intercept.  Use a histogram technique since phi looping
    // hasn't been resolved yet, and to avoid inefficiency at the phi wrapping edge
    _hphi.Reset();
    for(const auto& hhit : validHhits) 
    {
      float phiex = rhel.circleAzimuth(hhit->_pos.z());
      float dphi = deltaPhi(phiex,hhit->helixPhi());
      _hphi.Fill(dphi);
      _hphi.Fill(dphi-CLHEP::twopi);
      _hphi.Fill(dphi+CLHEP::twopi);
    }

    // take the average of the maximum bin +- 1
    int imax = _hphi.GetMaximumBin();
    float count(0.0);
    float fz0(0.0);
    for (int ibin=std::max((int)0,imax-1); ibin <= std::min((int)imax+1,(int)_nphibins); ++ibin)
    {
      count += _hphi.GetBinContent(ibin);
      fz0 += _hphi.GetBinContent(ibin)*_hphi.GetBinCenter(ibin);
    }

    if(count > _minnphi)
    {
      fz0 /= count;
      rhel._fz0 = deltaPhi(0.0,fz0);
      for (auto& hhit : hhits) resolvePhi(hhit,rhel);
      retval = true;
    }

    return retval;
  }



  void RobustHelixFit::fitFZ(HelixSeed& hseed) {
    // if required, initialize
    hseed._status.clear(TrkFitFlag::phizOK);
    if (!hseed._status.hasAllProperties(TrkFitFlag::phizInit))
    {
      if (initFZ(hseed))
	hseed._status.merge(TrkFitFlag::phizInit);
      else
	return;
    }

    ComboHitCollection& hhits = hseed._hhits;
    RobustHelix& rhel         = hseed._helix;

      std::vector<const ComboHit*> validHhits;
      validHhits.reserve(hhits.size());
      for (const auto& hhit: hhits) if (use(hhit)) validHhits.push_back(&hhit);
      if (validHhits.empty()) return;


      //iterate over lambda and loop resolution
      unsigned niter(0);
      bool changed(true);
      while(changed && niter < _maxniter)
      {
          changed = false;
          accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accf;

          for (auto ihit=validHhits.begin(); ihit != std::prev(validHhits.end()); ++ihit) 
          {
              for(auto jhit = std::next(ihit); jhit != validHhits.end(); ++jhit)
              {
	          float dz = (*jhit)->_pos.z() - (*ihit)->_pos.z();
	          float dphi = (*jhit)->helixPhi()-(*ihit)->helixPhi(); 
	          if (dz < _minzsep || fabs(dphi) < _mindphi) continue;

                  float lambda = dz/dphi;
	          if (goodLambda(rhel.helicity(),lambda)){
		    float wt = sqrtf((*ihit)->nStrawHits()*(*jhit)->nStrawHits());
//		    float wt = (*ihit)->nStrawHits() + (*jhit)->nStrawHits();
		    accf(lambda, weight=wt);
		  }
	      }
          } 
          rhel._lambda = extract_result<tag::weighted_median>(accf);

          // now extract intercept.  Here we solve for the difference WRT the previous value
          accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > acci;
          for (const auto & hhit : validHhits)
          {
	      float phiex = rhel.circleAzimuth(hhit->_pos.z());
	      float dphi = deltaPhi(phiex,hhit->helixPhi());
	      float wt = hhit->nStrawHits();
	      acci(dphi,weight = wt);// accumulate the difference WRT the current intercept
          }
          // enforce convention on azimuth phase
          float dphi = extract_result<tag::weighted_median>(acci);
          rhel._fz0 = deltaPhi(0.0,rhel.fz0()+ dphi);

          // resolve the hit loops again
          for (auto& hhit : hhits) changed |= resolvePhi(hhit,rhel);
          ++niter;
      }

      if (goodFZ(rhel)) hseed._status.merge(TrkFitFlag::phizOK);
  }



  // simple median fit.  No initialization required
  void RobustHelixFit::fitCircleMedian(HelixSeed& hseed) 
  {
     const float mind2 = _mindist*_mindist;
     const float maxd2 = _maxdist*_maxdist;
      
     ComboHitCollection& hhits = hseed._hhits;
     RobustHelix& rhel         = hseed._helix;
     accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accx, accy, accr;
     //float xcmean(0),ycmean(0),rcmean(0),sumWeights(0);


     // pick out a subset of hits. I can aford to be choosy
     std::vector<XYWVec> wpos;
     wpos.reserve(hhits.size()+1);

     for (const auto& hhit : hhits) {
       if (use(hhit) ){
	 wpos.push_back(XYWVec(hhit.pos(),hhit.nStrawHits()));
       }
     }
     XYZVec cog;
     if ( _usecc && hseed.caloCluster().isNonnull())
     {
       mu2e::GeomHandle<mu2e::Calorimeter> ch;
       const Calorimeter* calo = ch.get();
       cog = Geom::toXYZVec(calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e(hseed.caloCluster()->diskId(),hseed.caloCluster()->cog3Vector())));
       wpos.push_back(XYWVec(cog,_ccwt));
     }

     if (wpos.size()<3) return;

     // loop over all triples
     unsigned ntriple(0);
     size_t np = wpos.size();    
     for(size_t ip=0; ip<np-2; ++ip)
     { 
         float ri2 = wpos[ip].Mag2();
         for(size_t jp=ip+1; jp<np-1; ++jp)
         {                
	     float dist2ij = (wpos[ip]-wpos[jp]).Mag2();
	     if (dist2ij < mind2 || dist2ij > maxd2) continue;	  

             float rj2 = wpos[jp].Mag2();
             for(size_t kp=jp+1; kp<np; ++kp)
             {
	         float dist2ik = (wpos[ip]-wpos[kp]).Mag2();
	         float dist2jk = (wpos[jp]-wpos[kp]).Mag2();
                 if (dist2ik < mind2 ||  dist2jk < mind2 || dist2ik > maxd2 || dist2jk > maxd2) continue;
		 // Heron's formula
		 float area2 = (dist2ij*dist2jk + dist2ik*dist2jk + dist2ij*dist2ik) - 0.5*(dist2ij*dist2ij + dist2jk*dist2jk + dist2ik*dist2ik);

                 // this effectively measures the slope difference
	         float delta = (wpos[kp].x() - wpos[jp].x())*(wpos[jp].y() - wpos[ip].y()) -
		                (wpos[jp].x() - wpos[ip].x())*(wpos[kp].y() - wpos[jp].y());
//	         if (fabs(delta) < _mindelta) continue;
	         if(area2 < _minarea2)continue;

                 float rk2 = wpos[kp].Mag2();

                 // find circle center for this triple
	         float cx = 0.5* (
		     (wpos[kp].y() - wpos[jp].y())*ri2 +
		     (wpos[ip].y() - wpos[kp].y())*rj2 +
		     (wpos[jp].y() - wpos[ip].y())*rk2 ) / delta;
	         float cy = -0.5* (
		     (wpos[kp].x() - wpos[jp].x())*ri2 +
		     (wpos[ip].x() - wpos[kp].x())*rj2 +
		     (wpos[jp].x() - wpos[ip].x())*rk2 ) / delta;
		 XYVec cent(cx,cy);
	         float rho = sqrtf((wpos[ip]-cent).Mag2());
	         float rc = sqrtf(cent.Mag2());
	         float rmin = fabs(rc-rho);
	         float rmax = rc+rho;

        	 // test circle parameters for this triple: should be inside the tracker,
		 // optionally consistent with the target
		 if (rc > _rcmin && rc < _rcmax &&
		    rho > _rmin && rho < _rmax && rmax < _trackerradius && 
                     ( !_targetcon || rmin < _targetradius) )
                 {
		    ++ntriple;

		    float wt = cbrt(wpos[ip].weight()*wpos[jp].weight()*wpos[kp].weight());
//		    float wt = wpos[ip].weight() + wpos[jp].weight() + wpos[kp].weight();
//		    wt *= sqrtf(area2);
		    accx(cx,weight = wt);
		    accy(cy,weight = wt);
		    if(_tripler) accr(rho,weight = wt);
		    if (ntriple>_ntripleMax) {ip=np;jp=np;kp=np;}
		 } 
	     } 
         } 
     } 


     // median calculation needs a reasonable number of points to function
     if (ntriple > _ntripleMin)
     {
        float centx = extract_result<tag::weighted_median>(accx);
        float centy = extract_result<tag::weighted_median>(accy);
	XYVec center(centx,centy);
	if(!_tripler) {
	  if(!_errrwt) {
	    for(auto const& wp : wpos ) {
	      float rho = sqrtf((wp - center).Mag2());
	      accr(rho,weight = wp.weight()); 
	    }
	  } else {
// set weight according to the errors
	    for (const auto& hhit : hhits) {
	      if (use(hhit) ){
	      // compute the projection to the radius
		XYVec rvec = (XYVec(hhit.pos().x(),hhit.pos().y())-center);
		XYVec rdir = rvec.unit();
		float wdot = rdir.Dot(XYVec(hhit.wdir().x(),hhit.wdir().y()));
		float wdot2 = wdot*wdot;
		float tdot2 = 1.0 - wdot2;
		float err2 = wdot2*hhit.wireErr2() + tdot2*hhit.transErr2();
		float wt = 1.0/sqrtf(err2); // or 1/err2?
		float rho = sqrtf(rvec.Mag2());
		accr(rho,weight=wt);
	      }
	    }
	    if ( _usecc && hseed.caloCluster().isNonnull()){
	    }
	  }
	}
	float rho = extract_result<tag::weighted_median>(accr);
        rhel._rcent = sqrtf(center.Mag2());
        rhel._fcent = center.Phi();
        rhel._radius = rho;

        if (goodCircle(rhel)) hseed._status.merge(TrkFitFlag::circleOK);
     }
  
  }



  void RobustHelixFit::findAGE(HelixSeed const& hseed, XYZVec const& center,float& rmed, float& age)
  {     
      const ComboHitCollection& hhits = hseed._hhits;

      // fill radial information for all points, given this center
      std::vector<WVal> radii;
      float wtot(0.0);
      for(auto const& hhit : hhits)
      {
         if(use(hhit))
         {
	    // find radial information for this point
	    float rad = sqrtf(XYZVec(hhit._pos - center).perp2());
	    float wt = hitWeight(hhit);
	    radii.push_back(make_pair(rad,wt));
	    wtot += wt;
         }
      }

      // optionally add calo cluster
      if(_usecc && hseed.caloCluster().isNonnull())
      {
          mu2e::GeomHandle<mu2e::Calorimeter> ch;
          const Calorimeter* calo = ch.get();
          XYZVec cog = Geom::toXYZVec(calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e(hseed.caloCluster()->diskId(),hseed.caloCluster()->cog3Vector())));
          float rad = sqrtf(XYZVec(cog - center).perp2());
          radii.push_back(make_pair(rad,_ccwt));
      }

      // compute AGE
      if (radii.size() > _minnhit)
      {
        // find the median radius
        accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accr;
        for(unsigned irad=0;irad<radii.size();++irad)
	  accr(radii[irad].first, weight = radii[irad].second);

        rmed = extract_result<tag::weighted_median>(accr);
        // now compute the AGE (Absolute Geometric Error)
        age = 0.0;
        for(unsigned irad=0;irad<radii.size();++irad)
	  age += radii[irad].second*abs(radii[irad].first-rmed);

        // normalize
        age *= radii.size()/wtot;
      }
  }


  void RobustHelixFit::fillSums(HelixSeed const& hseed, XYZVec const& center,float rmed, AGESums& sums)
  {    
     ComboHitCollection const& hhits = hseed._hhits;
     sums.clear();

     float wtot(0.0);
     for (const auto& hhit : hhits)
     {
         if (!use(hhit)) continue;
	 // find radial information for this point
	 float rad = sqrtf(XYZVec(hhit._pos - center).perp2());
	 float wt = hitWeight(hhit);
	 // now x,y projections
	 float pcos = (hhit._pos.x()-center.x())/rad;
	 float psin = (hhit._pos.y()-center.y())/rad;
	 // 3 conditions: either the radius is inside the median, outside the median, or 'on' the median.  We define 'on'
	 // in terms of a window
	 if(fabs(rmed - rad) < _rwind  ){
	   sums._scc += wt*fabs(pcos);
	   sums._ssc += wt*fabs(psin);
	   ++sums._nc;
	 } else if (rad > rmed  ) {
	   sums._sco += wt*pcos;
	   sums._sso += wt*psin;
	   ++sums._no;
	 } else {
	   sums._sci += wt*pcos;
	   sums._ssi += wt*psin;
	   ++sums._ni;
	 }
	 wtot += wt;
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

  float RobustHelixFit::deltaPhi(float phi1, float phi2)
  {
     float dphi = fmod(phi2-phi1,CLHEP::twopi);
     if (dphi>CLHEP::pi) dphi -= CLHEP::twopi;
     if (dphi<-CLHEP::pi)dphi += CLHEP::twopi;
     return dphi;
  }

  bool RobustHelixFit::use(const ComboHit& hhit) const 
  {
     return (!hhit._flag.hasAnyProperty(_dontuseflag));
  }

  void RobustHelixFit::initPhi(ComboHit& hhit, const RobustHelix& rhel) const 
  {
     // ray from the circle center to the point
     hhit._hphi = XYZVec(hhit._pos - rhel.center()).phi();
  }

  bool RobustHelixFit::resolvePhi(ComboHit& hhit, const RobustHelix& rhel) const
  {
     // find phi expected
     float phiex = rhel.circleAzimuth(hhit._pos.z());
     int nloop = (int)rint((phiex-hhit.helixPhi())/CLHEP::twopi);
     hhit._hphi += nloop*CLHEP::twopi;

     static StrawHitFlag resphi(StrawHitFlag::resolvedphi);
     hhit._flag.merge(resphi);
     return nloop != 0;
  }

  bool RobustHelixFit::goodFZ(const RobustHelix& rhel)
  {
     return rhel.validHelicity() && fabs(rhel.lambda()) > _lmin && fabs(rhel.lambda()) < _lmax;
  }

  bool RobustHelixFit::goodCircle(const RobustHelix& rhel)
  {
     return rhel.radius() > _rmin && rhel.radius() < _rmax;
  }

  bool RobustHelixFit::goodHelix(const RobustHelix& rhel)
  {
     return goodCircle(rhel) && goodFZ(rhel);
  }

  float RobustHelixFit::hitWeight(const ComboHit& hhit) const 
  {
     float retval(hhit.nStrawHits());
     // add an option to evaluate the error relative to the current center FIXME
     return retval;
  }

  bool RobustHelixFit::goodLambda(Helicity const& h, float lambda) const {
    bool retval(false);
    switch (h._value) {
      case Helicity::neghel :
	retval = (lambda > -_lmax) && (lambda < -_lmin);
	break;
      case Helicity::poshel :
	retval = (lambda > _lmin) && (lambda < _lmax);
	break;
      default:
	break;
    }
    return retval;
  }
}
