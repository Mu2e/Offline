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
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>

#include "TH1F.h"

#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <cmath>

using namespace std;
using namespace boost::accumulators;


namespace {
  
  // struct for weighted positions
  class WPos : public CLHEP::Hep3Vector {
    public :
      WPos(CLHEP::Hep3Vector pos,double weight=1.0) : CLHEP::Hep3Vector(pos), _weight(weight) {}
      double weight() const { return _weight; }
    private :
      CLHEP::Hep3Vector _pos; //position
      double _weight; // weight for this position
  };

}
namespace mu2e
{
 
  typedef std::pair<double,double> WVal;

  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _cfit(static_cast<CircleFit>(pset.get<int>("CircleFitType",median))),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnhit(pset.get<unsigned>("minNHit",5)),
    _lambda0(pset.get<double>("lambda0",0.1)),
    _lstep(pset.get<double>("lstep",0.01)),
    _minlambda(pset.get<double>("minlambda",0.001)),
    _nphibins(pset.get<unsigned>("NPhiHistBins",25)),
    _phifactor(pset.get<double>("PhiHistRangeFactor",1.2)),
    _minnphi(pset.get<unsigned>("MinNPhi",5)),
    _maxniter(pset.get<unsigned>("maxniter",100)),
    _minzsep(pset.get<double>("minzsep",100.0)),
    _maxzsep(pset.get<double>("maxzsep",500.0)),
    _mindphi(pset.get<double>("mindphi",0.5)),
    _maxdphi(pset.get<double>("maxdphi",2.5)),
    _mindist(pset.get<double>("mindist",50.0)), // mm
    _maxdist(pset.get<double>("maxdist",500.0)), // mm
    _rmin(pset.get<double>("minR",150.0)), // mm
    _rmax(pset.get<double>("maxR",400.0)), // mm
    _mindelta(pset.get<double>("minDelta",500.0)),
    _lmin(pset.get<double>("minAbsLambda",100.0)),
    _lmax(pset.get<double>("maxAbsLambda",400.0)),
    _stereoinit(pset.get<bool>("stereoinit",false)),
    _stereofit(pset.get<bool>("stereofit",false)),
    _targetinit(pset.get<bool>("targetinit",true)),
    _targetinter(pset.get<bool>("targetintersect",false)),
    _usecc(pset.get<bool>("UseCaloCluster",false)),
    _ccwt(pset.get<double>("CaloClusterWeight",10.0)), // Cluster weight in units of non-stereo hits
    _stwt(pset.get<double>("StereoHitWeight",1.0)), // Stereo hit weight in units of non-stereo hits
    _hqwt(pset.get<bool>("HitQualityWeight",false)), // weight hits by 'quality' = MVA value
    _targetradius(pset.get<double>("targetradius",150.0)), // effective target radius (mm)
    _trackerradius(pset.get<double>("trackerradius",750.0)), // tracker out radius; include some buffer (mm)
    _rwind(pset.get<double>("RadiusWindow",10.0)), // window for calling a point to be 'on' the helix in the AGG fit (mm)
    _helicity(pset.get<int>("Helicity",Helicity::unknown)),
    _hphi("hphi","phi value",_nphibins,-_phifactor*CLHEP::pi,_phifactor*CLHEP::pi),
    _ntripleMax(pset.get<unsigned>("ntripleMax",100000))
  {
    if (_helicity._value == Helicity::unknown)
      throw cet::exception("RECO")<<"mu2e::RobustHelixFit: Invalid Helicity specified"<< std::endl;
    
    if (_stereofit)_useflag = StrawHitFlag(StrawHitFlag::stereo);
    if (_helicity._value < 0)
    {
       // swap order and sign!
       double lmin = -1.0*_lmax;
       _lmax = -1.0*_lmin;
       _lmin = lmin;
    }
     
  }

  RobustHelixFit::~RobustHelixFit()
  {}


  void RobustHelixFit::fitHelix(HelixSeed& hseed, const CaloClusterCollection* caloClusters)
  {
    hseed._status.clear(TrkFitFlag::helixOK);

    fitCircle(hseed, caloClusters);
    if (hseed._status.hasAnyProperty(TrkFitFlag::circleOK))
    {
      fitFZ(hseed);
      if (goodHelix(hseed._helix)) hseed._status.merge(TrkFitFlag::helixOK);
    }
  }
  
  void RobustHelixFit::fitCircle(HelixSeed& hseed, const CaloClusterCollection* caloClusters)
  {
     hseed._status.clear(TrkFitFlag::circleOK);
     switch ( _cfit )
     {
        case median : default :
	  fitCircleMedian(hseed, caloClusters);
	  break;
        case AGE :
	  fitCircleAGE(hseed, caloClusters);
	  break;
        case chisq :
	  throw cet::exception("Reco") << "mu2e::RobustHelixFit: Chisq fit not yet implemented " << std::endl;
	  break;
     }
  }

  void RobustHelixFit::fitCircleAGE(HelixSeed& hseed, const CaloClusterCollection* caloClusters) 
  {
    // this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
    // This algorithm requires a reasonable initial estimate: use the median fit
    // this algorithm needs extension to use the calorimeter cluster position FIXME!
    if (!hseed._status.hasAllProperties(TrkFitFlag::circleInit)) fitCircleMedian(hseed, caloClusters);
    
    if (hseed._status.hasAllProperties(TrkFitFlag::circleInit))
    {
      RobustHelix& rhel = hseed._helix;

      unsigned niter(0);
      double age;
      CLHEP::Hep3Vector center = rhel.center();
      double rmed = rhel.radius();
      // initialize step
      double lambda = _lambda0;
      // find median and AGE for the initial center
      findAGE(hseed,center,rmed,age);
      // loop while step is large
      CLHEP::Hep3Vector descent(1.0,0.0,0.0);
      while(lambda*descent.mag() > _minlambda && niter < _maxniter)
      {
	// fill the sums for computing the descent vector
	AGESums sums;
	fillSums(hseed,center,rmed,sums);
	// descent vector cases: if the inner vs outer difference is significant (compared to the median), damp using the median sums,
	// otherwise not.  These expressions take care of the undiferentiable condition on the boundary.
	double dx(sums._sco-sums._sci);
	double dy(sums._sso-sums._ssi);
	if(fabs(dx) < sums._scc)
	  dx += (sums._sco < sums._sci) ? -sums._scc : sums._scc;
	if(fabs(dy) < sums._ssc)
	  dy += (sums._sso < sums._ssi) ? -sums._ssc : sums._ssc;
	descent = CLHEP::Hep3Vector(dx,dy,0.0);
	// compute error function, decreasing lambda until this is better than the previous
	double agenew;
	CLHEP::Hep3Vector cnew = center + lambda*descent;
	findAGE(hseed,cnew,rmed,agenew);
	// if we've improved, increase the step size and iterate
	if(agenew < age){
	  lambda *= (1.0+_lstep);
	} else {
	  // if we haven't improved, keep reducing the step till we do
	  unsigned miter(0);
	  while(agenew > age && miter < _maxniter && lambda*descent.mag() > _minlambda){
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
      rhel._rcent = center.perp();
      rhel._fcent = center.phi();
      rhel._radius = rmed;
      // update flag
      if(goodCircle(rhel)) hseed._status.merge(TrkFitFlag::circleOK);
    }
  }


  void RobustHelixFit::forceTargetInter(CLHEP::Hep3Vector& center, double& radius)
  {    
     double rperigee = center.perp()-radius;    
     if (fabs(rperigee) > _targetradius)
     {
        // adjust both center position and radius till they touch the target, holding phi constant.  
        // This keeps the circle near the hits.  Sign matters!
        double dr;
        if(rperigee > 0)
	  dr = 0.5*(rperigee - _targetradius);  // change is 1/2 the difference
        else
	  dr = 0.5*(rperigee + _targetradius);
        radius += dr;

        // direction radially outwards from origin to the center, the center moves opposite the radius
        CLHEP::Hep3Vector pdir = CLHEP::Hep3Vector(center.x(),center.y(),0.0).unit();
        center -= dr*pdir;
     }
  }



  bool RobustHelixFit::initFZ(HelixHitCollection& hhits,RobustHelix& rhel)
  {      
      bool retval(false);

      rhel._lambda = 1.0e12; //infinite slope for now
      rhel._fz0 = 0.0;
      static TrkFitFlag circleOK(TrkFitFlag::circleOK);
      static TrkFitFlag helixOK(TrkFitFlag::helixOK);

      std::sort(hhits.begin(),hhits.end(),[](const HelixHit& p1, const HelixHit& p2){return p1._pos.z() < p2._pos.z();});    
      for(auto& hhit : hhits)initPhi(hhit,rhel); 


      std::vector<const HelixHit*> validHhits;
      validHhits.reserve(hhits.size());
      for (const auto& hhit: hhits) 
         if (use(hhit) && ( (!_stereoinit) || stereo(hhit))) validHhits.push_back(&hhit);
      if (validHhits.empty()) return retval;


      // make initial estimate of dfdz using 'nearby' pairs.  This insures they are on the same loop
      accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf;
      for (auto ihit=validHhits.begin(); ihit != std::prev(validHhits.end()); ++ihit) 
      {
          for(auto jhit = std::next(ihit); jhit != validHhits.end(); ++jhit)
          {
	      double dz = (*jhit)->_pos.z() - (*ihit)->_pos.z();
	      if (dz < _minzsep || dz > _maxzsep) continue;
              double dphi = (*jhit)->_phi-(*ihit)->_phi; 
	      if (dphi < _mindphi || dphi > _maxdphi) continue;

              double lambda = dz/dphi;
	      if (lambda < _lmin || lambda > _lmax) continue;
     	      double wt = hitWeight(**ihit)*hitWeight(**jhit);
              accf(lambda, weight=wt);
	  }
      } 

      if(boost::accumulators::extract::count(accf) < _minnhit) return retval;

      double lambda = extract_result<tag::weighted_median>(accf);

      if( lambda > _lmax || lambda < _lmin) return retval;
      rhel._lambda = lambda;

      // find phi at z intercept.  Use a histogram technique since phi looping
      // hasn't been resolved yet, and to avoid inefficiency at the phi wrapping edge
      _hphi.Reset();
      for(const auto& hhit : validHhits) 
      {
	  double phiex = rhel.circleAzimuth(hhit->_pos.z());
	  double dphi = deltaPhi(phiex,hhit->_phi);
	  _hphi.Fill(dphi);
	  _hphi.Fill(dphi-CLHEP::twopi);
	  _hphi.Fill(dphi+CLHEP::twopi);
      }

      // take the average of the maximum bin +- 1
      int imax = _hphi.GetMaximumBin();
      double count(0.0);
      double fz0(0.0);
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



  void RobustHelixFit::fitFZ(HelixSeed& hseed) 
  {
      HelixHitCollection& hhits = hseed._hhits;
      RobustHelix& rhel         = hseed._helix;

      // if required, initialize
      hseed._status.clear(TrkFitFlag::phizOK);
      if (!hseed._status.hasAllProperties(TrkFitFlag::phizInit))
      {
        if (initFZ(hhits,rhel))
	  hseed._status.merge(TrkFitFlag::phizInit);
        else
	  return;
      }


      std::vector<const HelixHit*> validHhits;
      validHhits.reserve(hhits.size());
      for (const auto& hhit: hhits) if (use(hhit)) validHhits.push_back(&hhit);
      if (validHhits.empty()) return;


      //iterate over lambda and loop resolution
      unsigned niter(0);
      bool changed(true);
      while(changed && niter < _maxniter)
      {
          changed = false;
          accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accf;

          for (auto ihit=validHhits.begin(); ihit != std::prev(validHhits.end()); ++ihit) 
          {
              for(auto jhit = std::next(ihit); jhit != validHhits.end(); ++jhit)
              {
	          double dz = (*jhit)->_pos.z() - (*ihit)->_pos.z();
	          double dphi = (*jhit)->_phi-(*ihit)->_phi; 
	          if (dz < _minzsep || fabs(dphi) < _mindphi) continue;

                  double lambda = dz/dphi;
	          if (lambda < _lmin || lambda > _lmax) continue;
     	          double wt = hitWeight(**ihit)*hitWeight(**jhit);
                  accf(lambda, weight=wt);
	      }
          } 
          rhel._lambda = extract_result<tag::weighted_median>(accf);

          // now extract intercept.  Here we solve for the difference WRT the previous value
          accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > acci;
          for (const auto & hhit : validHhits)
          {
	      double phiex = rhel.circleAzimuth(hhit->_pos.z());
	      double dphi = deltaPhi(phiex,hhit->_phi);
	      double wt = hitWeight(*hhit);
	      acci(dphi,weight = wt);// accumulate the difference WRT the current intercept
          }
          // enforce convention on azimuth phase
          double dphi = extract_result<tag::weighted_median>(acci);
          rhel._fz0 = deltaPhi(0.0,rhel.fz0()+ dphi);

          // resolve the hit loops again
          for (auto& hhit : hhits) changed |= resolvePhi(hhit,rhel);
          ++niter;
      }

      if (goodFZ(rhel)) hseed._status.merge(TrkFitFlag::phizOK);
  }



  // simple median fit.  No initialization required
  void RobustHelixFit::fitCircleMedian(HelixSeed& hseed, const CaloClusterCollection* caloClusters) 
  {
     const double mind2 = _mindist*_mindist;
     const double maxd2 = _maxdist*_maxdist;
      
     HelixHitCollection& hhits = hseed._hhits;
     RobustHelix& rhel         = hseed._helix;
     accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accx, accy, accr;
     //double xcmean(0),ycmean(0),rcmean(0),sumWeights(0);


     // pick out a subset of hits. I can aford to be choosy
     std::vector<WPos> pos;
     pos.reserve(hhits.size()+1);

     //get the calorimeter hit first
     //int icIdx = hseed.timeCluster()->caloFastIdx();
     //CLHEP::Hep3Vector caloPos(caloClusters->at(icIdx).pos().x(),caloClusters->at(icIdx).pos().y(),0);
     //pos.push_back(WPos(caloPos.perpPart(),_ccwt));



     for (const auto& hhit : hhits) 
         if (use(hhit) && (stereo(hhit) || (!_stereoinit)) )        	                
               pos.push_back(WPos(hhit._pos.perpPart(),hitWeight(hhit)));        

     
     if ( _usecc && hseed.caloCluster().isNonnull())
     {
        mu2e::GeomHandle<mu2e::Calorimeter> ch;
        const Calorimeter* calo = ch.get();
        CLHEP::Hep3Vector cog = calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e(hseed.caloCluster()->diskId(),hseed.caloCluster()->cog3Vector()));
        pos.push_back(WPos(cog.perpPart(),_ccwt));
     }

     if (pos.size()<3) return;

     // loop over all triples
     unsigned ntriple(0);
     size_t np = pos.size();    
     for(size_t ip=0; ip<np-2; ++ip)
     {      
         double ri2 = pos[ip].perp2();
         for(size_t jp=ip+1; jp<np-1; ++jp)
         {                
	     double dist2ij = pos[ip].diff2(pos[jp]);
	     if (dist2ij < mind2 || dist2ij > maxd2) continue;	  

             double rj2 = pos[jp].perp2();
             for(size_t kp=jp+1; kp<np; ++kp)
             {
	         double dist2ik = pos[ip].diff2(pos[kp]);
	         double dist2jk = pos[jp].diff2(pos[kp]);
                 if (dist2ik < mind2 ||  dist2jk < mind2 || dist2ik > maxd2 || dist2jk > maxd2) continue;

                 // this effectively measures the slope difference
	         double delta = (pos[kp].x() - pos[jp].x())*(pos[jp].y() - pos[ip].y()) -
		                (pos[jp].x() - pos[ip].x())*(pos[kp].y() - pos[jp].y());
	         if (fabs(delta) < _mindelta) continue;

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
		 if (rho > _rmin && rho < _rmax && rmax < _trackerradius && 
                     ( !_targetinit || rmin < _targetradius) )
                 {
		    ++ntriple;
		    double wt = pos[ip].weight()*pos[jp].weight()*pos[kp].weight();
		    accx(cx,weight = wt);
		    accy(cy,weight = wt);
		    accr(rho,weight = wt);
                    //xcmean  += cx*wt;
                    //ycmean  += cy*wt;
                    //sumWeights += wt;
                    //rcmean  += rho*wt;

                    //if (ntriple>50) stride=2;		  
		    if (ntriple>_ntripleMax) {ip=np;jp=np;kp=np;}
		 } 
	     } 
         } 
     } 


     // median calculation needs a reasonable number of points to function
     if (ntriple > _minnhit)
     {
        double centx = extract_result<tag::weighted_median>(accx);
        double centy = extract_result<tag::weighted_median>(accy);
        double rho = extract_result<tag::weighted_median>(accr);
        //double centx = xcmean/sumWeights;
        //double centy = ycmean/sumWeights;
        //double rho = rcmean/sumWeights;

        CLHEP::Hep3Vector center(centx,centy,0.0);
        rhel._rcent = center.perp();
        rhel._fcent = center.phi();
        rhel._radius = rho;

        hseed._status.merge(TrkFitFlag::circleInit);
        if (goodCircle(rhel)) hseed._status.merge(TrkFitFlag::circleOK);
     }
  
  }



  void RobustHelixFit::findAGE(HelixSeed const& hseed, CLHEP::Hep3Vector const& center,double& rmed, double& age)
  {     
      const HelixHitCollection& hhits = hseed._hhits;

      // fill radial information for all points, given this center
      std::vector<WVal> radii;
      double wtot(0.0);
      for(auto const& hhit : hhits)
      {
         if(use(hhit))
         {
	    // find radial information for this point
	    double rad = CLHEP::Hep3Vector(hhit._pos - center).perp();
	    double wt = hitWeight(hhit);
	    radii.push_back(make_pair(rad,wt));
	    wtot += wt;
         }
      }

      // optionally add calo cluster
      if(_usecc && hseed.caloCluster().isNonnull())
      {
          mu2e::GeomHandle<mu2e::Calorimeter> ch;
          const Calorimeter* calo = ch.get();
          CLHEP::Hep3Vector cog = calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e(hseed.caloCluster()->diskId(),hseed.caloCluster()->cog3Vector()));
          double rad = CLHEP::Hep3Vector(cog - center).perp();
          radii.push_back(make_pair(rad,_ccwt));
      }

      // compute AGE
      if (radii.size() > _minnhit)
      {
        // find the median radius
        accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > accr;
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


  void RobustHelixFit::fillSums(HelixSeed const& hseed, CLHEP::Hep3Vector const& center,double rmed, AGESums& sums)
  {    
     HelixHitCollection const& hhits = hseed._hhits;
     sums.clear();

     double wtot(0.0);
     for (const auto& hhit : hhits)
     {

         if (!use(hhit)) continue;
	 // find radial information for this point
	 double rad = CLHEP::Hep3Vector(hhit._pos - center).perp();
	 double wt = hitWeight(hhit);
	 // now x,y projections
	 double pcos = (hhit._pos.x()-center.x())/rad;
	 double psin = (hhit._pos.y()-center.y())/rad;
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

  double RobustHelixFit::deltaPhi(double phi1, double phi2)
  {
     double dphi = fmod(phi2-phi1,CLHEP::twopi);
     if (dphi>CLHEP::pi) dphi -= CLHEP::twopi;
     if (dphi<-CLHEP::pi)dphi += CLHEP::twopi;
     return dphi;
  }

  bool RobustHelixFit::use(const HelixHit& hhit) const 
  {
     return (!hhit._flag.hasAnyProperty(_dontuseflag))
      && (hhit._flag.hasAllProperties(_useflag) || _useflag.empty());
  }

  bool RobustHelixFit::stereo(const HelixHit& hhit) const 
  {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return hhit._flag.hasAllProperties(stereo);
  }

  void RobustHelixFit::setOutlier(HelixHit& hhit) const 
  {
     static StrawHitFlag outlier(StrawHitFlag::outlier);
     hhit._flag.merge(outlier);
  }

  void RobustHelixFit::initPhi(HelixHit& hhit, const RobustHelix& rhel) const 
  {
     // ray from the circle center to the point
     hhit._phi = CLHEP::Hep3Vector(hhit._pos - rhel.center()).phi();
  }

  bool RobustHelixFit::resolvePhi(HelixHit& hhit, const RobustHelix& rhel) const
  {
     // find phi expected
     double phiex = rhel.circleAzimuth(hhit._pos.z());
     int nloop = (int)rint((phiex-hhit._phi)/CLHEP::twopi);
     hhit._phi += nloop*CLHEP::twopi;

     static StrawHitFlag resphi(StrawHitFlag::resolvedphi);
     hhit._flag.merge(resphi);
     return nloop != 0;
  }

  bool RobustHelixFit::goodFZ(const RobustHelix& rhel)
  {
     return rhel.lambda() > _lmin && rhel.lambda() < _lmax;
  }

  bool RobustHelixFit::goodCircle(const RobustHelix& rhel)
  {
     return rhel.radius() > _rmin && rhel.radius() < _rmax;
  }

  bool RobustHelixFit::goodHelix(const RobustHelix& rhel)
  {
     return goodCircle(rhel) && goodFZ(rhel) && rhel.helicity() == _helicity;
  }

  double RobustHelixFit::hitWeight(const HelixHit& hhit) const 
  {
     double retval(1.0);
     if (hhit.flag().hasAnyProperty(StrawHitFlag::stereo)) retval = _stwt;
     if (_hqwt) retval*= std::max(0.0,(double)hhit._hqual);
     return retval;
  }


}
