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
#include "TTrackerGeom/inc/TTracker.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"

#include "Mu2eUtilities/inc/polyAtan2.hh"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "boost_fix/accumulators/statistics.hpp"
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
// root
// #include "TH1F.h"
// #include "Math/VectorUtil.h"
// #include "Math/Vector2D.h"
//c++
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <cmath>

using namespace std;
using namespace boost::accumulators;
using namespace ROOT::Math::VectorUtil;

namespace mu2e
{
 
 
  typedef std::pair<float,float> WVal;

  RobustHelixFit::RobustHelixFit(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _cinit(static_cast<CircleFit>(pset.get<int>("CircleInitType",median))),
    _cfit(static_cast<CircleFit>(pset.get<int>("CircleFitType",median))),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnsh(pset.get<unsigned>("minNStrawHits",10)),
    _minnhit(pset.get<unsigned>("minNHit",5)),
    _minxyresid(pset.get<float>("minXYResid",5.)),
    _refineXYFit(pset.get<int>("refineXYFit",0)),
    _refineZPhiFit(pset.get<int>("refineZPhiFit",0)),
    _chi2xymax(pset.get<float>("chi2xymax",5.)),
    _chi2zphimax(pset.get<float>("chi2zphimax",5.)),
    _lambda0(pset.get<float>("lambda0",0.1)),
    _lstep(pset.get<float>("lstep",0.01)),
    _minlambda(pset.get<float>("minlambda",0.001)),
    _mindfdz(pset.get<float>("minDfDz",0.002)),
    _maxdfdz(pset.get<float>("maxDfDz",0.01)),
    _nphibins(pset.get<unsigned>("NPhiHistBins",25)),
    _phifactor(pset.get<float>("PhiHistRangeFactor",1.2)),
    _minnphi(pset.get<unsigned>("MinNPhi",2)),//5)),//FIXME!
    _maxniter(pset.get<unsigned>("maxniter",100)),
    _minzsep(pset.get<float>("minzsep",100.0)),
    _maxzsep(pset.get<float>("maxzsep",500.0)),
    _mindphi(pset.get<float>("mindphi",0.5)),
    _maxdphi(pset.get<float>("maxdphi",2.5)),
    _sigmaPhi(pset.get<float>("sigmaPhi",0.49)),//1636)),// rad
    _mindist(pset.get<float>("mindist",100.0)), // mm
    _maxdist(pset.get<float>("maxdist",500.0)), // mm
    _maxdxy(pset.get<float>("maxdxy",100.0)),
    _maxXDPhi(pset.get<float>("maxXDPhi",5.0)), 
    _rmin(pset.get<float>("minR",160.0)), // mm
    _rmax(pset.get<float>("maxR",320.0)), // mm
    _rcmin(pset.get<float>("minCenterR",140.0)), // mm
    _rcmax(pset.get<float>("maxCenterR",410.0)), // mm
    //    _mindelta(pset.get<float>("minDelta",500.0)),
    _lmin(pset.get<float>("minAbsLambda",130.0)),
    _lmax(pset.get<float>("maxAbsLambda",320.0)),
    _targetcon(pset.get<bool>("targetconsistent",true)),
    _targetinter(pset.get<bool>("targetintersect",false)),
    _tripler(pset.get<bool>("TripleRadius",false)),
    _errrwt(pset.get<bool>("HitErrorWeight",false)),
    _usecc(pset.get<bool>("UseCaloCluster",false)),
    _ccwt(pset.get<float>("CaloClusterWeight",10.0)), // Cluster weight in units of non-stereo hits
    _targetradius(pset.get<float>("targetradius",100.0)), // effective target radius (mm)
    _trackerradius(pset.get<float>("trackerradius",700.0)), // tracker out radius; (mm)
    _rwind(pset.get<float>("RadiusWindow",10.0)), // window for calling a point to be 'on' the helix in the AGG fit (mm)
    _hphi("hphi","phi value",_nphibins,-_phifactor*CLHEP::pi,_phifactor*CLHEP::pi),
    _ntripleMin(pset.get<unsigned>("ntripleMin",5)),
    _ntripleMax(pset.get<unsigned>("ntripleMax",500))
  {
    float minarea(pset.get<float>("minArea",5000.0));
    _minarea2 = minarea*minarea;
  }

  RobustHelixFit::~RobustHelixFit()
  {}


  void RobustHelixFit::fitHelix(RobustHelixFinderData& HelixData) {
    HelixData._hseed._status.clear(TrkFitFlag::helixOK);

    fitCircle(HelixData);
    if (HelixData._hseed._status.hasAnyProperty(TrkFitFlag::circleOK))
      {
	fitFZ(HelixData);

	if (goodHelix(HelixData._hseed._helix)) HelixData._hseed._status.merge(TrkFitFlag::helixOK);
      }
  }

  bool RobustHelixFit::initCircle(RobustHelixFinderData& HelixData) {
    bool retval(false);

    switch ( _cinit ) {
    case median : default :
      fitCircleMedian(HelixData);
      retval = HelixData._hseed._status.hasAllProperties(TrkFitFlag::circleOK);
      break;
    }
    return retval;
  }

  void RobustHelixFit::fitCircle(RobustHelixFinderData& HelixData) {
    HelixData._hseed._status.clear(TrkFitFlag::circleOK);

    // if required, initialize
    bool init(false);
    if (!HelixData._hseed._status.hasAllProperties(TrkFitFlag::circleInit)) {
      init = true;
      if (initCircle(HelixData))
	HelixData._hseed._status.merge(TrkFitFlag::circleInit);
      else
	return;
    }
    // if we initialized and the initialization is the same as the fit type, we're done
    // make sure to refit (iterate) otherwise
    if(!init || _cfit != _cinit) {
      switch ( _cfit ) {
      case median : default :
	fitCircleMedian(HelixData);
	break;
      case mean :
	fitCircleMean(HelixData);
	break;
      case AGE :
	fitCircleAGE(HelixData);
	break;
      }
    }
  }

  void RobustHelixFit::fitCircleMean(RobustHelixFinderData& HelixData) {

  }

  void RobustHelixFit::fitCircleAGE(RobustHelixFinderData& HelixData) {
    // this algorithm follows the method described in J. Math Imagin Vis Dec. 2010 "Robust Fitting of Circle Arcs" (Volume 40, Issue 2, pp. 147-161)
    // this algorithm needs extension to use the calorimeter cluster position FIXME!

    RobustHelix& rhel = HelixData._hseed._helix;

    unsigned niter(0);
    float age;
    XYZVec center = rhel.center();
    float rmed = rhel.radius();
    // initialize step
    float lambda = _lambda0;
    // find median and AGE for the initial center
    findAGE(HelixData,center,rmed,age);
    // loop while step is large
    XYZVec descent(1.0,0.0,0.0);
    while(lambda*sqrtf(descent.mag2()) > _minlambda && niter < _maxniter)
      {
	// fill the sums for computing the descent vector
	AGESums sums;
	fillSums(HelixData,center,rmed,sums);
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
	findAGE(HelixData,cnew,rmed,agenew);
	// if we've improved, increase the step size and iterate
	if(agenew < age){
	  lambda *= (1.0+_lstep);
	} else {
	  // if we haven't improved, keep reducing the step till we do
	  unsigned miter(0);
	  while(agenew > age && miter < _maxniter && lambda*sqrtf(descent.mag2()) > _minlambda){
	    lambda *= (1.0-_lstep);
	    cnew = center + lambda*descent;
	    findAGE(HelixData,cnew,rmed,agenew);
	    ++miter;
	  }
	  // if this fails, reverse the descent drection and try again
	  if(agenew > age){
	    descent *= -1.0;
	    lambda *= (1.0 +_lstep);
	    cnew = center + lambda*descent;
	    findAGE(HelixData,cnew,rmed,agenew);
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
    rhel._fcent = polyAtan2(center.y(), center.x());//center.phi();
    rhel._radius = rmed;
    // update flag
    if(goodCircle(rhel)) HelixData._hseed._status.merge(TrkFitFlag::circleOK);
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



  bool RobustHelixFit::initFZ(RobustHelixFinderData& HelixData, int InitHiPhi){
    bool retval(false);
    // ComboHitCollection& hhits = HelixData._hseed._hhits;
    RobustHelix& rhel         = HelixData._hseed._helix;

    rhel._lambda = 1.0e12; //infinite slope for now
    rhel._fz0 = 0.0;
    static TrkFitFlag circleOK(TrkFitFlag::circleOK);
    static TrkFitFlag helixOK(TrkFitFlag::helixOK);

    ComboHit*    hit(0);

    if (InitHiPhi > 0) {
      for (unsigned i=0; i<HelixData._chHitsToProcess.size(); ++i){ 
	hit = &HelixData._chHitsToProcess[i];
	if (!use(*hit))             continue;
	initPhi(*hit,rhel); 
      }
    }

    // make initial estimate of dfdz using 'nearby' pairs.  This insures they are on the same loop
    accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accf;
    
    ComboHit*      hitP1(0), *hitP2(0);
    uint16_t       facezF1(0), facezF2(0);
    int            nHits(HelixData._chHitsToProcess.size());

    for (int f1=0; f1<nHits-1; ++f1){
      hitP1   = &HelixData._chHitsToProcess[f1];
      if (!use(*hitP1))            continue;
      facezF1 = hitP1->strawId().uniqueFace();
	    
      for (int f2=f1+1; f2<nHits; ++f2){
	hitP2 = &HelixData._chHitsToProcess[f2];
	if (!use(*hitP2))          continue;

	facezF2 = hitP2->strawId().uniqueFace();
	if ( facezF1 == facezF2 )  continue;

	float dz = hitP2->pos().z() - hitP1->pos().z();//facezF2->z - facezF1->z;
	if (dz < _minzsep || dz > _maxzsep)          continue;
	float dphi = deltaPhi(hitP1->_hphi, hitP2->_hphi);
	if (dphi*int(rhel.helicity()._value) < 0 || fabs(dphi) < _mindphi || fabs(dphi) > _maxdphi) continue;

	float lambda = dz/dphi;
	if(goodLambda(rhel.helicity(),lambda)){
	  float wt = sqrtf(hitP1->nStrawHits()*hitP2->nStrawHits());
	  //		float wt = (*ihit)->nStrawHits() + (*jhit)->nStrawHits();
	  accf(lambda, weight=wt);
	}
		
      }//end secondloop over the hits
    }//end first loop over the hits

    if(boost::accumulators::extract::count(accf) < _minnhit) return retval;

    float lambda = extract_result<tag::weighted_median>(accf);

    if(!goodLambda( rhel.helicity(),lambda) ) return retval;
    rhel._lambda = lambda;

    if (_diag){
      HelixData._diag.lambda_0 = rhel._lambda;
    }

    // find phi at z intercept.  Use a histogram technique since phi looping
    // hasn't been resolved yet, and to avoid inefficiency at the phi wrapping edge
    _hphi.Reset();
    for (int f=0; f<nHits; ++f){
      hitP1 = &HelixData._chHitsToProcess[f];
      if (!use(*hitP1) )             continue;   
      
      float phiex = rhel.circleAzimuth(hitP1->pos().z());
      float dphi  = deltaPhi(phiex,hitP1->helixPhi());
      _hphi.Fill(dphi);
      _hphi.Fill(dphi-CLHEP::twopi);
      _hphi.Fill(dphi+CLHEP::twopi);
    }//end loop over the hits


    // take the average of the maximum bin +- 1
    int imax = _hphi.GetMaximumBin();
    // if ( (imax <1) || (imax>_nphibins)) ###put excpetion throw
    float count(0.0);
    float fz0(0.0);
    for (int ibin=std::max((int)0,imax-1); ibin <= std::min((int)imax+1,(int)_nphibins); ++ibin)
      {
	count += _hphi.GetBinContent(ibin);
	fz0   += _hphi.GetBinContent(ibin)*_hphi.GetBinCenter(ibin);
      }
    
    if (_diag){
      HelixData._diag.nfz0counter = count;
    }

    if(count > _minnphi)
      {
	fz0 /= count;
	rhel._fz0 = deltaPhi(0.0,fz0);

	for (int f=0; f<nHits; ++f){
	  hitP1 = &HelixData._chHitsToProcess[f];
	  // if (!use(*hitP1) )                          continue;   
	  resolvePhi(*hitP1,rhel);
	}
	
	retval = true;
      }

    return retval;
  }


//----------------------------------------------------------------------------------------
// 2015-01-13  calculate track DphiDz using histogrammed distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  bool RobustHelixFit::initFZ_2(RobustHelixFinderData& HelixData) {

    float        phi, phi_ref(-1e10), z_ref, dphi, dz, hdfdz(0);

    float        hist[20], minX(0), maxX(0.01), stepX(0.0005), nbinsX(20); // make it 20 bins

    RobustHelix& rhel   = HelixData._hseed._helix;
    XYZVec       center = rhel.center();
    XYZVec       pos_ref;
//-----------------------------------------------------------------------------
// 2017-09-26 gianipez fixed a bug: in case the Helix phi-z fit didn't converge yet, 
// Helix._dfdz is set to -1e6, so we need to make a check here!
// this is a tempOrary fix that doesn't take into account the particle helicity. FIX ME!
//-----------------------------------------------------------------------------
    if (_debug > 5) {
      printf("[RobustHelixFinder::initFZ_2:BEGIN] x0 = %9.3f y0 = %9.3f Helix._radius = %9.3f Helix._nStrawHits = %3i",
	     center.x(), center.y(), HelixData._hseed._helix._radius,HelixData._nStrawHits);
    }

    int       nstations, nhits[30], nstations_with_hits(0);
    float     phiVec[30], zVec[30], weight(0);
    
    nstations = _tracker->nStations();//RobustHelixFinderData::kNStations;

    for (int i=0; i<nstations; i++) {
      phiVec[i] = 0;
      zVec  [i] = 0;
      nhits [i] = 0;
    }

    for (int i=0; i<nbinsX; i++) hist[i] = 0;

    //-----------------------------------------------------------------------------
    // Step 1: for each station with track candidate hits, calculate average phi per station
    //-----------------------------------------------------------------------------
    FaceZ_t*   facez(0);
    PanelZ_t*  panelz(0);
    
    for (int f1=0; f1<StrawId::_ntotalfaces; ++f1){
      facez = &HelixData._oTracker[f1];
      for (int p1=0; p1<FaceZ_t::kNPanels; ++p1){
	panelz  = &facez->panelZs[p1];   
	int  nhitsPanelF1  = panelz->nChHits();
	if (nhitsPanelF1 == 0)                                     continue;

	for (int i=0; i<nhitsPanelF1; ++i){   
	  ComboHit* hit = &HelixData._chHitsToProcess[panelz->idChBegin + i];

	  // if (hit->_flag.hasAnyProperty(_dontuseflag))                  continue;
	  if (!use(*hit))   continue;


	  int ist = hit->strawId().station();//_straw->id().getStation();                   // station number
	  phi     = polyAtan2(hit->pos().y() - center.y(),hit->pos().x() - center.x()); // atan2 returns its result in [-pi,pi], convert to [0,2pi]
	  if (phi < 0) phi += 2*M_PI;
	  zVec  [ist] += hit->pos().z();
	  //-----------------------------------------------------------------------------
	  // make sure there all hits within the station get close values of phi, although a 
	  // common 2pi ambiguity is still unresolved
	  //-----------------------------------------------------------------------------
	  if (nhits[ist] == 0) {
	    phiVec[ist] = phi;
	    nstations_with_hits += 1;
	  }
	  else {
	    while (phi-phiVec[ist] >  M_PI) phi -= 2*M_PI;
	    while (phi-phiVec[ist] < -M_PI) phi += 2*M_PI;
	
	    phiVec[ist] = (phiVec[ist]*nhits[ist]+phi)/(nhits[ist]+1);
	  }
	  nhits [ist] += 1;
	}
      }//end loop over the panels
    }

    for (int i=0; i<nstations; i++) {
      if (nhits[i] > 0) {
	zVec  [i] = zVec  [i]/nhits[i];
      }
    }

    if (_debug >5) {
      printf("[RobustHelixFinder::initFZ_2] StationID  nhits       z        phi\n");
      for (int i=0; i<nstations; i++) {
	if (nhits[i] > 0) printf("[RobustHelixFinder::initFZ_2] %5i %6i    %9.3f %8.5f\n", i,nhits[i],zVec[i],phiVec[i]);
      }
    }
//-----------------------------------------------------------------------------
// Step 2: determine the most likely value of phi
//-----------------------------------------------------------------------------
    for (int i=0; i<nstations; i++) {
      if (nhits[i] == 0)                                    continue; 

      phi_ref = phiVec[i];
      z_ref   = zVec  [i];

      for(int j=i+1; j<nstations; ++j) { 
	if (nhits[j] == 0)                                  continue;

	dphi = phiVec[j]-phi_ref;
	dz   = zVec[j] - z_ref;
	double dphidz = dphi/dz*double(rhel.helicity()._value);
	
	weight = nhits[i] + nhits[j];
//-----------------------------------------------------------------------------
// calculate N potential choices for 2*PI ambiguity resolution
//-----------------------------------------------------------------------------
	int n(0), nmax(0), nmin(0), nchoices = 0;

	double x = dphidz + n*2*M_PI/dz;

	if (x < _mindfdz) {
//-----------------------------------------------------------------------------
// for n=0, x < _mindfdz
//-----------------------------------------------------------------------------
	  while (x < _maxdfdz) {
	    if (x < _mindfdz) nmin = n+1;
	    nmax = n;
	    if ((x > _mindfdz) && (x < _maxdfdz)) nchoices += 1;
	    n += 1;
	    x += 2*M_PI/dz;
	  }
	}
	else if (x < _maxdfdz) {
//-----------------------------------------------------------------------------
// for n=0,   _xMin <= x < _xMax
//-----------------------------------------------------------------------------
	  while (x < _maxdfdz) {
	    nmax = n;
	    if ((x > _mindfdz) && (x < _maxdfdz)) nchoices += 1;
	    n += 1;
	    x += 2*M_PI/dz;
	  }

	  nmin = 0;
	  x    = dphidz+(nmin-1)*2*M_PI/dz;

	  while (x > _mindfdz) {
	    nchoices += 1;
	    nmin -= 1;
	    x    -= 2*M_PI/dz;
	  }
	}
	else {
//-----------------------------------------------------------------------------
// for n=0, x >= _xMax
//-----------------------------------------------------------------------------
	  while (x > _mindfdz) {
	    if (x > _maxdfdz) nmax = n-1;
	    nmin = n;
	    if ((x > _mindfdz) && (x < _maxdfdz)) nchoices += 1;
	    n -= 1;
	    x -= 2*M_PI/dz;
	  }
	}

	if (nchoices == 0)                                  continue;

	weight = 1.; // 1./nchoices;
//-----------------------------------------------------------------------------
// loop again over all choices and fill a histogram
// histogram is from 
//-----------------------------------------------------------------------------
	for (int n=nmin; n<=nmax; n++) { // 
	  double x = dphidz + n*2*M_PI/dz;
	  int bin = (x-minX)/stepX;
	  hist[bin] += weight;
	}
      }
    }
//-----------------------------------------------------------------------------
// the 'histogram' is filled, find a peak
//-----------------------------------------------------------------------------
    int ixmax = int(maxX/stepX);

    double swmax(0), sw, xmp(0);
    for (int ix=0; ix<ixmax-1; ix++) {
      sw = (hist[ix]+hist[ix+1]);
      if (sw > swmax) { 
	xmp = (stepX*(ix+0.5)*hist[ix] + stepX*(ix+1+0.5)*hist[ix+1])/sw;
	swmax = sw;
      }
    }
//-----------------------------------------------------------------------------
// Part 2: perform a more accurate estimate - straight line fit
//-----------------------------------------------------------------------------
    if (nstations_with_hits < 2) return false;//hdfdz = _mpDfDz;
    else                         hdfdz = xmp*rhel.helicity()._value;
//-----------------------------------------------------------------------------
// last step - determine phi0 = phi(z=0)
//-----------------------------------------------------------------------------
    double phi0(0), sdphi(0);
    int    sn(0);

    for (int i=0; i<nstations; i++) {
      if (nhits[i] == 0) continue;

      if (sn == 0) { // first station with hits gives the "2*PI normalization";
	phi0 = phiVec[i]-zVec[i]*hdfdz;
	sdphi = 0;
	sn    = 1;
      }
      else {
//-----------------------------------------------------------------------------
// for all points different from the first one need to choose the turn number
//-----------------------------------------------------------------------------
	dphi = phiVec[i]-(phi0+zVec[i]*hdfdz);
	double dphi_min = dphi;

	int n= 0;
	while (1) {
	  n += 1;
	  double dphi = phiVec[i]+2*M_PI*n-(phi0+zVec[i]*hdfdz);
	  if (fabs(dphi) < fabs(dphi_min)) dphi_min = dphi;
	  else break;
	}

	n=0;
	while (1) {
	  n -= 1;
	  double dphi = phiVec[i]+2*M_PI*n-(phi0+zVec[i]*hdfdz);
	  if (fabs(dphi) < fabs(dphi_min)) dphi_min = dphi;
	  else break;
	}
	
	sdphi += dphi_min;
	sn += 1;
      }
    }

    float hphi0 = phi0 + sdphi/sn;

    //now update the helix info
    rhel._fz0    = deltaPhi(0.0, hphi0);
    rhel._lambda = 1/hdfdz;
      
    if (_debug > 5) {
      printf("[RobustHelixFinder::initFZ_2] END: hdfdz = %9.5f hphi0 = %9.6f ", hdfdz, hphi0);
    }

    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facez     = &HelixData._oTracker[f];
      for (int p1=0; p1<FaceZ_t::kNPanels; ++p1){
	panelz  = &facez->panelZs[p1];   

	int nhitsPanelF1 = panelz->nChHits();
	if (nhitsPanelF1 == 0)                          continue;
      
	for (int ip=0; ip<nhitsPanelF1; ++ip){
	  ComboHit* hit = &HelixData._chHitsToProcess[panelz->idChBegin + ip];
	  // if (!use(*hitP1) )                          continue;   
	  resolvePhi(*hit,rhel);
	}
      }//end loop pver the panels
    }//end loop over the faces

    return goodFZ(rhel);
  }


  void RobustHelixFit::fitFZ(RobustHelixFinderData& HelixData) {
    // if required, initialize
    HelixData._hseed._status.clear(TrkFitFlag::phizOK);
    if (!HelixData._hseed._status.hasAllProperties(TrkFitFlag::phizInit))
      {
	if (initFZ(HelixData))
	  HelixData._hseed._status.merge(TrkFitFlag::phizInit);
	else
	  return;
      }

    // ComboHitCollection& hhits = HelixData._hseed._hhits;
    RobustHelix& rhel         = HelixData._hseed._helix;

    //iterate over lambda and loop resolution
    unsigned niter(0);
    bool changed(true);
    while(changed && niter < _maxniter)
      {
	changed = false;
	accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accf;

	ComboHit*      hitP1(0), *hitP2(0);
	uint16_t       facezF1(0), facezF2(0);
        int            nHits(HelixData._chHitsToProcess.size());

	for (int f1=0; f1<nHits-1; ++f1){
	  hitP1 = &HelixData._chHitsToProcess[f1];
	  if (!use(*hitP1))                               continue;
	  facezF1 = hitP1->strawId().uniqueFace();

	  for (int f2=f1+1; f2<nHits; ++f2){
	    hitP2 = &HelixData._chHitsToProcess[f2];
	    facezF2 = hitP2->strawId().uniqueFace();

	    if (!use(*hitP2) || (facezF1 == facezF2) )    continue;
		
	    float dz   = hitP2->pos().z() - hitP1->pos().z();//facezF2->z - facezF1->z;
	    float dphi = hitP2->helixPhi()- hitP1->helixPhi(); 
	    if (dz < _minzsep || fabs(dphi) < _mindphi) continue;

	    float lambda = dz/dphi;
	    if (goodLambda(rhel.helicity(),lambda)){
	      float wt = sqrtf(hitP1->nStrawHits()*hitP2->nStrawHits());
	      //		    float wt = (*ihit)->nStrawHits() + (*jhit)->nStrawHits();
	      accf(lambda, weight=wt);
	    }
		
	  }//end secondloop over faces 
	}//end first loop over faces

	rhel._lambda = extract_result<tag::weighted_median>(accf);

	// now extract intercept.  Here we solve for the difference WRT the previous value
	accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > acci;

	for (unsigned i=0; i<HelixData._chHitsToProcess.size(); ++i){ 
	  hitP1 = &HelixData._chHitsToProcess[i];
	  if (!use(*hitP1))   continue;
	  
	  float phiex = rhel.circleAzimuth(hitP1->_pos.z());
	  float dphi  = deltaPhi(phiex,hitP1->helixPhi());
	  float wt    = hitP1->nStrawHits();
	  acci(dphi,weight = wt);// accumulate the difference WRT the current intercept
	}


	// enforce convention on azimuth phase
	float dphi = extract_result<tag::weighted_median>(acci);
	rhel._fz0 = deltaPhi(0.0,rhel.fz0()+ dphi);

	// resolve the hit loops again
	for (unsigned i=0; i<HelixData._chHitsToProcess.size(); ++i){ 
	  hitP1    = &HelixData._chHitsToProcess[i];
	  changed |= resolvePhi(*hitP1,rhel);
	}

	++niter;
      }

    if (goodFZ(rhel)) {
      HelixData._hseed._status.merge(TrkFitFlag::phizOK);
      if (_refineZPhiFit) refineFitZPhi(HelixData);
    }
    if (_diag){
      HelixData._diag.lambda_1 = rhel._lambda;
    }
  }




  void RobustHelixFit::fitFZ_2(RobustHelixFinderData& HelixData, int UseInteligentWeights) {
    // if required, initialize
    HelixData._hseed._status.clear(TrkFitFlag::phizOK);
    if (!HelixData._hseed._status.hasAllProperties(TrkFitFlag::phizInit))
      {
    	if (initFZ(HelixData, 0))
    	  HelixData._hseed._status.merge(TrkFitFlag::phizInit);
    	else
    	  return;
      }

    RobustHelix& rhel         = HelixData._hseed._helix;

    if (goodFZ(rhel)) {
      HelixData._hseed._status.merge(TrkFitFlag::phizOK);
     
      //refine the circle fit
      //perform a reduced chi2 fit using only 1 ComboHit-per-face
      //the best hit is defined computing the residual
      if (_refineZPhiFit) refineFitZPhi(HelixData, UseInteligentWeights);
    }
    if (_diag){
      HelixData._diag.lambda_1 = rhel._lambda;
    }
  }

  //reduced linear fit to evaluate lambda and phi0
  void  RobustHelixFit::refineFitZPhi(RobustHelixFinderData& HelixData, int UseInteligentWeights){
    ::LsqSums4 szphi;
    float      phi,phi_ref,wt, dphi, z;
    int        minNReducedChi2Points(15);//FIXME!
    int        nZPHISh(0);

    RobustHelix& rhel   = HelixData._hseed._helix;
    XYVec        center(rhel.center().x(),rhel.center().y());
    float        radius = rhel.radius();
    float        dfdz   = 1/rhel._lambda;
    float        phi0   = rhel._fz0;

    ComboHit*    hitP1(0);
    FaceZ_t*     facezP1(0);
    int          nhitsFace1(0);
    
    float        dzMin_second_arch(500.), last_z_first_arch(9999.);
    int          nPoints_second_arch(0), nPoints_first_arch(0), minNPoints_second_arch(5);//OPTIMIZE ME!

    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facezP1    = &HelixData._oTracker[f];

      HitInfo_t    indexBestHit;

      nhitsFace1 = facezP1->nChHits();
      if (nhitsFace1 == 0)                          continue;

      float  xdphi(1e10), xdphiBestHit(_maxXDPhi);
      for (int ip=0; ip<nhitsFace1; ++ip){
	hitP1 = &HelixData._chHitsToProcess[facezP1->idChBegin + ip];
	z     = hitP1->pos().z();

	//set by default as outlier. Only the best hit will be used
	phi = hitP1->helixPhi();
	if (phi < 0) phi = phi + 2*M_PI;//      = TVector2::Phi_0_2pi(phi);
	//	dz  = z - zlast;
	
	phi_ref = z*dfdz + phi0;		// predicted value of phi
	dphi    = phi_ref - phi;		// signed residual
	// resolve 2PI ambiguity
	while (dphi > M_PI) {
	  phi += 2*M_PI;
	  dphi = phi_ref - phi;
	}
	while (dphi < -M_PI) {
	  phi -= 2*M_PI;
	  dphi = phi_ref - phi;
	}

	// store the corrected value of phi
	hitP1->_hphi = phi;

	if (!use(*hitP1) )                          continue;

	dphi        = fabs(dphi);
	wt          = 1./(_sigmaPhi*_sigmaPhi);
	if (UseInteligentWeights == 1) wt = evalWeightZPhi(*hitP1,center,radius);
	xdphi       = dphi*sqrtf(wt);


	//remove points with large residuals
	if (xdphi <= xdphiBestHit){
	  indexBestHit.face          = f;
	  indexBestHit.panel         = hitP1->strawId().uniquePanel();
	  indexBestHit.panelHitIndex = facezP1->idChBegin + ip;
	  xdphiBestHit               = xdphi;      
	  hitP1->_zphiWeight         = wt;
	}
	
	hitP1->_flag.merge(StrawHitFlag::outlier);

      }//end loop over the hits within a face
	    
      //now add the best hit if found!
      if (indexBestHit.face>=0) {
	hitP1    = &HelixData._chHitsToProcess[indexBestHit.panelHitIndex];

	//remove the outlier flag
	hitP1->_flag.clear(StrawHitFlag::outlier);

	//include the point
	szphi.addPoint(z,hitP1->_hphi,hitP1->_zphiWeight);

	//increase the StrwaHit counter
	nZPHISh += hitP1->nStrawHits();

	if (nPoints_first_arch == 0) {	//initialize last_z_first_arch
	  last_z_first_arch = z;
	  ++nPoints_first_arch;
	}else if ( (z - last_z_first_arch) >= dzMin_second_arch ){//we reached the second arch of the helix
	  ++nPoints_second_arch;
	}else  {
	  last_z_first_arch = z;
	  ++nPoints_first_arch;
	}	

	//if we collected enough points update the center and the radius
	if ( (szphi.qn() >= minNReducedChi2Points) &&
	     (nPoints_second_arch >= minNPoints_second_arch) &&
	     (fabs(dfdz - szphi.dfdz()) < 8e-4)){//in case delta-electron hits are present, with need to pay attention to the change of slope!
	  dfdz  = szphi.dfdz();
	  phi0  = szphi.phi0();
	  //	  zlast = z;
	}
      }
	    
    }//end loop over the faces

    //if we collected enough points update the results
    //should we also check the chi2?
    if (nZPHISh >= _minnsh){

      rhel._fz0    = szphi.phi0();
      rhel._lambda = 1/szphi.dfdz();
	
      HelixData._szphi   = szphi;
      HelixData._nZPhiSh = nZPHISh;

      if (_diag){
	if (UseInteligentWeights == 1){
	  HelixData._diag.lambdaszphi_1 = 1./szphi.dfdz();
	  HelixData._diag.chi2dszphi_1  = szphi.chi2DofLine();
	  HelixData._diag.nshszphi_1    = nZPHISh;
	} if (UseInteligentWeights == 0){
	  HelixData._diag.lambdaszphi_0 = 1./szphi.dfdz();
	  HelixData._diag.chi2dszphi_0  = szphi.chi2DofLine();
	  HelixData._diag.nshszphi_0    = nZPHISh;
	}
      }
    }

    
    
  }


  // simple median fit.  No initialization required
  void RobustHelixFit::fitCircleMedian(RobustHelixFinderData& HelixData) 
  {
    const float mind2 = _mindist*_mindist;
    const float maxd2 = _maxdist*_maxdist;
      
    // ComboHitCollection& hhits = HelixData._hseed._hhits;
    RobustHelix* rhel         = &HelixData._hseed._helix;
    accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, float > accx, accy, accr;

    // loop over all triples
    unsigned      ntriple(0);
 
    ComboHit*     hitP1(0), *hitP2(0), *hitP3(0);
    int           facezP1(-1), facezP2(-1), facezP3(-1);
    int           nHits(HelixData._chHitsToProcess.size());

    for (int f1=0; f1<nHits-2; ++f1){
      hitP1 = &HelixData._chHitsToProcess[f1];
      if (!use(*hitP1) )                          continue;
      XYWVec& wposP1 = HelixData._chHitsWPos[f1];
      facezP1 = wposP1.face();
      float   ri2    = wposP1.Mag2();

      for (int f2=f1+1; f2<nHits-1; ++f2){
	hitP2 = &HelixData._chHitsToProcess[f2];
	XYWVec& wposP2  = HelixData._chHitsWPos[f2];
	facezP2 = wposP2.face();
	if (!use(*hitP2) || (facezP1 == facezP2)) continue;

	float   dist2ij = (wposP1 - wposP2).Mag2();
	if (dist2ij < mind2 || dist2ij > maxd2) continue;	  

	float rj2 = wposP2.Mag2();

	for (int f3=f2+1; f3<nHits; ++f3){
	  hitP3 = &HelixData._chHitsToProcess[f3];
	  XYWVec& wposP3 = HelixData._chHitsWPos[f3];
	  facezP3 = wposP3.face();
	  if (!use(*hitP3) || (facezP2 == facezP3) ) continue;

	  float dist2ik = (wposP1-wposP3).Mag2();
	  float dist2jk = (wposP2-wposP3).Mag2();
	  if (dist2ik < mind2 || dist2jk < mind2 ||
	      dist2ik > maxd2 || dist2jk > maxd2)   continue;

	  // Heron's formula
	  float area2 = (dist2ij*dist2jk + dist2ik*dist2jk + dist2ij*dist2ik) - 0.5*(dist2ij*dist2ij + dist2jk*dist2jk + dist2ik*dist2ik);
	  if(area2 < _minarea2)              continue;
	  // this effectively measures the slope difference
	  float delta = (wposP3.x() - wposP2.x())*(wposP2.y() - wposP1.y()) -
	    (wposP2.x() - wposP1.x())*(wposP3.y() - wposP2.y());

	  float rk2 = wposP3.Mag2();

	  // find circle center for this triple
	  float cx = 0.5* (
			   (wposP3.y() - wposP2.y())*ri2 +
			   (wposP1.y() - wposP3.y())*rj2 +
			   (wposP2.y() - wposP1.y())*rk2 ) / delta;
	  float cy = -0.5* (
			    (wposP3.x() - wposP2.x())*ri2 +
			    (wposP1.x() - wposP3.x())*rj2 +
			    (wposP2.x() - wposP1.x())*rk2 ) / delta;
	  XYVec cent(cx,cy);
	  float rho = sqrtf((wposP1-cent).Mag2());
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

	      float wt = cbrtf(wposP1.weight()*wposP2.weight()*wposP3.weight()); 
	      accx(cx,weight = wt);
	      accy(cy,weight = wt);
	      if(_tripler) accr(rho,weight = wt);
	      if (ntriple>_ntripleMax) {
		f1=nHits-2; f2=nHits-1; f3=nHits;
	      }
	    }
	}//end loop for f3 Faces
      }//end loop for f2 Faces
    }//end loop for f1 Faces
    
    // median calculation needs a reasonable number of points to function
    if (ntriple > _ntripleMin)
      {        
	float centx = extract_result<tag::weighted_median>(accx);
        float centy = extract_result<tag::weighted_median>(accy);
	XYVec center(centx,centy);
	if(!_tripler) {
	  if(!_errrwt) {		   
	    for (unsigned i=0; i<HelixData._chHitsToProcess.size(); ++i){
	      hitP1 = &HelixData._chHitsToProcess[i];
	      if (!use(*hitP1) )                          continue;	  
	      float rho = sqrtf(pow(hitP1->pos().x() - centx,2) + pow(hitP1->pos().y() - centy,2));//(wp - center).Mag2());
	      accr(rho,weight = hitP1->nStrawHits()); 
	    }
	  } else {
	    // set weight according to the errors
	    for (unsigned i=0; i<HelixData._chHitsToProcess.size(); ++i){
	      hitP1 = &HelixData._chHitsToProcess[i];
	      if (!use(*hitP1) )                          continue;	  
	      XYVec rvec = (XYVec(hitP1->pos().x(),hitP1->pos().y())-center);
	      float wt   = evalWeightXY(*hitP1,center);
	      float rho  = sqrtf(rvec.Mag2());
	      // if (rho*wt > _xyHitCut)         continue;//FIXME! need an histogram to implement this cut
	      accr(rho,weight=wt);
	    }

	    if ( _usecc && HelixData._hseed.caloCluster().isNonnull()){
	    }
	  }
	}
	float rho = extract_result<tag::weighted_median>(accr);
        rhel->_rcent = sqrtf(center.Mag2());
        rhel->_fcent = polyAtan2(center.y(), center.x());//center.Phi();
        rhel->_radius = rho;

	//refine the circle fit
	//perform a chi2 fit using only 1 ComboHit-per-face
	//the best hit is defined computing the residual
	if (_refineXYFit){
	  //first perform the chi2 fit assuming all hits have same error
	  refineFitXY(HelixData,0);

	  //now that radius and center are more accurate, repeat the fit using the orientation of the hit
	  //to estimate more accurately the expected uncertanty
	  refineFitXY(HelixData);
	}

	// 	//check the number of strawHits associated with the track
	// 	if (nStrawHits < _minnhit) return;

        if (goodCircle(*rhel)) HelixData._hseed._status.merge(TrkFitFlag::circleOK);
      }
  
    if (_diag){
      if ( HelixData._diag.circleFitCounter == 0){
	HelixData._diag.ntriple_0 = ntriple;
	HelixData._diag.radius_0  = rhel->_radius;
      } 

      HelixData._diag.ntriple_1 = ntriple;
      HelixData._diag.radius_1  = rhel->_radius;
      
      ++HelixData._diag.circleFitCounter;
    }
  }
  
  void  RobustHelixFit::refineFitXY(RobustHelixFinderData&HelixData, int WeightMode){

    ::LsqSums4 sxy;
    if (_targetcon) {
      if (WeightMode == 1) 
//-------------------------------------------------------------------------------
// add stopping target center with a position error of 100 mm/sqrt(12) ~ 30mm => wt = 1/900
//-------------------------------------------------------------------------------
	sxy.addPoint(0.,0.,1./900.);
      else{
	float wtTarget=1.;
	sxy.addPoint(0.,0.,wtTarget);
      }
    }
	
    float      wt(1.), resid;
    int        minNReducedChi2Points(15);//FIXME!
    int        nXYSh(0);

    FaceZ_t*      facezP1(0);
    ComboHit*     hitP1(0);
    int           nhitsPanelF1(0);

    //get the radius
    RobustHelix* rhel  = &HelixData._hseed._helix;
    float        rho   = rhel->radius();
    XYVec        center(XYVec(rhel->centerx(), rhel->centery()));

    float        dzMin_second_arch(500.), last_z_first_arch(9999.);
    int          nPoints_second_arch(0), nPoints_first_arch(0), minNPoints_second_arch(5);//OPTIMIZE ME!

    for (int f=0; f<StrawId::_ntotalfaces; ++f){
      facezP1    = &HelixData._oTracker[f];
      HitInfo_t  indexBestComboHit;

      float      minResid(_minxyresid);//version 1:
      if (WeightMode == 0 ) minResid = _maxdxy;

      nhitsPanelF1 = facezP1->nChHits();
      if (nhitsPanelF1 == 0)                          continue;
	
      for (int ip=0; ip<nhitsPanelF1; ++ip){
	hitP1 = &HelixData._chHitsToProcess[facezP1->idChBegin + ip];
	if (!use(*hitP1) )                          continue;

	//set by default as outlier. Only the best hit will be used
	// hitP1->_flag.merge(StrawHitFlag::outlier);//FIXME!

	XYVec rvec = (XYVec(hitP1->pos().x(),hitP1->pos().y())-center);

	if (WeightMode == 1 ) wt         = evalWeightXY(*hitP1, center);

	resid      = fabs(sqrtf(rvec.Mag2()) - rho)*sqrtf(wt);
	if (resid < minResid) {
	  indexBestComboHit.face          = f;
	  indexBestComboHit.panel         = hitP1->strawId().uniquePanel();
	  indexBestComboHit.panelHitIndex = facezP1->idChBegin + ip;
	    
	  minResid          = resid;
	  hitP1->_xyWeight  = wt;
	}
	
	hitP1->_flag.merge(StrawHitFlag::outlier);
      }//end loop over the hits in a Face
	    
      //now add the best hit if found!
      if (indexBestComboHit.face >=0 ) {
	hitP1    = &HelixData._chHitsToProcess[indexBestComboHit.panelHitIndex];

	//remove the outlier flag
	hitP1->_flag.clear(StrawHitFlag::outlier);//FIXME!

	//include the point
	sxy.addPoint(hitP1->pos().x(),hitP1->pos().y(),hitP1->_xyWeight);

	//increase the StrwaHit counter
	nXYSh += hitP1->nStrawHits();

	if (nPoints_first_arch == 0) {	//initialize last_z_first_arch
	  last_z_first_arch = hitP1->pos().z();//facezP1->z;
	  ++nPoints_first_arch;
	}else if ( (hitP1->pos().z()/*facezP1->z*/ - last_z_first_arch) >= dzMin_second_arch ){//we reached the second arch of the helix
	  ++nPoints_second_arch;
	}else  {
	  last_z_first_arch = hitP1->pos().z();//facezP1->z;
	  ++nPoints_first_arch;
	}	

	//if we collected enough points update the center and the radius
	if ( (sxy.qn() >= minNReducedChi2Points) &&
	     (nPoints_second_arch >= minNPoints_second_arch) ){
	  center.SetX(sxy.x0());
	  center.SetY(sxy.y0());
	  rho    = sxy.radius();
	}
      }
	    
    }//end loop over the faces

    //if we collected enough points update the results
    //should we also check the chi2?
    if (nXYSh >= _minnsh){
      //update the parameters
      center.SetX(sxy.x0());
      center.SetY(sxy.y0());
      rho    = sxy.radius();
      
      rhel->_rcent = sqrtf(center.Mag2());
      rhel->_fcent = polyAtan2(center.y(), center.x());//center.Phi();
      rhel->_radius = rho;
	  
      HelixData._sxy   = sxy;
      HelixData._nXYSh = nXYSh;
      if (_diag){
	if ( HelixData._diag.circleFitCounter == 0){
	  HelixData._diag.rsxy_0     = rho;
	  HelixData._diag.chi2dsxy_0 = sxy.chi2DofCircle();
	  HelixData._diag.nshsxy_0   = nXYSh;
	}else {
	  HelixData._diag.rsxy_1     = rho;
	  HelixData._diag.chi2dsxy_1 = sxy.chi2DofCircle();
	  HelixData._diag.nshsxy_1   = nXYSh;

	}
      }
    }
  }



  void RobustHelixFit::findAGE(RobustHelixFinderData const& HelixData, XYZVec const& center,float& rmed, float& age)
  {     
    const ComboHitCollection& hhits = HelixData._hseed._hhits;

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
    if(_usecc && HelixData._hseed.caloCluster().isNonnull())
      {
	// mu2e::GeomHandle<mu2e::Calorimeter> ch;
	// const Calorimeter* calo = ch.get();
	XYZVec cog = Geom::toXYZVec(_calorimeter->geomUtil().mu2eToTracker(_calorimeter->geomUtil().diskFFToMu2e(HelixData._hseed.caloCluster()->diskId(),HelixData._hseed.caloCluster()->cog3Vector())));
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


  void RobustHelixFit::fillSums(RobustHelixFinderData const& HelixData, XYZVec const& center,float rmed, AGESums& sums)
  {    
    ComboHitCollection const& hhits = HelixData._hseed._hhits;
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
    float dx = hhit._pos.x() - rhel.center().x();
    float dy = hhit._pos.y() - rhel.center().y();
    hhit._hphi = polyAtan2(dy,dx);//XYZVec(hhit._pos - rhel.center()).phi();
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

  bool RobustHelixFit::goodZPhiFit(LsqSums4& lsqsum)
  {
    return fabs(1./lsqsum.dfdz()) > _lmin && fabs(1./lsqsum.dfdz()) < _lmax;
  }

  bool RobustHelixFit::goodCircle(const RobustHelix& rhel)
  {
    return rhel.radius() > _rmin && rhel.radius() < _rmax;
  }

  bool RobustHelixFit::goodHelix(const RobustHelix& rhel)
  {
    return goodCircle(rhel) && goodFZ(rhel);
  }
  
  bool RobustHelixFit::goodHelixChi2(RobustHelixFinderData& helixData)
  {
    RobustHelix* rhel = &helixData._hseed._helix;
    return goodCircle(*rhel) && goodFZ(*rhel) && (helixData._sxy.chi2DofCircle() <= _chi2xymax) && (helixData._szphi.chi2DofLine() <= _chi2zphimax);
  }
  
  void RobustHelixFit::defineHelixParams(RobustHelixFinderData& helixData)
  {
    //this function is needed in case we used refineFitZPhi to set the value of phi0 properly
    static const float pi(M_PI), halfpi(pi/2.0);

    float amsign =  helixData._hseed._helix._helicity == Helicity::poshel ? 1. : -1.;
    
    float  dx     = amsign*helixData._hseed._helix.center().x();
    float  dy     = -amsign*helixData._hseed._helix.center().y();
    float  fcent  = polyAtan2(dy, dx);
    float  dphi   = deltaPhi(helixData._hseed._helix._fz0+amsign*halfpi,fcent);

       helixData._hseed._helix._fz0 = fcent - halfpi*amsign - dphi;
    //    helixData._hseed._helix._fz0 = deltaPhi(0.0, helixData._hseed._helix._fz0);
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

  float RobustHelixFit::evalWeightXY(const ComboHit& Hit, XYVec& Center){
    // XYVec rvec = (XYVec(Hit.pos().x(),Hit.pos().y())-Center);
    // XYVec rdir = rvec.unit();
    // float wdot = rdir.Dot(XYVec(Hit.wdir().x(),Hit.wdir().y()));
    // float wdot2 = wdot*wdot;
    // float tdot2 = 1.0 - wdot2;
    // float err2 = wdot2*Hit.wireErr2() + tdot2*Hit.transErr2();
    // float wt = 1/err2;// or 1.0/sqrtf(err2); // or 1/err2?

    float    transErr = 5./sqrt(12.);
    //scale the error based on the number of the strawHits that are within teh ComboHit
    if (Hit.nStrawHits() > 1) transErr *= 1.5;
    float    transErr2 = transErr*transErr;

    static const XYZVec _zdir(0.0,0.0,1.0);
    XYZVec _sdir  = _zdir.Cross(Hit._wdir);
   
    float x   = Hit.pos().x();
    float y   = Hit.pos().y();
    float dx  = x-Center.x();
    float dy  = y-Center.y();
    float dxn = dx*_sdir.x()+dy*_sdir.y();

    float costh2 = dxn*dxn/(dx*dx+dy*dy);
    float sinth2 = 1-costh2;

    // float e2     = _ew*_ew*sinth2+rs*rs*costh2;
    float e2     = Hit.wireErr2()*sinth2+transErr2*costh2;
    float wt     = 1./e2;
      
    return wt;
  }

  float RobustHelixFit::evalWeightZPhi(const ComboHit& Hit, XYVec& Center, float Radius){
    float    transErr = 5./sqrt(12.);
    //scale the error based on the number of the strawHits that are within teh ComboHit
    if (Hit.nStrawHits() > 1) transErr *= 1.5;
    float    transErr2 = transErr*transErr;

    float x  = Hit.pos().x();
    float y  = Hit.pos().y();
    float dx = x-Center.x();
    float dy = y-Center.y();
    
    static const XYZVec _zdir(0.0,0.0,1.0);
    XYZVec _sdir  = _zdir.Cross(Hit._wdir);
    
    float dxn    = dx*_sdir.x()+dy*_sdir.y();

    float costh2 = dxn*dxn/(dx*dx+dy*dy);
    float sinth2 = 1-costh2;

    // float e2     = _ew*_ew*costh2+rs*rs*sinth2;
    float e2     = Hit.wireErr2()*costh2+transErr2*sinth2;
    float wt     = Radius*Radius/e2;
    //    wt           *= 0.025;//_weightZPhi;

    return wt;
  }

}
