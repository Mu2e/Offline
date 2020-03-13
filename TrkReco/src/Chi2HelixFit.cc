//
// Object to perform helix fit to straw hits
//
// $Id: $
// $Author: $
// $Date: $
//
// mu2e
#include "TrkReco/inc/Chi2HelixFit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TrackerGeom/inc/Tracker.hh"

#include "RecoDataProducts/inc/CaloCluster.hh"

#include "Mu2eUtilities/inc/polyAtan2.hh"

// root
// #include "TH1F.h"
// #include "Math/VectorUtil.h"
// #include "Math/Vector2D.h"
//c++
#include <vector>
#include <utility>
#include <string>
#include <cmath>

using namespace std;
using namespace ROOT::Math::VectorUtil;

namespace mu2e
{
 
 
  typedef std::pair<float,float> WVal;

  Chi2HelixFit::Chi2HelixFit(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _dontuseflag(pset.get<std::vector<std::string>>("DontUseFlag",vector<string>{"Outlier"})),
    _minnsh(pset.get<unsigned>("minNStrawHits",10)),
    _minxyresid(pset.get<float>("minXYResid",5.)),
    _chi2xymax(pset.get<float>("chi2xymax",5.)),
    _chi2zphimax(pset.get<float>("chi2zphimax",5.)),
    _mindfdz(pset.get<float>("minDfDz",0.002)),
    _maxdfdz(pset.get<float>("maxDfDz",0.01)),
    _sigmaPhi(pset.get<float>("sigmaPhi",0.49)),//1636)),// rad
    _maxdxy(pset.get<float>("maxdxy",100.0)),
    _maxXDPhi(pset.get<float>("maxXDPhi",5.0))
  { }

  Chi2HelixFit::~Chi2HelixFit()
  {}

  bool Chi2HelixFit::initChi2Circle(RobustHelixFinderData& HelixData, bool TargetCon) {
    bool retval(false);

    fitChi2CircleMedian(HelixData, TargetCon);
    retval = HelixData._hseed._status.hasAllProperties(TrkFitFlag::circleOK);
    return retval;
  }

  void Chi2HelixFit::fitChi2Circle(RobustHelixFinderData& HelixData, bool TargetCon) {
    HelixData._hseed._status.clear(TrkFitFlag::circleOK);

    // if required, initialize
    bool init(false);
    if (!HelixData._hseed._status.hasAllProperties(TrkFitFlag::circleInit)) {
      init = true;
      if (initChi2Circle(HelixData, TargetCon))
	HelixData._hseed._status.merge(TrkFitFlag::circleInit);
      else
	return;
    }
    if (!init) fitChi2CircleMedian(HelixData, TargetCon);
  }

//----------------------------------------------------------------------------------------
// 2015-01-13  calculate track DphiDz using histogrammed distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  bool Chi2HelixFit::initFZ_2(RobustHelixFinderData& HelixData) {

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
    
    nstations = StrawId::_nstations;//RobustHelixFinderData::kNStations;

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

    return _rhfit->goodFZ(rhel);
  }


  void Chi2HelixFit::initFitChi2FZ(RobustHelixFinderData& HelixData) {
    _rhfit->fitFZ(HelixData);
    
    if (_rhfit->goodFZ(HelixData._hseed._helix)) {
      //      HelixData._hseed._status.merge(TrkFitFlag::phizOK);
      refineFitZPhi(HelixData);
    }
    if (_diag){
      HelixData._diag.lambda_1 = HelixData._hseed._helix._lambda;
    }
  }




  void Chi2HelixFit::fitChi2FZ(RobustHelixFinderData& HelixData, int UseInteligentWeights) {
    // if required, initialize
    HelixData._hseed._status.clear(TrkFitFlag::phizOK);
    if (!HelixData._hseed._status.hasAllProperties(TrkFitFlag::phizInit))
      {
    	if (_rhfit->initFZ(HelixData, 0))
    	  HelixData._hseed._status.merge(TrkFitFlag::phizInit);
    	else
    	  return;
      }

    RobustHelix& rhel         = HelixData._hseed._helix;

    if (_rhfit->goodFZ(rhel)) {
      HelixData._hseed._status.merge(TrkFitFlag::phizOK);
     
      //refine the circle fit
      //perform a reduced chi2 fit using only 1 ComboHit-per-face
      //the best hit is defined computing the residual
      refineFitZPhi(HelixData, UseInteligentWeights);
    }
    if (_diag){
      HelixData._diag.lambda_1 = rhel._lambda;
    }
  }

  //reduced linear fit to evaluate lambda and phi0
  void  Chi2HelixFit::refineFitZPhi(RobustHelixFinderData& HelixData, int UseInteligentWeights){
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
  void Chi2HelixFit::fitChi2CircleMedian(RobustHelixFinderData& HelixData, bool TargetCon) 
  {
    _rhfit->fitCircleMedian(HelixData, false);
    HelixData._hseed._status.clear(TrkFitFlag::circleOK);
      
    //first perform the chi2 fit assuming all hits have same error
    refineFitXY(HelixData, TargetCon, 0);
    
    //now that radius and center are more accurate, repeat the fit using the orientation of the hit
    //to estimate more accurately the expected uncertanty
    refineFitXY(HelixData, TargetCon);
    
    RobustHelix* rhel = &HelixData._hseed._helix;

    if (_rhfit->goodCircle(*rhel)) HelixData._hseed._status.merge(TrkFitFlag::circleOK);
      
    
    if (_diag){
      if ( HelixData._diag.circleFitCounter == 1){
	HelixData._diag.radius_0  = rhel->_radius;
      } 

      HelixData._diag.radius_1  = rhel->_radius;
      
      ++HelixData._diag.circleFitCounter;
    }
  }
  
  void  Chi2HelixFit::refineFitXY(RobustHelixFinderData&HelixData, bool Targetcon, int WeightMode){

    ::LsqSums4 sxy;
    if (Targetcon) {
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
      if (Targetcon) {
	HelixData._nXYCh = sxy.qn() - 1;
      }else {
	HelixData._nXYCh = sxy.qn();      
      }

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

  float Chi2HelixFit::deltaPhi(float phi1, float phi2)
  {
    float dphi = fmod(phi2-phi1,CLHEP::twopi);
    if (dphi>CLHEP::pi) dphi -= CLHEP::twopi;
    if (dphi<-CLHEP::pi)dphi += CLHEP::twopi;
    return dphi;
  }

  bool Chi2HelixFit::use(const ComboHit& hhit) const 
  {
    return (!hhit._flag.hasAnyProperty(_dontuseflag));
  }

  bool Chi2HelixFit::resolvePhi(ComboHit& hhit, const RobustHelix& rhel) const
  {
    // find phi expected
    float phiex = rhel.circleAzimuth(hhit._pos.z());
    int nloop = (int)rint((phiex-hhit.helixPhi())/CLHEP::twopi);
    hhit._hphi += nloop*CLHEP::twopi;

    static StrawHitFlag resphi(StrawHitFlag::resolvedphi);
    hhit._flag.merge(resphi);
    return nloop != 0;
  }
  
  bool Chi2HelixFit::goodZPhiFit(LsqSums4& lsqsum)
  {
    return fabs(1./lsqsum.dfdz()) > _rhfit->lambdaMin() && fabs(1./lsqsum.dfdz()) < _rhfit->lambdaMax();
  }

  bool Chi2HelixFit::goodHelixChi2(RobustHelixFinderData& helixData)
  {
    RobustHelix* rhel = &helixData._hseed._helix;
    return _rhfit->goodCircle(*rhel) && _rhfit->goodFZ(*rhel) && (helixData._sxy.chi2DofCircle() <= _chi2xymax) && (helixData._szphi.chi2DofLine() <= _chi2zphimax);
  }
  
  void Chi2HelixFit::defineHelixParams(RobustHelixFinderData& helixData)
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

  float Chi2HelixFit::hitWeight(const ComboHit& hhit) const 
  {
    float retval(hhit.nStrawHits());
    // add an option to evaluate the error relative to the current center FIXME
    return retval;
  }

  float Chi2HelixFit::evalWeightXY(const ComboHit& Hit, XYVec& Center){
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

  float Chi2HelixFit::evalWeightZPhi(const ComboHit& Hit, XYVec& Center, float Radius){
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
