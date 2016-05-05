///////////////////////////////////////////////////////////////////////////////
// helix fit to straw hits
//
// $Id: HelixFitHack.cc,v 1.13 2014/06/06 21:35:08 murat Exp $
// $Author: murat $ 
// $Date: 2014/06/06 21:35:08 $
//
//  use of HackData:
//  ----------------
//  [00:13] : Giani
//  [14:15] : Pasha - parameters of fit with non-equal weights
//
// the following has to come before other BaBar includes
///////////////////////////////////////////////////////////////////////////////
#include "BTrk/BaBar/BaBar.hh"
#include "CalPatRec/inc/HelixFitHack.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// framework
#include "fhiclcpp/ParameterSet.h"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
// Root
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TArc.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TFile.h"
#include "TVector2.h"
#include "TMarker.h"
#include "TCanvas.h"
// C++
#include <vector>
#include <string>
#include <algorithm>

#include "CalPatRec/inc/HelixFitHack.hh"
#include "CalPatRec/inc/THackData.hh"

using CLHEP::HepVector;
using CLHEP::Hep3Vector;

namespace mu2e {

//   // comparison functor for ordering points
//   struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
//     bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
//   };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<XYZPHack,XYZPHack,bool> {
    bool operator()(XYZPHack const& p1, XYZPHack const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

//-----------------------------------------------------------------------------
  void HelixFitHack::helixParams(HelixFitHackResult const& helix,
				 CLHEP::HepVector&         pvec,
				 CLHEP::HepVector&         perr) const 
  {
    HelixDefHack const& mytrk = helix._hdef;
    static const double pi(M_PI);
    //    static const double twopi(2*pi);
    static const double halfpi(pi/2.0);
    // the helix fit introduces a radial bias due to an asymmetry in the detector (more phase space for
    // noise hits outside the circle than inside.  correct for it.
    double radius = helix._radius + _rbias;
    pvec = HepVector(5,0);
    // omega is the inverse transverse radius of the particle's circular motion.  Its
    // signed by the particle angular momentum about the cirle center.
    // This CANNOT be deduced geometrically, so must be supplied as an ad-hoc assumption
    double amsign = copysign(1.0,-mytrk.particle().charge()*bz());
    pvec[HelixTraj::omegaIndex] = amsign/radius;
    // phi0 is the azimuthal angle of the particle velocity vector at the point
    // of closest approach to the origin.  It's sign also depends on the angular
    // momentum.  To translate from the center, we need to reverse coordinates
    pvec[HelixTraj::phi0Index] = atan2(-amsign*helix._center.x(),amsign*helix._center.y());
    // d0 describes the distance to the origin at closest approach.
    // It is signed by the particle angular momentum WRT the origin.
    // The Helix fit radial bias is anti-correlated with d0; correct for it here.
    pvec[HelixTraj::d0Index] = amsign*(helix._center.perp() - helix._radius - 2*_rbias);
    // the dip angle is measured WRT the perpendicular.  It is signed by the particle Z momentum    
    pvec[HelixTraj::tanDipIndex] = amsign/(radius*helix._dfdz);
    // must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
    // resolve the loop ambiguity such that the POCA is closest to z=0.
    double dphi = deltaPhi(helix._fz0+amsign*halfpi,pvec[HelixTraj::phi0Index]);
    // choose z0 (which loop) so that f=0 is as close to z=0 as possible
    pvec[HelixTraj::z0Index] = dphi*pvec[HelixTraj::tanDipIndex]/pvec[HelixTraj::omegaIndex];
    // estimated covariance based on average performance.  These should be parameters, FIXME!!!
    perr = HepVector(5,0);
    perr[HelixTraj::d0Index]     = 34.0;
    perr[HelixTraj::phi0Index]   = 0.02;
    perr[HelixTraj::omegaIndex]  = 0.0002;
    perr[HelixTraj::tanDipIndex] = 0.05;
    perr[HelixTraj::z0Index]     = 15.0;

  }

//-------------------------------------------------------------------------//
//  2014-12-26 Gianipez added the following method for asking t the helix
// finder if the hit "index" has already been marked as used by a previous 
// by a previous search
//-----------------------------------------------------------------------------
  int   HelixFitHack::isHitUsed(int index) {
    if ( (_goodPointsTrkCandidate < _minPointsTrkCandidate) ||
	 //(_chi2TrkCandidate > _maxChi2TrkCandidate) ||
	 _markCandidateHits == 0) return 0;
    
    if (index >= 400) {
      printf("[HelixFitHack::isHitUsed] requested index = %i range exceeded the range allowed\n", index);
      return 1;
    }

    int rc = (_indicesTrkCandidate[index] > 0);
    return rc;
  }

//-----------------------------------------------------------------------------
  HelixFitHack::HelixFitHack(fhicl::ParameterSet const& pset) :
    fHitChi2Max(pset.get<double>("hitChi2Max"       )),
    _diag      (pset.get<int>   ("diagLevel"        )),
    _debug     (pset.get<int>   ("debugLevel"       )),
    _debug2    (pset.get<int>   ("debugLevel2"      )),
    _minnhit   (pset.get<int>   ("minNHit"          )),
    _minzsep   (pset.get<double>("minzsep",50.0)),
    _maxzsep   (pset.get<double>("maxzsep",500.0)),
    _maxdz     (pset.get<double>("maxdz",35.0)),
    _maxdot    (pset.get<double>("maxdot",0.9)),
    _mpDfDz    (pset.get<double>("mostProbableDfDz")),
    _maxDfDz   (pset.get<double>("maxDfDz",0.01)),
    _minDfDz   (pset.get<double>("minDfDz",5e-04)),
    _sigmaPhi  (pset.get<double>("sigmaPhi")),
    _weightXY  (pset.get<double>("weightXY")),
    _weightZPhi(pset.get<double>("weightZPhi")),
    _maxXDPhi  (pset.get<double>("maxXDPhi",5.)),
    _distPatRec(pset.get<double>("distPatRec")),
    _rbias     (pset.get<double>("radialBias",0.0)),
    _efac      (pset.get<double>("ErrorFactor",1.0)),
    _rhomin    (pset.get<double>("rhomin",350.0)),
    _rhomax    (pset.get<double>("rhomax",780.0)),
    _mindist   (pset.get<double>("mindist",50.0)),
    _maxdist   (pset.get<double>("maxdist",500.0)),
    _pmin      (pset.get<double>("minP",50.0)),
    _pmax      (pset.get<double>("maxP",150.0)),
    _tdmin     (pset.get<double>("minAbsTanDip",0.3)),
    _tdmax     (pset.get<double>("maxAbsTanDip",2.0)),
    _rcmin     (pset.get<double>("rcmin",200.0)),
    _rcmax     (pset.get<double>("rcmax",350.0)),
    _forcep    (pset.get<bool>  ("forceP",false)),
    _xyweights (pset.get<bool>  ("xyWeights",false)),
    _zweights  (pset.get<bool>  ("zWeights",false)),
    _filter    (pset.get<bool>  ("filter",true)),
    _plotall   (pset.get<bool>  ("plotall",false)),
    _usetarget (pset.get<bool>  ("usetarget",true)),
    _bz        (0.0),
    _x0        (-9999.),
    _y0        (-9999.), 
    _phi0      (-9999.), 
    _radius    (-9999.), 
    _dfdz      (-9999.),
    _goodPointsTrkCandidate(-9999),
    _minPointsTrkCandidate (pset.get<int>("minPointsTrkCandidate")),
    _chi2TrkCandidate      (1e10),
    _maxChi2TrkCandidate   (pset.get<double>("maxChi2TrkCandidate")),
    _markCandidateHits     (pset.get<int>("markCandidateHits")),
    _chi2xyMax             (pset.get<double>("chi2xyMax")),
    _chi2zphiMax           (pset.get<double>("chi2zphiMax")),
    _dfdzErr               (pset.get<double>("dfdzErr")){

    XYZPHack::_efac = _efac;
    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    XYZPHack::_useflag = StrawHitFlag(bitnames);

    for (int i=0; i<400; ++i){
      _indicesTrkCandidate[i] = -9999;
      _distTrkCandidate[i]    = -9999;
      _dzTrkCandidate[i]      = -9999;
    }

    _chi2nFindZ    = 0.0,
    _eventToLook   = -1;
    
    _hDfDzRes      = NULL;//new TH1F("hDfDzRes","dfdz residuals",
    //  20, _minDfDz, _maxDfDz);
  }


//-----------------------------------------------------------------------------
  HelixFitHack::~HelixFitHack() {
    if (_hDfDzRes){
      delete _hDfDzRes;
      _hDfDzRes = 0;
    }
  }

//-----------------------------------------------------------------------------
  double HelixFitHack::bz() const {
    if(_bz == 0.0){
      // find the magnetic field Z component at the origin
      GeomHandle<BFieldManager> bfmgr;
      GeomHandle<DetectorSystem> det;
      // change coordinates to mu2e
      CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
      CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
      CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
      _bz = field.z();
    }
    return _bz;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  bool HelixFitHack::findHelix(HelixFitHackResult& Helix, const CalTimePeak* TimePeak) {

    fTimePeak = TimePeak;

    HelixDefHack const& mytrk = Helix._hdef;
//-----------------------------------------------------------------------------
//  compute the allowed radial range for this fit
//-----------------------------------------------------------------------------
    double pb = fabs((CLHEP::c_light*1e-3)/(bz()*mytrk.particle().charge()));
    _rmin = _pmin/(pb*sqrt(1.0+_tdmax*_tdmax));
    _rmax = _pmax/(pb*sqrt(1.0+_tdmin*_tdmin));
//-----------------------------------------------------------------------------
//  particle charge, field, and direction affect the pitch range
//-----------------------------------------------------------------------------
    _dfdzsign = copysign(1.0,-mytrk.particle().charge()*mytrk.fitdir().dzdt()*bz());

    if(_dfdzsign > 0.0){
      _smin = 1.0/(_rmax*_tdmax);
      _smax = 1.0/(_rmin*_tdmin);
    } else {
      _smax = -1.0/(_rmax*_tdmax);
      _smin = -1.0/(_rmin*_tdmin);
    }
//-----------------------------------------------------------------------------
// initialize internal array of hits, print if requested
//-----------------------------------------------------------------------------
    fillXYZP(mytrk);
				        // 2013-10-18 P.Murat: print 
    int n = _xyzp.size();

    XYZPHack* pt;
    if (_debug > 0) {
      printf("[HelixFitHack::findHelix] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g\n",
	     0.,0.,0.,-1.);
    }
    int banner_printed(0);
    for(int i=0; i<n; ++i) {
      pt = &_xyzp[i];

      if (_debug > 0) {
	if (banner_printed == 0) {
	  printf("[HelixFitHack::findHelix] flag    Used  X       Y       Z      _debug: %i5\n",_debug);
	  printf("[HelixFitHack::findHelix] -------------------------------------\n");
	  banner_printed = 1;
	}
	printf("[HelixFitHack::findHelix] %08x %2i %12.5f %12.5f %12.5f \n",
	       *((int*) &pt->_flag), pt->use(), pt->_pos.x(), pt->_pos.y(), pt->_pos.z()
	       );
      }
    }
//-----------------------------------------------------------------------------
// call down
//-----------------------------------------------------------------------------
    bool retval = findHelix(Helix);

    return retval;
  }
//---------------------------------------------------------------------------
// reset track paramters
//---------------------------------------------------------------------------
  void HelixFitHack::resetTrackParamters() {
// indices on the xyzp vector of: the straw hit seeding the search, the second strawhit used for recalculating the dfdz value
// fLastIndex which is a helpful paramter for searching goodpoints during the pattern recognition
    fSeedIndex      = -9999;   
    fCandidateIndex = -9999;
    fLastIndex      = -9999;   
    fUseDefaultDfDz = 0;
					// follow helix paramters
    _x0     = -9999.;
    _y0     = -9999.;
    _phi0   = -9999.;
    _radius = -9999.;
    _dfdz   = -9999.;
//-----------------------------------------------------------------------------
// quality paramters used for doing comparison between several track candidates
//-----------------------------------------------------------------------------
    _goodPointsTrkCandidate = -99999;
    _chi2TrkCandidate       = 1e10;
//-----------------------------------------------------------------------------
// vector which holds indices of the strawhit far from the predicted position
//-----------------------------------------------------------------------------
    for(int i=0; i<400; ++i){
      _indicesTrkCandidate[i] = -9999;
      _distTrkCandidate   [i] = -9999;
      _dzTrkCandidate     [i] = -9999;
      fPhiCorrected       [i] = -9999;
    }
  }

//-----------------------------------------------------------------------------
// called internally; in the diagnostics mode save several states of _xyzp
//-----------------------------------------------------------------------------
  bool HelixFitHack::findHelix(HelixFitHackResult& Helix) {
    bool retval(false);

//-----------------------------------------------------------------------------
// 2014-11-09 gianipez: reset the track candidate parameters if a new time peak is used! 
// so the previous candidate should not be compared to the new one at this level
//-----------------------------------------------------------------------------
    resetTrackParamters();
//-----------------------------------------------------------------------------
// save results in the very beginning
//-----------------------------------------------------------------------------
    if (_diag > 0) saveResults(_xyzp,Helix,0);

    if (_filter) {
      filterDist();
      if (_diag > 0) saveResults(_xyzp,Helix,1);
    }

    doPatternRecognition(Helix);
//---------------------------------------------------------------------------
// 2014-11-11 gianipez changed the following if() statement to test the
// possibility of spead up the pattern recognition in presence of background
// how the number of good points may be different from the number used in sums ?
//---------------------------------------------------------------------------
    if (_debug != 0) {
      printf("[HelixFitHack::findHelix] Helix._sxy.qn() = %5.0f goodPointsTrkCandidate = %i\n", 
	     Helix._sxy.qn(), _goodPointsTrkCandidate);
    }
    
    if (_goodPointsTrkCandidate< int(_minnhit) ) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,1); // small number of hits
    }
    else if ((Helix._radius < _rmin) || (Helix._radius > _rmax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,2); // initialization failure
    }
    else if ( (Helix._sxyw.qn() < int(_minnhit)) || (Helix._sxyw.chi2DofCircle() > _chi2xyMax)){
      Helix._fit = TrkErrCode(TrkErrCode::fail,3); // xy reconstruction failure
    }
    else if ( (Helix._srphi.qn() < int(_minnhit)) || (Helix._srphi.chi2DofLine() > _chi2zphiMax) ){
      Helix._fit = TrkErrCode(TrkErrCode::fail,4); // phi-z reconstruction failure
    }else {
				// success
      retval = true;
    }
    return retval;
  }

//----------------------------------------------------------------------------------------
// 2015-01-13  calculate the dfdz of the track using distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  int HelixFitHack::findDfDz(HelixFitHackResult& Helix, 
			      int seedIndex, int *indexVec) {
    
    double phi, phi_ref(-1e10), z, z_ref, dphi, dz, dzOverHelPitch;
    
    _hDfDzRes = new TH1F("hDfDzRes","dfdz residuals",
		      20, _minDfDz, _maxDfDz);
    _hPhi0Res = new TH1F("hPhi0Res", "phi0 residuals",
			 20, 0., 2.*M_PI);
    mu2e::GeomHandle<mu2e::TTracker> ttHandle;
    const mu2e::TTracker* tracker = ttHandle.get();

    CLHEP::Hep3Vector center = Helix._center;
    CLHEP::Hep3Vector pos_ref, pos;

    //2015 - 03 -30 G. Pezzu changed the value of tollMax.
    // using the initial value of dfdz we can set it more accuratelly:
    // tollMax = half-helix-step = Pi / dfdz
    double tollMin(100.), tollMax;//(800.);
    tollMax = 2.*M_PI / Helix._dfdz;
    
    //avoid the use of the couple of points which have dz %(mod) HelPitch less than 0.8
    // still the value need to be optimized, this is just out of lots of debugging
    double dzOverHelPitchCut = 0.7;

    if (_debug >5){
      printf("[HelixFitHack::findDfDz] x0 = %9.3f y0 = %9.3f radius = %9.3f dfdz = %9.6f straw-hits = %9.5f dzPitch = %8.6f\n",
	     center.x(), center.y(), Helix._radius, Helix._dfdz, (Helix._sxy.qn() - 1), tollMax);// -1 to remove the EMC cluster contribute
      //      printInfo(Helix);
    }

    int        np, ist, nstations;

    np        = _xyzp.size();
    nstations = tracker->nStations();
    
    double    phiVec[30], zVec[30];
    int       nhits [30];

    for (int i=0; i<nstations; i++) {
      phiVec[i] = 0;
      zVec  [i] = 0;
      nhits [i] = 0;
    }
//-----------------------------------------------------------------------------    
// Part 1: use only contiguous parts of the trajectory
//-----------------------------------------------------------------------------
    for (int i=seedIndex; i<np; i++) {
      if ( (! _xyzp[i].isOutlier()) && (indexVec[i] >0 )) {
	// didn't find an accessor returning the station number, hack
	ist = _xyzp[i]._straw->id().getPlane()/2;
	pos = _xyzp[i]._pos;
	phi = CLHEP::Hep3Vector(pos - center).phi();
	phi = TVector2::Phi_0_2pi(phi);
	zVec  [ist] += pos.z();
	if (nhits == 0) phiVec[ist] = phi;
	else {
	  while (phi-phiVec[ist] >  M_PI) phi -= 2*M_PI;
 	  while (phi-phiVec[ist] < -M_PI) phi += 2*M_PI;

	  phiVec[ist] = (phiVec[ist]*nhits[ist]+phi)/(nhits[ist]+1);
	}
	nhits [ist] += 1;
      }
    }

    for (int i=0; i<nstations; i++) {
      if (nhits[i] > 0) {
	zVec  [i] = zVec  [i]/nhits[i];
      }
    }


    if (_debug >5){
      printf("[HelixFitHack::findDfDz]  idStation  nhits       z       phi    \n");
      for (int i=0; i<nstations; i++) {
	if (nhits[i] > 0) {
	  printf("[HelixFitHack::findDfDz]  %8i %3i %10.3f %9.6f   \n",
		 i, nhits[i], zVec[i],  phiVec[i]);
	}
      }
    }

 
    int i0(-1), first_point(1);

    //add the cluster phi
    double zCl   = fTimePeak->ClusterZ();
    pos          = Hep3Vector(fTimePeak->ClusterX(), 
			      fTimePeak->ClusterY(), 
			      fTimePeak->ClusterZ());
    double phiCl = CLHEP::Hep3Vector(pos - center).phi();
    phiCl = TVector2::Phi_0_2pi(phiCl);//if (phiCl<0.0) phiCl +=twopi;
    

    for (int i=0; i<nstations; i++) {
      if (nhits[i] == 0) goto NEXT_POINT;
				        // find station corresponding to the first point
      if (first_point) {
	i0 = i;
	first_point = 0;
      }

      phi_ref = phiVec[i];
      z_ref   = zVec  [i];

      for(int j=i+1; j<nstations; ++j){
	if (nhits[j] == 0) continue;
	phi = phiVec[j];
	z   = zVec  [j];
	dz  = z - z_ref;

	dzOverHelPitch = dz/tollMax - int(dz/tollMax);

	if ( (phi_ref > -9999 ) && 
	     (dzOverHelPitch < dzOverHelPitchCut ) &&
	     //	     (  dz < tollMax  ) &&
	     (  dz > tollMin  )){
	  dphi    = phi - phi_ref;
 	  while (dphi >  M_PI) dphi -= 2*M_PI;
 	  while (dphi < -M_PI) dphi += 2*M_PI;
	                      //add 2 pi for taking into account the fact we are in the second loop
	                      //FIX ME: what to do if we are in the third loop?
	  if (  dz > tollMax  ) dphi += 2*M_PI*int(dz/tollMax);

	  double dphidz = dphi/dz;
	  while (dphidz < 0.) {
	    dphi   += dphi+2.*M_PI;
	    dphidz  = dphi/dz;
	  }
	  _hDfDzRes->Fill(dphidz);
	  
	  double tmpphi0 = phi_ref - dphidz*z_ref;
	  tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);
	  
	  if (_debug >5){
	    printf("[HelixFitHack::findDfDz] z_ref = %9.3f z = %9.3f phi_ref = %9.5f phi = %9.5f dz = %10.3f dz/HelPitch = %10.3f df/dz = %9.5f phi0 = %9.6f\n",
		   z_ref, z, phi_ref, dphi-phi_ref, dz, dzOverHelPitch, dphi/dz, tmpphi0);
	  }
	  
	  //in case dfdz is out of limits
	  // set tmpphi0 as negative
	  if ( (dphidz < _minDfDz) || (dphidz >  _maxDfDz)){
	    tmpphi0 = -1;
	  }
	  _hPhi0Res->Fill(tmpphi0);
	  
	  
	}
      }

      //now use the calorimeter cluster phi
      dz             = zCl - z_ref;
      dzOverHelPitch = dz/tollMax - int(dz/tollMax);

      if ( (phi_ref > -9999 ) && 
//	     (  dz < tollMax  ) &&
	   (  dz > tollMin  )){
	  dphi    = phiCl - phi_ref;
// 	  while (dphi >  M_PI) dphi -= 2*M_PI;
// 	  while (dphi < -M_PI) dphi += 2*M_PI;
	  dphi  = TVector2::Phi_0_2pi(dphi);

	                      //add 2 pi for taking into account the fact we are in the second loop
	                      //FIX ME: what to do if we are in the third loop?
	  if (  dz > tollMax  ) dphi += 2*M_PI*int(dz/tollMax);

	  double dphidz = dphi/dz;
	  while (dphidz < 0.) {
	    dphi   += 2.*M_PI;
	    dphidz = dphi/dz;
	  }
	  
	  double tmpphi0 = phi_ref - dphidz*z_ref;
	  tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

	  if (_debug >5){
	    printf("[HelixFitHack::findDfDz] z_ref = %9.3f z = %9.3f phi_ref = %9.5f phi = %9.5f dz = %10.3f dz/HelPitch = %10.3f df/dz = %9.5f phi0 = %9.6f\n",
		   z_ref, zCl, phi_ref, dphi-phi_ref, dz, dzOverHelPitch, dphidz, tmpphi0);
	  }
	  
	  if (dzOverHelPitch < dzOverHelPitchCut ) {
	    _hDfDzRes->Fill(dphidz);
	    if ( (dphidz < _minDfDz) || (dphidz >  _maxDfDz)){
	      tmpphi0 = -1;
	    }
	    _hPhi0Res->Fill(tmpphi0);
	  }

	
      }
      
    NEXT_POINT:;
    }
    // 2015 - 04- 02 G. Pezzu changed the way the maximum is searched
    // since sometimes 2pi ambig creates two peaks in the histogram
    // we want to use the second, because it is the correct one

    //   int      maxBin  = _hDfDzRes->GetMaximumBin();
    double  maxContent = _hDfDzRes->GetMaximum() - 0.001;
    int      maxBin    = _hDfDzRes->FindLastBinAbove(maxContent);//GetMaximumBin();
    _hdfdz             = _hDfDzRes->GetBinCenter(maxBin);//_hDfDzRes->GetMean();
    double dfdzmean    = _hDfDzRes->GetMean();
    int    nentries    = _hDfDzRes->GetEntries();
    int    overflows   = _hDfDzRes->GetBinContent(0)  + _hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1);

  
//                                             //calculate the mean phi0
//     for (int i=0; i<nstations; i++) {
//       if (nhits[i] == 0) continue;

//       phi_ref = phiVec[i];
//       dz      =  zVec  [i];
//       double phi = phi_ref - _hdfdz*dz;
//       phi        = TVector2::Phi_0_2pi(phi);

//       _hPhi0Res->Fill(phi);
//     }
//                                             //add the cluster contribute
//     phi = phiCl - _hdfdz*zCl;
//     phi = TVector2::Phi_0_2pi(phi);

//     _hPhi0Res->Fill(phi);

    // maxBin              = _hPhi0Res->GetMaximumBin();
    maxContent          = _hPhi0Res->GetMaximum() - 0.001;
    maxBin              = _hPhi0Res->FindLastBinAbove(maxContent);//GetMaximumBin();
  
    double mpvphi0      = _hPhi0Res->GetBinCenter(maxBin);//_hPhi0Res->GetMean();
    double menaphi0     = _hPhi0Res->GetMean();
    int    nentriesphi  = _hPhi0Res->GetEntries();
    //    int    overflowsphi = _hPhi0Res->GetBinContent(0)  + _hPhi0Res->GetBinContent(_hPhi0Res->GetNbinsX()+1);
    _hphi0 = mpvphi0;
    if (_debug > 5){
      
      printf("[HelixFitHack::findDfDz] nentries = %i mpvDfDz = %9.6f meanDphiDz = %9.6f under: %5.0f over: %5.0f \n",
	     nentries, _hdfdz, dfdzmean, 
	     _hDfDzRes->GetBinContent(0),_hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1)
	     );
      for (int i=0; i<_hDfDzRes->GetNbinsX(); i++) {
	printf(" %5.0f",_hDfDzRes->GetBinContent(i+1));
      }
      printf("\n");

      printf("[HelixFitHack::findPhi0] nentries = %i mpvPhi0 = %9.6f meanPhi0 = %9.6f under: %5.0f over: %5.0f \n",
	     nentriesphi, mpvphi0,  menaphi0,
	     _hPhi0Res->GetBinContent(0),_hPhi0Res->GetBinContent(_hPhi0Res->GetNbinsX()+1)
	     );
      for (int i=0; i<_hPhi0Res->GetNbinsX(); i++) {
	printf(" %5.0f",_hPhi0Res->GetBinContent(i+1));
      }
      printf("\n");
    }
    delete     _hDfDzRes;
    _hDfDzRes = 0;
    delete     _hPhi0Res;
    _hPhi0Res = 0;
//-----------------------------------------------------------------------------
// Part 2: try to perform a more accurate estimate - straight line fit
//-----------------------------------------------------------------------------
    double z0, phi0, dphidz, pred;

    z0   = 0.;//zVec  [i0];
    phi0 = _hphi0;//phiVec[i0];

    dphidz = _hdfdz;

    //    _hphi0 = phi0 + _hdfdz*z0;
    _sdfdz = -1;
   
    if (_debug >5) {
      double tmpphi0=phi0+dphidz*z0;
      printf("[HelixFitHack::findDfDz] ------------ Part 2: phi0 = %9.6f dfdz = %9.6f\n", tmpphi0, dphidz);
    }

    //--------------------------------------------------------------------------------
    // 2015-03-25 G. Pezzu changed theway the 2Phi ambiguity is resolved
    //--------------------------------------------------------------------------------
    
    double xdphi;

    ::LsqSums4 srphi;
    srphi.clear();
    
   
    double weight = 1./(_sigmaPhi*_sigmaPhi);
  
    double zLast(z0), zdist;
    
    if (_debug >5) {
      xdphi = fabs(dphi)/(2./30.);//_sigmaPhi;
      double tmpphi0=phi0+dphidz*z0;
      printf("[HelixFitHack::findDfDz]  i       z        dz      phiVec     phi    srphi.dfdz    dphi    xdphi     dfdz    \n");
      printf("[HelixFitHack::findDfDz] phi0 = %9.6f dfdz = %9.6f\n", tmpphi0, dphidz);
    }
       
    dz   = zCl - z0;
    dphi = dz*dphidz + phi0 - phiCl;
    while (dphi > M_PI){
      phiCl += 2*M_PI;
      dphi  -= 2*M_PI; // dz*dphidz + phi0 - phiCl;
    }
    while (dphi < -M_PI){
      phiCl -= 2*M_PI;
      dphi  += 2*M_PI; // = dz*dphidz + phi0 - phiCl;
    }
      
    double errCl = 2./30.;
    double weightCl = 1./(errCl*errCl);

    xdphi = fabs(dphi)/errCl;

    //2015-04-21 Gianipez added the condition (xdphi < 2.*_maxXDPhi) for adding or not 
    //the calorimeter point to the fitter. In case the particle scattered in the ened of the tracker
    //the calorimeter point is dangerous. 
    if (xdphi < 2.*_maxXDPhi){ 
      srphi.addPoint(zCl, phiCl, weightCl);
    }

    if (_debug >5) {
      printf("[HelixFitHack::findDfDz] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
	     0, zCl, dz, phiCl, phiCl, srphi.dfdz(), dphi, xdphi, dphidz);
    }
    
    //2015-07-06 Gianipez added the following line for avoiding infinite loops
    if ( i0 < 0 ) goto NEXT_STEP;

    for (int i=i0; i<nstations; i++) {
      if (nhits[i] > 0) {
	z    = zVec[i];
	dz   = z-z0;
	pred = phi0 + dz*dphidz;
	phi  = phiVec[i];
	dphi = pred - phi;

	while (dphi > M_PI){
	  phi += 2*M_PI;
	  dphi = pred - phi;
	}
	while (dphi < -M_PI){
	  phi -= 2*M_PI;
	  dphi = pred - phi;
	}

	xdphi = fabs(dphi)/_sigmaPhi;

	if (xdphi < 2.*_maxXDPhi){
	  srphi.addPoint(z, phi, weight);

	  zdist = z - zLast;
	  
	  if ( (srphi.qn() >= 3.) && (zdist > 500.)){
	    z0     = 0.;
	    phi0   = srphi.phi0();
	    dphidz = srphi.dfdz();
	    //	  zLast  = z; 
	  }
	}
	

	if (_debug >5) {

	  double tmpDfDz = srphi.dfdz();//, Helix._srphi.chi2DofLine());
	  printf("[HelixFitHack::findDfDz] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
		 i, z, dz, phiVec[i], phi, tmpDfDz, dphi, xdphi, dphidz);
	}
	
      }
    }

  NEXT_STEP:;

    if (srphi.qn() >= 3.){
     
      _hdfdz = srphi.dfdz();//sigxy/sigxx;
      _hphi0 = srphi.phi0();//ymean - xmean*sigxy/sigxx;
      _sdfdz = srphi.chi2DofLine();
    }else {
      _hphi0 = phi0 + _hdfdz*z0;
      _sdfdz = -1;
    }
    
    if (_debug >5) {
      printf("[HelixFitHack::findDfDz] END: _hdfdz = %9.5f _hphi0 = %9.6f chi2 = %9.3f ", _hdfdz, 
	     _hphi0, -1.);
      printf(" FIT: srphi.dfdz() = %9.5f srphi.phi0() = %9.6f chi2 = %9.3f qn = %6.0f\n", srphi.dfdz(), 
	     srphi.phi0(), srphi.chi2DofLine(), srphi.qn());
    }

    if ( (nentries - overflows) == 0) _hdfdz = _mpDfDz;

    return 1;
    // if ( (nentries - overflows) > 8 ) {
//       return 1;
//     } else {
//       _hdfdz = _mpDfDz;
//       return 0;
//     }
  }
    



//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  bool HelixFitHack::doLinearFitPhiZ(HelixFitHackResult& Helix    , 
				     int                 SeedIndex, 
				     int                *IndexVec,
				     int                 UseInteligentWeight) {
    // sort points by z
    //    std::sort(_xyzp.begin(),_xyzp.end(),zcomp());

    bool   success(false);
    int    nPointsRemoved(0);
    int    N = _xyzp.size();
    int    idVec[N];
    for (int i=0; i<N; ++i){
      idVec[i] = IndexVec[i];
    }
    
    double phi_corrected[N];
    double phi(0.0), z(0.0), weight(0.0);
    double dphi, err, xdphi;
    
    CLHEP::Hep3Vector center = Helix._center;
    CLHEP::Hep3Vector pos;

    // calculate the mean pitch assuming conversion electron 
    // hypothesis (so assuming a pitch angle of ~ 0.65)

    ///    double pitch = twopi*std::pow(Helix._radius,2.)*Helix._dfdz; //*tan(0.65);

  
//-----------------------------------------------------------------------------
// gianipez: procedure for aligning the phi vector
//-----------------------------------------------------------------------------
    ::LsqSums4 srphi;
    int        iworst, count(0), indexWorst;
    double     chi2,chi2min, dfdz, deltaPhi, dphi_max(0), phi_ref, weightWorst(0);
  

    CLHEP::Hep3Vector helCenter( Helix._sxy.x0(), Helix._sxy.y0(), 0), strawDir;
    double            radius(Helix._sxy.radius());
    const char        banner[] = "doLinearFitPhiZ";

    //    if (_debug > 5) printf("[HelixFitHack::doLinearFitPhiZ] _xyzp.size() = %u\n", N);
//--------------------------------------------------------------------------------
// set EMC cluster info
//-----------------------------------------------------------------------------
      //initilize the dfdz for the search
    dfdz = Helix._dfdz;
    double phi0 = Helix._fz0;
    
    double zCl   = fTimePeak->ClusterZ();
    pos          = Hep3Vector(fTimePeak->ClusterX(), 
			      fTimePeak->ClusterY(), 
			      fTimePeak->ClusterZ());
    double phiCl = CLHEP::Hep3Vector(pos - center).phi();
    phiCl = TVector2::Phi_0_2pi(phiCl);//if (phiCl<0.0) phiCl +=twopi;

    deltaPhi = zCl*dfdz + phi0 - phiCl;
    while (deltaPhi > M_PI){
      phiCl += 2*M_PI;
      deltaPhi = zCl*dfdz + phi0 - phiCl;
    }
    while (deltaPhi < -M_PI){
      phiCl -= 2*M_PI;
      deltaPhi = zCl*dfdz + phi0 - phiCl;
    }
//-----------------------------------------------------------------------------
// add the cluster phi to the LSq sum
// weight_cl of 10 corresponds to an uncertainty of 0.1 in phi(cluster), 
// which means sigma(R-phi) of about 2-3 cm, about right
// hit weight is determined by the assumed phi error of _sigmaPhi=0.1
//-----------------------------------------------------------------------------
    if (UseInteligentWeight == 0){
      weight = 1./(_sigmaPhi*_sigmaPhi);
    }
    
    xdphi  = fabs(deltaPhi)/_sigmaPhi;
    
    Helix._srphi.clear();
    if ( xdphi < 2.*_maxXDPhi ) {
      double weight_cl = 10.0;
      Helix._srphi.addPoint(zCl,phiCl,weight_cl);
    }
 
    count = 0;
    double zlast, dz;
    
    if (_debug > 5) {
      printf("[HelixFitHack::doLinearFitPhiZ]    flag    A   sh-id   hit-id    z     ");
      printf("   phi         dphi        xdphi     zlast      dz       dphidz      szphidfdz      chi2  \n");
      printf("[HelixFitHack::doLinearFitPhiZ] phi0 = %10.6f dfdz = %10.6f chi2N = %10.3f\n",
	     Helix._fz0,  Helix._dfdz, 0.);
      printf("[HelixFitHack::doLinearFitPhiZ] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
	     0, 1, 0, 0,  zCl, phiCl, deltaPhi, xdphi, 0., 0., dfdz, 0., 0.);
     
    }

    zlast = 0;
    for(int i=SeedIndex; i<N; ++i){
      pos                  = _xyzp[i]._pos;
      z                    = pos.z();
      strawDir             = _xyzp[i]._sdir;

      phi      = CLHEP::Hep3Vector(pos - center).phi();
      phi      = TVector2::Phi_0_2pi(phi);
      dz       = z - zlast;

                                    //calculate the predicted value of phi
      phi_ref = z*dfdz + phi0;
                                    //calculate the signed residual
      dphi    = phi_ref - phi;
                                    //resolve the 2PI ambiguity
      while (dphi > M_PI){
	phi += 2*M_PI;
	dphi = phi_ref - phi;
      }
      while (dphi < -M_PI){
	phi -= 2*M_PI;
	dphi = phi_ref - phi;
      }
      
      // store the corrected value of phi
      phi_corrected[i] = phi;
      fPhiCorrected[i] = phi;

      dphi             = fabs(dphi);
      err              = _sigmaPhi;

      if (UseInteligentWeight == 1){
	weight = calculatePhiWeight(pos, strawDir, helCenter, radius, 0, banner);
	err    = 1./sqrt(weight);
      }

      xdphi = dphi/err;

      if (_debug > 5) printf("[HelixFitHack::doLinearFitPhiZ] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
			     *((int*) &_xyzp[i]._flag), idVec[i] < 0 ? 0 : idVec[i],
			     _xyzp[i]._strawhit->strawIndex().asInt()/*int(_xyzp[i]._ind)*/, i, 
			     z, phi_corrected[i], dphi,xdphi,zlast,dz,
			     dfdz, Helix._srphi.dfdz(), Helix._srphi.chi2DofLine());
      
      if (_xyzp[i].isOutlier())                        continue;
      if (  idVec[i] < 1     )                         continue;

      if ( xdphi > _maxXDPhi ) {
	idVec[i] = 0;
	++nPointsRemoved;
	                                               continue;
      }
      Helix._srphi.addPoint(z, phi_corrected[i], weight);
      ++count;

      if (count == 1) {
	zlast = z;
	dz    = 0.;
      }
      
      if ( (count>=5) && 
	   (xdphi < 2.)){
	if ( (fabs(dfdz - Helix._srphi.dfdz()) < 8.e-4) ){//  || //require that the new value of dfdz is 
	                                                    //close to the starting one. update dfdz only if:
	  if ( (Helix._srphi.dfdz() > 0.) && //{                    // 1. the points browsed are more the half
	       (dz >=500 ) ){
	    phi0 = Helix._srphi.phi0();                     // 2. and require dfdz to be positivie! scattered hits or 
	    dfdz = Helix._srphi.dfdz();                     //    delta hits could have moved dfdz to negative value!
	    zlast   = z;
	  }
	}
      }
    }

    fPhiCorrectedDefined = 1;

    if (_debug > 5) {
      printf("[HelixFitHack::doLinearFitPhiZ] phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f points removed = %4i\n", 
	     Helix._srphi.phi0(),Helix._srphi.dfdz(), Helix._srphi.chi2DofLine(), nPointsRemoved);
    }
//-----------------------------------------------------------------------------
// perform a cleanup in RZ
//-----------------------------------------------------------------------------
    if ( Helix._srphi.chi2DofLine() > _chi2zphiMax) {
    NEXT_ITERATION:;
      iworst      = -1;
      indexWorst  = -1;
      chi2min     = 1e10;
      weightWorst = -1;
      for(int ixyzp=SeedIndex; ixyzp < N; ++ixyzp){
	if (_xyzp[ixyzp].isOutlier())                        continue;
	if (  idVec[ixyzp] < 1  )                            continue;
	srphi.init(Helix._srphi);
	pos      = _xyzp[ixyzp]._pos;
	strawDir = _xyzp[ixyzp]._sdir;

	z        = pos.z();
	phi      = phi_corrected[ixyzp];

	if (UseInteligentWeight == 1){
	  weight = calculatePhiWeight(pos, strawDir, helCenter, radius, 0, banner);
	}

	srphi.removePoint(z, phi, weight);
	chi2 = srphi.chi2DofLine();
	//printf("[HelixFitHack::doLinearFitPhiZ] chi2 = %5.3e chi2min = %5.3e\n", chi2, chi2min);
	if (chi2 < chi2min) {
	  iworst      = ixyzp;
	  indexWorst  = _xyzp[ixyzp]._ind;
	  chi2min     = chi2;
	  weightWorst = weight;
	}
      }
  
      if ((iworst >= 0) && (Helix._srphi.qn() > 10.)) {
	idVec[iworst] = 0;
	pos = _xyzp[iworst]._pos;
	z   = pos.z();
	phi = phi_corrected[iworst];
    
	Helix._srphi.removePoint(z, phi, weightWorst);
	chi2min = Helix._srphi.chi2DofLine();
	if (_debug > 5) {
	  printf("[HelixFitHack::doLinearFitPhiZ_removed] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);
	}
      }

    CHECK_RESIDUALS:;
      dphi_max    = _maxXDPhi;
      iworst      = -1;
      weightWorst = -1;

      for(int i=SeedIndex; i < N; ++i){
	if (_xyzp[i].isOutlier())     continue;
	if (  idVec[i] < 1  )         continue;

	pos      = _xyzp[i]._pos;
	z        = pos.z();
	strawDir = _xyzp[i]._sdir;

	phi      = z* Helix._srphi.dfdz() + Helix._srphi.phi0();
	dphi     = fabs(phi_corrected[i] - phi);

	err      = _sigmaPhi;

	if (UseInteligentWeight == 1){
	  weight = calculatePhiWeight(pos, strawDir, helCenter, radius, 0, banner);
	  err    = 1./sqrt(weight);
	}
	
	xdphi = dphi/err;

	if ( xdphi > dphi_max) {
	  iworst      = i;
	  indexWorst  = _xyzp[i]._ind;
	  dphi_max    = xdphi;
	  weightWorst = weight;
	}
      }

      //remove the point
      if(iworst>=0 && Helix._srphi.qn() > 10.){
	idVec[iworst] = 0;//.setOutlier();
	pos = _xyzp[iworst]._pos;
	z   = pos.z();
	phi = phi_corrected[iworst];//CLHEP::Hep3Vector(pos - center).phi();
	
	Helix._srphi.removePoint(z, phi, weightWorst);
	chi2min = Helix._srphi.chi2DofLine();
	if (_debug > 5) {
	  printf("[HelixFitHack::doLinearFitPhiZ_removed] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);
	}
	goto CHECK_RESIDUALS;
      }
      
      
      if(Helix._srphi.qn()<=10) chi2min = Helix._srphi.chi2DofLine();
  
      if ( (chi2min >= _chi2zphiMax) ||
	   (iworst>=0 )) {
	
//--------------------------------------------------------------------------------
// 2016-04-27 gianipez: why should I not check the chi2 f I have 10 hits?
//--------------------------------------------------------------------------------
	if (Helix._srphi.qn() > 10.) {
	  goto NEXT_ITERATION;
	}
      }
    }

    //2015-04-21 Gianipez changed the threshold from 3 to _minnhit. there is no reason 
    // which should allow to keep a result which selects a number of points lower than the threshold!
    //    if ( (Helix._srphi.qn() >= 3) && (Helix._srphi.chi2DofLine() < _chi2zphiMax) ){
    if ( (Helix._srphi.qn() >= _minnhit) && (Helix._srphi.chi2DofLine() < _chi2zphiMax) ){
      success = true;
    }
    //----------------------------------------------------------------------//
    if ( Helix._srphi.dfdz() < 0.) {
      success = false;
    }else if (success){                               //update helix results
      Helix._fz0  = Helix._srphi.phi0();
      Helix._dfdz = Helix._srphi.dfdz();
    }

    if (SeedIndex ==0){
      THackData* hack;
      hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
      hack->fData[6]  = Helix._srphi.phi0();
      hack->fData[7]  = Helix._srphi.dfdz()*hack->fData[9];
      hack->fData[8]  = Helix._srphi.dfdz();//_dfdz;//_mpDfDz; 
      hack->fData[13] = Helix._srphi.chi2DofLine();

      if (success){
	int h=0;
	for (int i=SeedIndex; i<N; ++i){
	  if (_xyzp[i].isOutlier())     continue;
	  if (  idVec[i] < 1  )        continue;
	  pos      = _xyzp[i]._pos;
	  z        = pos.z();
	  phi      = z* Helix._dfdz + Helix._fz0;
	  deltaPhi = phi_corrected[i] - phi;
	  
	  hack->fResid[h] = deltaPhi;
	  ++h;
	}
      }
    }
    
 
    
    if (_debug > 5) {
      printf("[HelixFitHack::doLinearFitPhiZ] Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f\n", 
	     Helix._srphi.phi0(),Helix._srphi.dfdz(), Helix._srphi.chi2DofLine() );
      printf("[HelixFitHack::doLinearFitPhiZ] srphi: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f\n", 
	     srphi.phi0(), srphi.dfdz(), srphi.chi2DofLine() );
    }
    
    _chi2nFindZ = Helix._srphi.chi2DofLine();
    if(_chi2nFindZ < 0.0)  _eventToLook = 0; // FIXME!!

    if (_debug > 5) {
      printf("[HelixFitHack::doLinearFitPhiZ] retval = %d\n",success ? 1:0);
      printf("[HelixFitHack::doLinearFitPhiZ]    flag    A   sh-id       z         phi      phi-dfdz*z-phi0\n");

      for(int i=N-1; i>=SeedIndex; --i){
      
	pos                  = _xyzp[i]._pos;
	z                    = pos.z();
	
	phi = z* Helix._dfdz + Helix._fz0;
	
	deltaPhi = phi_corrected[i] - phi;
	
	printf("[HelixFitHack::doLinearFitPhiZ] %08x %2i %6i %12.5f %12.5f %12.5f\n",
	       *((int*) &_xyzp[i]._flag), idVec[i] < 0 ? 0 : idVec[i], 
	       _xyzp[i]._strawhit->strawIndex().asInt(),  z, phi_corrected[i], deltaPhi);
      }

    }
    
    if (success){
      for (int i=0; i<N; ++i){
	IndexVec[i] = idVec[i];
      }
    }
    return success;
  }

//-----------------------------------------------------------------------------
// 12-09-2013 gianipez modified this procedure to avoid the doubling of the 
// same stereohitposition 
//-------------------------------------------------------------------------
  void HelixFitHack::fillXYZP(HelixDefHack const& mytrk){
      
    //clear xyzp vector
    _xyzp.clear();
    
    const Tracker& tracker = getTrackerOrThrow();
    
    const std::vector<hitIndex> shIndices = mytrk.strawHitIndices();
    
    int size = shIndices.size();

    //--------------------------------------------------------------------------------
    if (mytrk.strawHitPositionCollection() != 0) {
      int loc;
      StrawHitFlag flag;
      for (int i=0; i<size; ++i) { 
	loc                = shIndices[i]._index;
	flag               = mytrk.strawHitFlagCollection()->at(loc);
	StrawHit const& sh = mytrk.strawHitCollection()->at(loc);
	Straw const& straw = tracker.getStraw(sh.strawIndex());
	StrawHitPosition const& shp = mytrk.strawHitPositionCollection()->at(loc);
//-----------------------------------------------------------------------------
// don't reuse straw hits 
//-----------------------------------------------------------------------------
	if (! flag.hasAnyProperty(StrawHitFlag::calosel)) {
	  XYZPHack pos(loc,sh,shp,straw,flag);
	  if (_debug > 0) {
	    if(i == 0){
	      printf("----->HelixFitHack::fillXYXP\n");
	      printf("   strawIndex       straw ID      stereo index        pos          _debug: %i5\n",_debug);
	    }					
	    printf("%10d %10d %10d %10.3f %10.3f %10.3f \n", 
		   (int) loc, sh.strawIndex().asInt(), shp.stereoHitIndex(),
		   shp.pos().x(), shp.pos().y(), shp.pos().z());
	  }
	  _xyzp.push_back(pos);
	}
      }
    }
//----------------------------------------------------------------------
// 2014-11-06 gianipez added the following line for ordering the xyzp
// strawhits along their z coordinate
//----------------------------------------------------------------------
    std::sort(_xyzp.begin(), _xyzp.end(), [ ]( const XYZPHack& lhs,
					       const XYZPHack& rhs )
	      {
		return lhs._pos.z() < rhs._pos.z();
	      } );
  }


  

//----------------------------------------------------------------------------
//2015-01-17 G. Pezzullo: the following procedure looks the hit with 
// z-coordinate smaller then the seeding one and calculates distance from
// prediction in order to check if they are good or outliers
//----------------------------------------------------------------------------
  void HelixFitHack::rescueHitsBeforeSeed(HelixFitHackResult&  mytrk){

    double      weight, radius, phi0, dfdz, x0, y0;
    dfdz        = mytrk._dfdz;
    phi0        = mytrk._fz0 + dfdz*(_xyzp[fSeedIndex]._pos.z());
    x0          = mytrk._center.x();
    y0          = mytrk._center.y();
    radius      = mytrk._radius;

    double      dx,dy,phi,dx2, dy2, max_dist;
    Hep3Vector  shPos, hePos;
    
    double deltaZ(0.); // , deltaX(0.), deltaY(0.);
    double distXY(0.0);
    double dist(0.0), dist2(0.0); //help parameter for storing strawhit position residual
    int    i_last(fSeedIndex), rescuedPoints(0);
    
    TString banner="HelixFitHack::rescueHitsBeforeSeed";

    if (_debug > 0) {
      printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.5f dfdz = %5.6f chi2 = %5.3f \n", banner.Data(),
	     x0, y0, radius, phi0, dfdz , mytrk._sxy.chi2DofCircle());
      printf("[%s] SeedIndex = %i N-points = %5.3f\n",  banner.Data(), fSeedIndex, mytrk._sxy.qn()-1);
    }
    
    for (int i=fSeedIndex-1; i>=0; --i){
      if (_xyzp[i].isOutlier())                              goto NEXT_POINT;
      weight = 1.;
      shPos  = _xyzp[i]._pos;
      deltaZ = shPos.z() - _xyzp[i_last]._pos.z();
      phi    = phi0 + (deltaZ)*dfdz;                     // tanLambda/radius;

      hePos  = Hep3Vector(x0 + radius*std::cos(phi),
			  y0 + radius*std::sin(phi),
			  shPos.z());

      dx  = hePos.x() - shPos.x();
      dx2 = dx*dx;
      dy  = hePos.y() - shPos.y();
      dy2 = dy*dy;
      
      dist2 = dx2 + dy2;
      dist  = std::sqrt(dist2);

      if( _debug>0){
	
	printf("[%s]   measured     %10.3f  %10.3f  %10.3f  %8i \n", banner.Data(), 
	       shPos.x(), shPos.y(), shPos.z(), i);
	printf("[%s]  predicted     %10.3f  %10.3f  %10.3f  %8i \n", banner.Data(), 
	       hePos.x(), hePos.y(), hePos.z(), i);
	printf("[%s] X0 = %5.3f Y0 = %5.3f r = %5.3f dfdz = %5.5f  dist-from-prediction = %5.3f  dist-from-seedXY = %5.3f dz-from-seed = %5.3f\n", banner.Data(), 
	       x0, y0, radius, dfdz, dist, distXY, deltaZ);
      }
 
      max_dist = _distPatRec + _dfdzErr*fabs(deltaZ);
      if( dist <= max_dist ){

	//store index of last good point
	i_last = i;
	phi0   =  phi;//CLHEP::Hep3Vector(shPos - CLHEP::Hep3Vector(x0, y0, 0.)).phi();

	//add point to the helixfithack result objet
	mytrk._sxy.addPoint(shPos.x(),shPos.y(), weight);
      
	//store the index of the good point found
	_indicesTrkCandidate[i] = 1;

	//store distance along z-axis from the last point found relying on thr helix
	_dzTrkCandidate     [i] = deltaZ;

	//store distance from predition
	_distTrkCandidate   [i] = dist;

	//update helix parameters
	x0      = mytrk._sxy.x0();
	y0      = mytrk._sxy.y0();
	radius  = mytrk._sxy.radius();

	++rescuedPoints;
	
	if( _debug>0){
	  printf("[%s] rescued %08x %2i %12.5f %12.5f %12.5f \n", banner.Data(),
		 *((int*) &_xyzp[i]._flag), _indicesTrkCandidate[i], 
		 _xyzp[i]._pos.x(), _xyzp[i]._pos.y(), _xyzp[i]._pos.z()
		 );
	}
      }

   NEXT_POINT:;
    }

    //update mytrk info
    
    mytrk._center.set(x0, y0, 0.0);
    mytrk._radius = radius;

    _goodPointsTrkCandidate = mytrk._sxy.qn() - 1;//removing the EMC cluster!

    banner += "-results";
    if (_debug > 5) {
      printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.3f dfdz = %5.6f chi2 = %5.3f \n", banner.Data(),
	     x0, y0, radius,  mytrk._fz0, dfdz , mytrk._sxy.chi2DofCircle());
      printf("[%s] SeedIndex = %i N rescued points = %i\n",  banner.Data(), fSeedIndex, rescuedPoints);
    }

    THackData* hack;
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    hack->fData[16] =   rescuedPoints;  
  }

//--------------------------------------------------------------------------------
  void HelixFitHack::printInfo(HelixFitHackResult& mytrk){
    const char banner [] = "HelixFitHack::printInfo";
    double dr(0), dx, dy;

    if (_debug > 0) {
      printf("[%s] N(points): %3.0f x0: %12.5f y0: %12.5f r: %12.5f chi2c: %12.5f phi0: %5.5f dfdz: %5.6f chi2l: %5.3f\n", 
	     banner, 
	     mytrk._sxy.qn() - 1,
	     mytrk._sxy.x0(),mytrk._sxy.y0(),mytrk._sxy.radius(),mytrk._sxy.chi2DofCircle(), 
	     mytrk._fz0, mytrk._dfdz , mytrk._srphi.chi2DofLine());

      int np = _xyzp.size();
      for (int i=0; i<np; i++) {
	dx = _xyzp[i]._pos.x() - mytrk._sxy.x0();
	dy = _xyzp[i]._pos.y() - mytrk._sxy.y0();
	dr = sqrt(dx*dx+dy*dy) - mytrk._sxy.radius();
	printf("[%s] %08x %6i %3i %6i %12.5f %12.5f %12.5f %10.3f\n",banner,
	       *((int*) &_xyzp[i]._flag),  _indicesTrkCandidate[i], i, _xyzp[i]._strawhit->strawIndex().asInt(),
	       _xyzp[i]._pos.x(), _xyzp[i]._pos.y(), _xyzp[i]._pos.z(), dr
	       );
      }
    }
  }

//-----------------------------------------------------------------------------
// this routine simply checks '_indicesTrkCandidate' array and for negative 
// indices sets the 'outlier' flag to the corresponding 'xyzp'
// no actual check of residuals is performed
//-----------------------------------------------------------------------------
  void HelixFitHack::filterUsingPatternRecognition(HelixFitHackResult& mytrk) {
    int np = _xyzp.size();
    if (fSeedIndex < 0) return;
    Hep3Vector pSeed = _xyzp[fSeedIndex]._pos;
    //if no second strawhit has been found, use the target center in the transverse plane
    Hep3Vector pCand(0.0, 0.0, 0.0);

    if (fCandidateIndex >= 0){
      pCand = _xyzp[fCandidateIndex]._pos;
    }
    Hep3Vector shPos;
    double dist, dz;
    int    nActive(0);
    //    for (int i=fSeedIndex; i<np; ++i){
    for (int i=0; i<np; ++i){
      if (_debug>0){
	//	if( i==fSeedIndex+1) {
	if( i==0) {
	  dist = 0;
	  dz   = pCand.z() - pSeed.z();
	  printf("[HelixFitHack::filterUsingPatternRecognition]  filterUsingPatternRecognition() will set asOutlier the following hits using helix parameters\n");
	  printf("[HelixFitHack::filterUsingPatternRecognition] X0 = %5.3f Y0 = %5.3f r = %5.3f chi2N = %5.5f phi0 = %5.5f dfdz = %5.5f chi2N = %5.5f straw-hits = %i\n", 
		 mytrk._sxy.x0(), mytrk._sxy.y0(), _radius, mytrk._sxy.chi2DofCircle(), 
		 mytrk._srphi.phi0(), mytrk._srphi.dfdz(), mytrk._srphi.chi2DofLine(), 
		 _goodPointsTrkCandidate );// +1 for counting also the seeding hit
	  printf("[HelixFitHack::filterUsingPatternRecognition]   point  type     X           Y           Z      xyzp-index    hit index    dist       Dz    \n");
	  printf("[HelixFitHack::filterUsingPatternRecognition] ----------------------------------------------------------\n");
	  printf("[HelixFitHack::filterUsingPatternRecognition]    seeding     %8.3f  %8.3f  %10.3f  %10i \n",
		 pSeed.x(), pSeed.y(), pSeed.z(), fSeedIndex);
	  printf("[HelixFitHack::filterUsingPatternRecognition]   candidate    %8.3f  %8.3f  %10.3f  %10i \n",
		 pCand.x(), pCand.y(), pCand.z(), fCandidateIndex);
	}
      }
      shPos = _xyzp[i]._pos;
      dist  = _distTrkCandidate[i];
      dz    = _dzTrkCandidate  [i];
      if ( _indicesTrkCandidate[i] <= 0 ){
	_xyzp[i].setOutlier();
	
	if (_debug>0){
	  printf("[HelixFitHack::filterUsingPatternRecognition]   outlier      %10.3f  %10.3f  %10.3f  %10i  %7i  %10.3f  %10.3f \n", 
		 shPos.x(), shPos.y(), shPos.z(), i, _xyzp[i]._strawhit->strawIndex().asInt(), dist, dz);
	}
      }else {
	++nActive;
	if (_debug>0){
	  printf("[HelixFitHack::filterUsingPatternRecognition]    active      %10.3f  %10.3f  %10.3f  %10i  %7i  %10.3f  %10.3f \n", 
		 shPos.x(), shPos.y(), shPos.z(), i, _xyzp[i]._strawhit->strawIndex().asInt(), dist, dz);
	}
      }
    } 
    if (_debug>0){
      printf("[HelixFitHack::filterUsingPatternRecognition]    ended: N active point = %i over N-hits = %i\n", nActive, np); 
    }
    _goodPointsTrkCandidate = nActive;
    
    THackData* hack;
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    hack->fData[11] = nActive;

  }
  
//-----------------------------------------------------------------------------
  void HelixFitHack::filterDist() {
    using namespace boost::accumulators;
    static const double pi(M_PI);
    static const double twopi(2*pi);

    double mphi(-9999.);

    if (fTimePeak->ClusterT0() > 0.) {
      mphi = atan2(fTimePeak->ClusterY(),fTimePeak->ClusterX());
    }
    
    if (_debug > 0) {
      printf("[HelixFitHack::filterDist] x0 = %12.5f y0 = %12.5f r = %12.5f chi2 = %12.5g mphi = %12.5g\n",
	     -1., -1., -1., -1., mphi);
    }

    int np = _xyzp.size();
    for (int ip=0; ip < np; ip++) {
      if(_xyzp[ip].use()){
	double dphi = _xyzp[ip]._phi - mphi;

	if (dphi >  pi) dphi -= twopi;
	if (dphi < -pi) dphi += twopi;

	if(fabs(dphi) > pi/2){
	  _xyzp[ip].setOutlier();
	}
      }
      if (_debug > 0) {
	printf("[HelixFit::filterDist]: %08x %3i %12.5f %12.5f %12.5f %12.5f \n",
	       *((int*) &_xyzp[ip]._flag),
	       ip,
	       _xyzp[ip]._pos.x(), _xyzp[ip]._pos.y(),_xyzp[ip]._pos.z(),
	       _xyzp[ip]._phi);
      }
    }
  }


//-----------------------------------------------------------------------------
  double HelixFitHack::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }

//-----------------------------------------------------------------------------
  void HelixFitHack::doPatternRecognition(HelixFitHackResult& Helix) {
    int np = _xyzp.size();
    //int SeedIndex(-1);//, targetIndex(-1);
    int mode(0);                   //flag for indexing event where a second strawhit has been included in the patern reco

    double chi2;
    int countGoodPoints(0);

    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: BEGIN\n");
    
    if (_debug2 == 0){
      _debug2 = _debug; 
      _debug  = 0;
    }

    for (int i=0; i<np; i++) {
      if (_xyzp[i].isOutlier()) goto NEXT_POINT;
      if (_xyzp[i].isCalosel()) goto NEXT_POINT;
//----------------------------------------------------------------------
// 2014-12-26 gianipez: don't start the search from an already used hit
// used in previous search
//-----------------------------------------------------------------------------
      if ( isHitUsed(i) == 1 )  goto NEXT_POINT;

      if (_debug > 10) printf("[HelixFitHack::doPatternRecognition]: calling findTrack i=%3i\n",i);
      if ( (np -i) > _goodPointsTrkCandidate){
	findTrack(i, chi2, countGoodPoints, Helix, mode, false); 
      }
//------------------------------------------------------------------------------
// 2015-01-22 P.Murat: what happens when the very first candidate is good enough ?
//                     where is the comparison of the found candidate with the best previous one ?
//-----------------------------------------------------------------------------
    NEXT_POINT:;
    }
//-----------------------------------------------------------------------------
// 2014-11-09 gianipez: if no track was found requiring the recalculation of dfdz
// look for a track candidate using the default value of dfdz and the target center
//-----------------------------------------------------------------------------
    if ( fUseDefaultDfDz == 0){
      for (int i=0; i<np; i++) {
	if (_xyzp[i].isOutlier()) goto NEXT_P;
	if (_xyzp[i].isCalosel()) goto NEXT_P;
	if (_debug >5 ) printf("[HelixFitHack::doPatternRecognition]: fUseDefaultDfDz=0, calling findTrack i=%3i\n",i);
	if ( (np -i) > _goodPointsTrkCandidate){
	  findTrack(i, chi2, countGoodPoints, Helix, mode, true);
	}
      NEXT_P:;
      } 
    }

    //2015-01-14 G. Pezzullo added the findDfDz procedure
    if (_debug > 5) printf("[HelixFitHack::doPatternRecognition]: ------------ calling findDfDz\n");

    if (fSeedIndex >= 0) {
      findDfDz(Helix,fSeedIndex,_indicesTrkCandidate);
      if (_debug >5) printf("[HelixFitHack::doPatternRecognition]: findDfDz ----> phi0 = %5.5f dfdz = %5.5f \n",
			      _hphi0, _hdfdz);
    } 
    else {
      //if SeedIndex is < 0 it means that no cancidate has been found
      // cases when it occurs are usually the one where the cluster is not in the trajectory or 
      // la very low number of hits is in the time peak
      // maybe we should set a threshold on the time peak size to avoid such?
      int np = _xyzp.size();
      int vIndices[np];
      for (int i=0; i<np; ++i){
	vIndices[i] = 1;
      }
      findDfDz(Helix, 0, vIndices);
      if (_debug > 5) {
	printf("[HelixFitHack::doPatternRecognition]: findDfDz called using SeedIndex = 0 and using all hits (expect outliers!) \n");
	printf("[HelixFitHack::doPatternRecognition]: findDfDz ----> phi0 = %5.5f dfdz = %5.5f \n",
	       _hphi0, _hdfdz);
      }
    }
//-----------------------------------------------------------------------------
// 2016-01-29 P.Murat:
// at this point, with dfdz in hand, initialize 'fPhiCorrected' - it is used 
// in 'rescueHits' and is not initialized upon first entry
// findDfDZ calculates '_hdfdz' and '_hphi0', use those
//-----------------------------------------------------------------------------
    CLHEP::Hep3Vector center = Helix._center;
    CLHEP::Hep3Vector pos_ref, pos;
    double z, phi, phi_ref, dphi;
    
    for(int i=0; i<np; ++i){
      pos = _xyzp[i]._pos;
      z   = pos.z();
      
      phi = CLHEP::Hep3Vector(pos - center).phi();
      phi = TVector2::Phi_0_2pi(phi);
                                    // predicted value of phi
      phi_ref = z*_hdfdz + _hphi0;
                                    // signed residual
      dphi    = phi_ref - phi;
                                    // resolve the 2PI ambiguity
      while (dphi > M_PI) {
	phi += 2*M_PI;
	dphi = phi_ref - phi;
      }
      while (dphi < -M_PI) {
	phi -= 2*M_PI;
	dphi = phi_ref - phi;
      }
                                    // store the corrected value of phi
      fPhiCorrected[i] = phi;
    }
					// don't know
    fPhiCorrectedDefined = 0;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
    int useMPVdfdz = 1;
    for (int i=0; i<np; i++) {
      if (_xyzp[i].isOutlier()) goto NEXT_HIT;
      if (_xyzp[i].isCalosel()) goto NEXT_HIT;
	if (_debug > 5) printf("[HelixFitHack::doPatternRecognition]: useMPVdfdz=1, calling findTrack i=%3i\n",i);
	if ( (np -i) > _goodPointsTrkCandidate){
	  findTrack(i, chi2, countGoodPoints, Helix, mode, false, useMPVdfdz);
	}
    NEXT_HIT:;
    } 

    if (_debug == 0){
      _debug  = _debug2; 
      _debug2 = 0;
    }

    TString banner;
    bool rc;
    int  rc1, refineHelParamRes, rs;
    int  usePhiResid;
    int  useInteligentWeight(1);

    if ( (fSeedIndex<0) || (_goodPointsTrkCandidate < 5) )                     goto  PATTERN_RECOGNITION_END;

    // 2015-01-17 G. Pezzullo added the following procedure to rescue points with z-coordinate 
    // less than the seed hit
    if (_debug != 0) {
      printf("[HelixFitHack::doPatternRecognition]: calling rescueHitsBeforeSeed\n");
      printInfo(Helix);
    }
    rescueHitsBeforeSeed(Helix);
//-----------------------------------------------------------------------------
// finally, assume that the found helix is already close enough and refine  
// the helix parameters accounting for different weights
//-----------------------------------------------------------------------------
    if (_debug != 0)  printInfo(Helix);
    banner="refineHelixParameters";
    refineHelParamRes = refineHelixParameters(Helix, 0, _indicesTrkCandidate, _debug, banner);
    if ( refineHelParamRes >= 0){
      Helix._center.setX(Helix._sxyw.x0());
      Helix._center.setY(Helix._sxyw.y0());
      Helix._radius    = Helix._sxyw.radius();
      Helix._sxy.init(Helix._sxyw);
      if (_debug != 0)  printInfo(Helix);
    }
//---------------------------------------------------------------------------------------
// use the results of the helix search to see if points along the track can be rescued
//---------------------------------------------------------------------------------------
    if (Helix._srphi.qn() == 0) {
      usePhiResid = 0;
    }else {
      usePhiResid = 1;
    }
    rescueHits(Helix, 0, _indicesTrkCandidate, usePhiResid);

    if (_debug != 0)  printInfo(Helix);
//--------------------------------------------------------------------------------------------------------------
// 2015-03-25 G. Pezzu added the following call to findDfDz(...) in order to help the fitter on finding
// the more reliable value of dfdz which is needed for resolving the 2pi ambiguity.
// Since in the previous step we could have rescued few points, that would give us an help!
//--------------------------------------------------------------------------------------------------------------
                                                    //re-evaluate the df/dz and phi0 using also new hits rescued
                                                    // and new helix parameters
    rs = findDfDz(Helix, 0, _indicesTrkCandidate);

    if (rs == 1) { //update Helix values
      Helix._dfdz = _hdfdz;
      Helix._fz0  = _hphi0;
    }

    rc = doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate, useInteligentWeight);
    //    rc = doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate);


    if (rc){
      usePhiResid = 1;
      rescueHits(Helix, 0, _indicesTrkCandidate, usePhiResid);

      if (_debug != 0)  printInfo(Helix);
      banner="refineHelixParameters-after-doLinearFitPhiZ";
      rc1 = refineHelixParameters(Helix, 0, _indicesTrkCandidate, _debug, banner);
      if (rc1 >=0 ){
	Helix._center.setX(Helix._sxyw.x0());
	Helix._center.setY(Helix._sxyw.y0());
	Helix._radius    = Helix._sxyw.radius();
	Helix._sxy.init(Helix._sxyw);
	if (_debug != 0)  printInfo(Helix);
      }
    }
    
    // 2014-11-09 gianipez changed the cleanup process. now it is faster and cleaner
    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: calling filterUsingPatternRecognition\n");

                                  //update hack data with last results
    THackData* hack;
    hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
    hack->fData[14] = Helix._radius;
    hack->fData[15] = Helix._sxy.chi2DofCircle();
    hack->fData[6]  = Helix._fz0;
    hack->fData[7]  = Helix._dfdz*Helix._radius;
    hack->fData[8]  = Helix._dfdz;

    filterUsingPatternRecognition(Helix);

  PATTERN_RECOGNITION_END:;
//-----------------------------------------------------------------------------
// if running in the diagnostics mode, save state of the Xyzp (this is a deep copy)
//-----------------------------------------------------------------------------
    if (_debug != 0) printf("[HelixFitHack::doPatternRecognition]: END\n");
  }


//--------------------------------------------------------------------------------
// define the function used for projecting the strawhit error along the radial 
// direction of the helix-circle
//--------------------------------------------------------------------------------
  double  HelixFitHack::calculateWeight(Hep3Vector& HitPos   ,
					Hep3Vector& StrawDir ,
					Hep3Vector& HelCenter,
					double      Radius   ,
					int         Print    ,
					const char* Banner   ) {

    double    rs(2.5);   // straw radius, mm
    double    ew(30.0);  // assumed resolution along the wire, mm

    double x  = HitPos.x();
    double y  = HitPos.y();
    double dx = x-HelCenter.x();
    double dy = y-HelCenter.y();

    double costh  = (dx*StrawDir.x()+dy*StrawDir.y())/sqrt(dx*dx+dy*dy);
    double sinth2 = 1-costh*costh;

    double e2     = ew*ew*sinth2+rs*rs*costh*costh;
    double wt     = 1./e2;
                                                    //scale the weight for having chi2/ndof distribution peaking at 1
    //    wt /= 1.57;
    //    wt /= 2.83;
    wt *= _weightXY;

    if (Print > 0) {
      double dr = calculateRadialDist(HitPos,HelCenter,Radius);
      printf("[HelixFitHack::%s] %10.3f %10.3f %10.5f %10.5f %10.5f %10.5f %12.5e %10.3f\n",
	     Banner,x,y,dx,dy,costh,sinth2,e2,dr);
    }
    
    return wt;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  double  HelixFitHack::calculatePhiWeight(Hep3Vector HitPos   ,
					   Hep3Vector StrawDir ,
					   Hep3Vector HelCenter,
					   double     Radius   ,
					   int        Print    ,
					   TString    Banner   ) {
    double    rs( 2.5);  // mm
    double    ew(30.0);  // mm - erro along the wire   double x  = HitPos.x();

    double x  = HitPos.x();
    double y  = HitPos.y();
    double dx = x-HelCenter.x();
    double dy = y-HelCenter.y();

    double costh  = (dx*StrawDir.x()+dy*StrawDir.y())/sqrt(dx*dx+dy*dy);
    double sinth2 = 1-costh*costh;

    double e2     = ew*ew*costh*costh+rs*rs*sinth2;
    double wt     = Radius*Radius/e2;
    wt           *= _weightZPhi;

    if (Print > 0) {
      double dr = calculateRadialDist( HitPos, HelCenter, Radius);
      printf("[HelixFitHack::%s] %10.3f %10.3f %10.5f %10.5f %10.5f %10.5f %12.5e %10.3f\n", Banner.Data(), x, y, dx, dy, costh, sinth2, e2, dr);
    }
    
    return wt;
  }

//--------------------------------------------------------------------------------
// calculate the radial distance of a srtawhit form the helix prediction
//--------------------------------------------------------------------------------
  double  HelixFitHack::calculateRadialDist (const Hep3Vector& HitPos   , 
					     const Hep3Vector& HelCenter, 
					     double            Radius   ) {
    double dx = HitPos.x()-HelCenter.x();
    double dy = HitPos.y()-HelCenter.y();
    double dr = sqrt(dx*dx+dy*dy)-Radius;
  
    return dr;
  }


//-----------------------------------------------------------------------------
  void   HelixFitHack::doWeightedCircleFit (::LsqSums4&     TrkSxy   , 
					    int             SeedIndex, 
					    int*            IdVec    ,
					    Hep3Vector&     HelCenter, 
					    double&         Radius   , 
					    double*         Weights  ,
					    int             Print    , 
					    TString         Banner   ) {
    Hep3Vector hitPos, strawDir;
    double     wt;
    int        np = _xyzp.size();

    //-----------------------------------------------------------------------------
    // add cluster with a position error of 10 mm => wt = 1/100
    //-----------------------------------------------------------------------------
    TrkSxy.clear();
    TrkSxy.addPoint(fTimePeak->ClusterX(), fTimePeak->ClusterY(), 1./100.);

    //-------------------------------------------------------------------------------
    // add stopping target center with a position error of 100 mm / sqrt(12) ~ 30mm
    // so wt = 1/900
    //-------------------------------------------------------------------------------
    TrkSxy.addPoint(0., 0., 1./900.);


    for (int i=SeedIndex; i<np; i++) {
      if ( _xyzp[i].isOutlier())           goto NEXT_POINT;
   
      hitPos     = _xyzp[i]._pos;
      strawDir   = _xyzp[i]._sdir;

      wt         = calculateWeight(hitPos, strawDir, HelCenter, Radius, Print > 0 ? IdVec[i] : 0, Banner);
      Weights[i] = wt;

      if ( IdVec[i] < 1)                  goto NEXT_POINT;
 
      TrkSxy.addPoint(hitPos.x(), hitPos.y(), Weights[i]);

    NEXT_POINT: ;
    }

    //update helix info
    Radius  = TrkSxy.radius();
    HelCenter.setX( TrkSxy.x0());
    HelCenter.setY( TrkSxy.y0());
  }


//-----------------------------------------------------------------------------
  void    HelixFitHack::searchWorstHitWeightedCircleFit(int             SeedIndex,
							int*            IdVec,
							Hep3Vector&     HelCenter, 
							double&         Radius, 
							double*         Weights,
							int&            Iworst, 
							double&         HitChi2Worst) 
  {
    HitChi2Worst  = fHitChi2Max;
    Iworst        = -1;

    int        np = _xyzp.size();
    double     wt, e2, dr, hitChi2;
    Hep3Vector hitPos;

    for (int i=SeedIndex; i<np; i++) {
      if (_xyzp[i].isOutlier())           continue;

      wt = Weights[i];
      e2 = 1./wt;

      hitPos = _xyzp[i]._pos;
      dr = calculateRadialDist( hitPos, HelCenter, Radius);

      hitChi2 = dr*dr/e2; 
    
      // store info aout the radial residual 
      if (SeedIndex == 0){
	_distTrkCandidate[i] = hitChi2;
      }
					// avoid the use of hit rejected by the helix search
      if (IdVec[i] < 1) continue;

      if (hitChi2 > HitChi2Worst) {
	HitChi2Worst  = hitChi2;
	Iworst  = i;
      }
      
    }
    
  }

//--------------------------------------------------------------------------------
// returns the index of the hit which provides the highest contribute to the chi2
//--------------------------------------------------------------------------------
void    HelixFitHack::doCleanUpWeightedCircleFit(::LsqSums4&     TrkSxy, 
						 int             SeedIndex,
						 int*            IdVec,
						 Hep3Vector&     HelCenter, 
						 double&         Radius, 
						 double*         Weights,
						 int&            Iworst)
{
    Iworst     = -1;

    ::LsqSums4 sxy;
    double     chi2_min = 1e6;
    int        np = _xyzp.size();
    double     wt, chi2, x, y;
    
    for (int i=SeedIndex; i<np; i++) {
      if (_xyzp[i].isOutlier())           goto NEXT_P;

      // avoid the use of hit rejected by the helix search
      if ( IdVec[i] < 1)                  goto NEXT_P;

      sxy.init(TrkSxy);

      x  = _xyzp[i]._pos.x();
      y  = _xyzp[i]._pos.y();
 
      wt = Weights[i];
      sxy.removePoint(x, y, wt);
    
      chi2  = sxy.chi2DofCircle();

      if (chi2 < chi2_min){
	chi2_min    = chi2;
	Iworst      = i;
      }
    NEXT_P: ;
    }
  }
  
//-----------------------------------------------------------------------------
// 1. where is the cluster - at this point it is no longer needed
// assume the 
//-----------------------------------------------------------------------------
  int HelixFitHack::refineHelixParameters(HelixFitHackResult& Trk,
					  int                 SeedIndex, 
					  int*                IndexVec, 
					  int                 Print, 
					  TString             Banner) {
    double     x, y, r;
    double     hitChi2Worst;

    ::LsqSums4 sxyw;
    int        pointsRemoved(0);
    int        iworst(-1);
    double     wtWorst;
    double     chi2, chi2_min;
    int        np = _xyzp.size();

    //crete a temporary vector housing the indexVec information
    int       idVec[np];
    for (int i=0; i<np; ++i){
      idVec[i] = IndexVec[i];
    }

    Hep3Vector hitPos, strawDir, helCenter;

    double weights[np];
    int    success = -1;
  
    //initialize helix info
    r  = Trk._sxy.radius();
    helCenter.setX( Trk._sxy.x0());
    helCenter.setY( Trk._sxy.y0());
    

    if (_debug > 5) {
      printf("[HelixFitHack::refineHelixParameters] starts x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	     Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle());
      printf("[HelixFitHack::refineHelixParameters] i       X        Y        dx        dy         costh        sinth2         e2     radial-dist\n");
    }
    
    doWeightedCircleFit (sxyw, SeedIndex, idVec,  helCenter,  r,  weights, Print, Banner);
 
    //now recalcute the weights using the most recent 
    //helix parameters
    //-----------------------------------------------------------------------------
    doWeightedCircleFit (sxyw, SeedIndex, idVec,  helCenter,  r,  weights, Print, Banner);
 
    searchWorstHitWeightedCircleFit(SeedIndex, idVec, helCenter, r, weights,
				    iworst, hitChi2Worst);
    
    chi2     = sxyw.chi2DofCircle();
    chi2_min = chi2;

    if ((chi2 <= _chi2xyMax) && ( hitChi2Worst  <= 25.)) {
      success = 0;
      goto F_END;		 
    }
    //-----------------------------------------------------------------------------
    // cleanup is needed
    //-----------------------------------------------------------------------------
    
    if (_debug > 5) {
      printf("[HelixFitHack::refineHelixParameters] cleaunup starts x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	     sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle());
    }

  NEXT_ITERATION:;
    doCleanUpWeightedCircleFit( sxyw, SeedIndex, idVec,
				helCenter,  r, weights, iworst);
    
    if (iworst >= 0){
      x       = _xyzp[iworst]._pos.x();
      y       = _xyzp[iworst]._pos.y();
      wtWorst = weights[iworst];
      
      //remove point from the track    
      sxyw.removePoint(x, y, wtWorst);
      
      if (_debug > 5) {
	printf("[HelixFitHack::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
	       sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle(), _chi2xyMax, iworst);
      
      }
      
      //mark point as outlier
      idVec[iworst] = 0;
      
     //update helix info
      r  = sxyw.radius();
      helCenter.setX( sxyw.x0());
      helCenter.setY( sxyw.y0());
  
      // now update helix 
      doWeightedCircleFit (sxyw, SeedIndex, idVec,  helCenter,  r,  weights, 0, Banner);
       
      //update the chi2 value
      chi2_min = sxyw.chi2DofCircle();
   
      ++pointsRemoved;
    }
    
    //-----------------------------------------------------------------------------
    // recalculate the worst radial residual
    //-----------------------------------------------------------------------------
  CHECK_RESIDUALS: ;

    searchWorstHitWeightedCircleFit(SeedIndex, idVec, helCenter, r, weights,
				    iworst, hitChi2Worst);
//-----------------------------------------------------------------------------
// if a hit contributes chi2 > 25, remove it and go back looking for the next such hit
//-----------------------------------------------------------------------------
    if (iworst >= 0) {
      x       = _xyzp[iworst]._pos.x();
      y       = _xyzp[iworst]._pos.y();
      wtWorst = weights[iworst];
     
      //remove point from the track    
      sxyw.removePoint(x, y, wtWorst);
      if (_debug > 5) {
	printf("[HelixFitHack::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
	       sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle(), _chi2xyMax, iworst);
      
      }
      //mark point as outlier
      idVec[iworst] = 0;
      
      //update helix info
      r  = sxyw.radius();
      helCenter.setX( sxyw.x0());
      helCenter.setY( sxyw.y0());

      // now update helix 
      doWeightedCircleFit (sxyw, SeedIndex, idVec,  helCenter,  r,  weights, 0, Banner);
     
      //update the chi2 value
      chi2_min = sxyw.chi2DofCircle();
   
      ++pointsRemoved;
      goto CHECK_RESIDUALS;
    }

   
    if ( ( chi2_min    >= _chi2xyMax ) && 
	 ( iworst >= 0 ) ){
      //-----------------------------------------------------------------------------
      // still bad chi2, repeat the cleanup cycle
      //-----------------------------------------------------------------------------
      if (sxyw.qn() > 10.) {
	goto NEXT_ITERATION;
      }
    }else if ( chi2_min < _chi2xyMax){
      success = 1;
    }
  
  F_END:;
    if (_debug > 5) {
      printf("[HelixFitHack::refineHelixParameters] success = %i\n", success);
      printf("[HelixFitHack::refineHelixParameters] points removed = %i  x0 = %5.3f y0 = %5.3f r = %5.3f chi2 = %5.5f\n", 
	     pointsRemoved, sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle());
    }

    if ( (success == 0) ||
	 (success == 1)){
    //-----------------------------------------------------------------------------
    // update circle parameters
    //-----------------------------------------------------------------------------
      Trk._sxyw.init(sxyw);
      Trk._cw.setX(Trk._sxyw.x0());
      Trk._cw.setY(Trk._sxyw.y0());
      Trk._rw    = Trk._sxyw.radius();
      Trk._chi2w = Trk._sxyw.chi2DofCircle();

      for (int i=0; i<np; ++i){
	IndexVec[i] = idVec[i];
      }
    }
   //  if (SeedIndex ==0){
//       THackData* hack;
//       hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
//       hack->fData[14] = Trk._rw ;
//       hack->fData[15] = Trk._chi2w ;
//     }
    
    return success;
  }


//-----------------------------------------------------------------------------
// The previous helix search may have  thrown away point which instead have
// small radial residual, so this function is devoted for rescueing these
//-----------------------------------------------------------------------------
  void HelixFitHack::rescueHits(HelixFitHackResult& Trk            ,
				int                 SeedIndex      , 
				int*                IndexVec       ,
				int                 UsePhiResiduals) 
  {
    const char  banner[] = "rescueHits";
    double      wt, e2, x, y, r;
    double      phiwt(-9999.);

    Hep3Vector  hitPos, strawDir, helCenter;

    int         np = _xyzp.size();
    double      weights[np];

    double      dfdz, phi0, dphi, dphiChi2(0.0), phi_pred;

    ::LsqSums4  sxy;
    int         n_added_points(0);
    int         ibest;
    double      wtBest, phiwtBest;
    double      chi2, chi2_min, dr, hitChi2, drChi2;

    //set  dfdz and phi0
    dfdz = Trk._dfdz;
    phi0 = Trk._fz0;

    //update helix info
    r  = Trk._sxy.radius();
    helCenter.setX( Trk._sxy.x0());
    helCenter.setY( Trk._sxy.y0());

    doWeightedCircleFit (Trk._sxy, SeedIndex, IndexVec,  helCenter,  r,  weights);
  
    if (_debug > 5) {
      printf("[HelixFitHack::%s] starts x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f phi0 = %9.6f dfdz = %9.6f chi2 = %8.3f\n",
	     banner,
	     Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(),
	     phi0, dfdz, Trk._srphi.chi2DofLine());
    }
  
      //now perform some clean up if needed

    chi2 = Trk._sxy.chi2DofCircle();
    if (chi2 >= _chi2xyMax) {
      if (_debug > 5) {
	printf("[HelixFitHack::%s] chi2 = %8.3f already at limit! no point can be added\n",banner,chi2);
      }
      goto F_END;		 
    }
    //-----------------------------------------------------------------------------
    // now add points
    //-----------------------------------------------------------------------------

    if (_debug > 5) {
      printf("[HelixFitHack::%s] x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	     banner,Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle());
      printf("[HelixFitHack::%s] i       X        Y        dx        dy         costh        sinth2         e2     radial-dist\n",
	     banner);
    }
  NEXT_ITERATION:;
    if (UsePhiResiduals == 1) {
      chi2_min    = 2.*fHitChi2Max;
    } else{
      chi2_min    = fHitChi2Max;
    }
    dphiChi2    = 0.;
    ibest       = -1;
    wtBest      = -1;
    phiwtBest   = -1;

    for (int i=SeedIndex; i<np; i++) {
      if (_xyzp[i].isOutlier())           continue ;

      // avoid the use of hit already set as active
      if ( IndexVec[i] == 1)             continue;

      hitPos    = _xyzp[i]._pos;
      strawDir  = _xyzp[i]._sdir;


      dr = calculateRadialDist(hitPos,helCenter,r);
      wt = calculateWeight    (hitPos,strawDir,helCenter,r,_debug,banner);

      drChi2  = (dr*dr)*wt;

      // 2015-03-25 G.Pezzu added the request of "distance" also for the phi
    
      if ((UsePhiResiduals == 1) && (fPhiCorrectedDefined)) {
	phi_pred = hitPos.z()*dfdz + phi0;
	dphi     = phi_pred - fPhiCorrected[i];
	phiwt    = calculatePhiWeight(hitPos, strawDir, helCenter, r, 0, banner);
	dphiChi2 = dphi*dphi*phiwt;                
	hitChi2  = drChi2 + dphiChi2;
      } else {
	hitChi2  = drChi2;
      }

      if (_debug > 5) {
	printf("[HelixFitHack::%s] sigmaphi2 = %10.3f drChi2 = %5.3f dphiChi2 = %5.3f chi2 = %5.3f\n",
	       banner, 1./phiwt, drChi2, dphiChi2, hitChi2);
      }
//-----------------------------------------------------------------------------
// require chi2 < fHitChi2Max, identify the closest point
//-----------------------------------------------------------------------------
      if ( (hitChi2 < chi2_min) && (drChi2 < fHitChi2Max) && (dphiChi2 < fHitChi2Max)) {

//-----------------------------------------------------------------------------
// check if XY-chi2 and ZPhi-chi2 are less than chi2xyMax and chi2zphiMax
//-----------------------------------------------------------------------------
	x = hitPos.x();
	y = hitPos.y();
	Trk._sxy.addPoint  (x, y, wt);
	Trk._srphi.addPoint(hitPos.z(), fPhiCorrected[i], phiwt);
	
	if (Trk._sxy.chi2DofCircle() < _chi2xyMax){
	  if (UsePhiResiduals == 1){
	    if (Trk._srphi.chi2DofLine() < _chi2zphiMax){
	      chi2_min    = hitChi2;
	      ibest       = i;
	      wtBest      = wt;
	      phiwtBest   = phiwt;
	    }
	  }else {
	    chi2_min    = hitChi2;
	    ibest       = i;
	    wtBest      = wt;
	  }
	}
	
	Trk._sxy.removePoint  (x, y, wt);
	Trk._srphi.removePoint(hitPos.z(), fPhiCorrected[i], phiwt);

      }
    }
    
    if (ibest >= 0){
      x  = _xyzp[ibest]._pos.x();
      y  = _xyzp[ibest]._pos.y();
                                       //add point from the track    
      Trk._sxy.addPoint(x, y, wtBest);
      
      if (UsePhiResiduals == 1){
	Trk._srphi.addPoint(_xyzp[ibest]._pos.z(), fPhiCorrected[ibest], phiwtBest);

      }
      
      if (_debug > 5) {
	printf("[HelixFitHack::%s] x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %6.3f chi2Maxxy = %6.3f index point added = %i straw-id = %6i hitChi2 = %6.3f x = %8.3f y = %8.3f z = %9.3f \n ",
	       banner,
	       Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(), _chi2xyMax, ibest, 
	       _xyzp[ibest]._strawhit->strawIndex().asInt(), chi2_min, 
	       x, y, _xyzp[ibest]._pos.z());
      }

      //mark point as active
      IndexVec[ibest] = 1;

      r  = Trk._sxy.radius();
      helCenter.setX( Trk._sxy.x0());
      helCenter.setY( Trk._sxy.y0());
                                        // now update helix 
      doWeightedCircleFit (Trk._sxy, SeedIndex, IndexVec,  helCenter,  r,  weights);
  
      ++n_added_points;
	                              goto NEXT_ITERATION;
    }
   
    //----------------------------------------------------------------------
    //now update information about the radial residual of the hits
    //----------------------------------------------------------------------
    for (int i=SeedIndex; i<np; i++) {
      if (_xyzp[i].isOutlier())           continue;
      
      wt      = weights[i];
      e2      = 1./wt;
      hitPos  = _xyzp[i]._pos;
      dr      = calculateRadialDist( hitPos, helCenter, r);
      hitChi2 = dr*dr/e2; 
      
      // store info aout the radial residual 
      if (SeedIndex == 0){
	_distTrkCandidate[i] = hitChi2;
      }
    }
    
 //-----------------------------------------------------------------------------
 // update circle parameters
 //-----------------------------------------------------------------------------
    Trk._center.setX(Trk._sxy.x0());
    Trk._center.setY(Trk._sxy.y0());
    Trk._radius  = Trk._sxy.radius();
    Trk._chi2    = Trk._sxy.chi2DofCircle();

  F_END:;
    if (_debug > 5 ) {
      printf("[HelixFitHack::%s] points added = %i chi2 = %5.5f\n",banner,n_added_points,Trk._chi2);
    }
   
    if (SeedIndex ==0){
      THackData* hack;
      hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
      //      hack->fData[14] = Trk._radius ;
      //hack->fData[15] = Trk._chi2 ;
      hack->fData[16] += n_added_points;  
    }
  }



//-----------------------------------------------------------------------------
  void HelixFitHack::findTrack(int                 SeedIndex          , 
			       double&             Chi2               , 
			       int&                CountGoodPoints    ,
			       HelixFitHackResult& Helix              ,
			       int&                Mode               , 
			       bool                UseDefaultDfDz     ,
			       int                 UseMPVDfDz         ) {
//-----------------------------------------------------------------------------
// Legend of the points used for performing the pattern recongition:
// ----------------------------------------------------------------
// p0 : center of the helix, assume p0.z = 0
// p1 : center of the stopping target
// p2 : point in the position SeedIndex on the vector Xyzp
// p3 : postion of the EMC cluster
//-----------------------------------------------------------------------------
    double     radius, phi0, tanLambda;
    double     dx,dy,phi,dx2, dy2;
    Hep3Vector p0, p1, p2, p3; 
    Hep3Vector shPos, hePos;
//-----------------------------------------------------------------------------
// initialize output paramters
//-----------------------------------------------------------------------------
    CountGoodPoints = 0;
    Mode            = 0;

    p1 = Hep3Vector(0., 0., 0.);   //5971. - 10200.
    p2 = _xyzp[SeedIndex]._pos;
    p3 = Hep3Vector(fTimePeak->ClusterX(),fTimePeak->ClusterY(), fTimePeak->ClusterZ());

    //    double radius,phi0,tanLambda;
    
    //----------------------------------------------------------------------//
    //index for setting the sostitute-point for the target center
    int np = _xyzp.size();
    int  mode0GoodPoints(0), mode1GoodPoints(0); //mode0GoodPoints is the numebr of points belonging to a trajectory
                                                 //when dfdz is not re-calculated using the function calculateDfDz()
                                                 //mode1GoodPoints is the numebr of point belonging to a trajectory 
                                                 //when dfdz is re-computed using calculateDfDz()
    int  rescuedPoints(0);
//-----------------------------------------------------------------------------
// initialize a vector to store the straw hit candidate for dfdz recalculation
// create a temporary array for storing indices of the good point which belong 
// to a track candidate
//-----------------------------------------------------------------------------
    int candidateList[np];
    int markIndexList[np];

    double dzList  [np];
    double distList[np];

    for(int i=0; i<np; ++i){
      candidateList[i] = -9999;
      markIndexList[i] = -9999;
      dzList[i]        = -9999;
      distList[i]      = -9999;
    }
						
    markIndexList[SeedIndex] = 1;
//---------------------------------------------------------------------
// define constrains on the z coordinate of the strawhit candidate 
// for re-calculating dfdz
// if the candidate and the seeding strawhit are too close, the dfdz 
// calculated could be affected by several effects which lead on the 
// wrong estimation.
// We are asking that the candidate straw must be at a distance along 
// the z axes greateer than tollMin and less than tollMax. 
// These parameters still need to be optimized
//-----------------------------------------------------------------------------
    double tollMin(100.), tollMax(500.);


    //2014-03-10 gianipez changed the values of the following tollerance
    // for X and Y distance to avoid delta electron for corrupting the pattern-reco
    //     double tollXYdist(100.);//require a dist in the transverse plane > 10 cm
    

    int goodPoint(-999); // index of the strawhit candidate for dfdz and helix paramters recalculation
    //2014-01-29 gianipez added the followign line

    //few paramterss used for calculating strawhit position residuals
    double weight(0.0);
    double deltaZ(0.), deltaZfromSeed(0.), deltaX(0.), deltaY(0.);
    double distXY(0.0);
    double dist(0.0), dist2(0.0); //help parameter for storing strawhit position residual

  //   double distGoodPoint(_distPatRec);             //threshold on the squared distance of a strawhit position with 
//                                                    //respect to the predicted position


    double z0,z1,phi_0,phi_1;
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
    double dfdz = _mpDfDz;//tanLambda/radius;
    //    double chi2min(1e10);
    double dfdz_end, phi0_end, /*phi0_seed,*/ radius_end;
    //    double x0_end, y0_end;
    //    int    mode_end(-1);
    

					// two flags are needed:
    bool removeTarget(true);            // avoid the recalculation of dfdz 
					// and helix parameters in case when
                                        // others strawhit candidates are found 

    bool isStored(false);               // help parameter for indicating if 
					// a straw hit has been already used 
                                        // for dfdz and helix recalculation
    Chi2 = 0.0;
//----------------------------------------------------------------------
// calculate helix paramters using the center of the stopping target,
// the EMC cluster which seeded the CalTimePeak and the seeding strawhit.
// The z coordinate of the target center is set to 0 because in the formula
// inside calculateTrackParameters(...) is not used its z coordinate
//-----------------------------------------------------------------------------
    p1 = Hep3Vector(0., 0., 0.);
    calculateTrackParameters(p0,radius,phi0,tanLambda,
			     p1,p2,p3,
			     Helix,
			     false);
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
    if (UseMPVDfDz ==1 ){
      dfdz = _hdfdz;//_mpDfDz;
    }

    fLastIndex = -9999;
    //--------------------------------------------------//

    ::LsqSums4 sxy;
    ::LsqSums4 srphi;

    sxy.addPoint(p2.x(), p2.y());         //seeding strawhit
    sxy.addPoint(p3.x(), p3.y());         //EMC cluster position
    sxy.addPoint(0., 0., 0.1);                 //Target center in the transverse plane

    TString banner="HelixFitHack::findTrack";
    if (UseMPVDfDz ==1) banner += "-UseMPVDfDz";

    double dz_max;
    TString name;
    name =  banner;
    name += "-loop";
    
    int i_last = SeedIndex;
      
    //define two paramters used for storing the index
    //of the straw hit closest to z=0
    int     i_z0(-1);
    double  dz_z0(1e10);
  
    for (int i=SeedIndex+1; i<np; i++) {
      name =  banner;
      name += "-loop";
      if (_xyzp[i].isOutlier()) goto NEXT_POINT;
      weight = 1.;
//----------------------------------------------------------------------//
// 2014-12-26 Gianipez added the request that the hit has not already 
// been used by a previous search
//-----------------------------------------------------------------------------
      if (isHitUsed(i) == 1) {
	if( _debug > 10){
	  //	  printf("[HelixFitHack::findTrack-loop]  XYZP-hit number = %i skipped\n", i);
	  printf("[%s]  XYZP-hit number = %i skipped\n", name.Data(), i);
	}
	goto NEXT_POINT;
      }
      shPos = _xyzp[i]._pos;
    
      deltaZ = shPos.z() - _xyzp[i_last]._pos.z();

      deltaZfromSeed = shPos.z() - p2.z();

      phi    = phi0 + (deltaZ)*dfdz;//tanLambda/radius;

      hePos  = Hep3Vector(p0.x() + radius*std::cos(phi),  // predicted hit position
			  p0.y() + radius*std::sin(phi),
			  shPos.z());
					                  // calculate residuals in XY
      dx  = hePos.x() - shPos.x();
      dx2 = dx*dx;
      dy  = hePos.y() - shPos.y();
      dy2 = dy*dy;
      
      dist2 = dx2 + dy2;
      dist  = std::sqrt(dist2);

      //calculate the distance of the straw hit point to the seeding strawhit on the transverse plane
      deltaY = std::fabs(p2.y() - shPos.y());
      deltaX = std::fabs(p2.x() - shPos.x());
      distXY = std::sqrt(deltaY*deltaY + deltaX*deltaX);
   
      if( _debug > 10){
	if( i==SeedIndex+1) {
	  printf("[%s]  findTrack() starts with helix parameters derived from these points \n", name.Data());
	  printf("[%s]   point  type      X         Y         Z       xyzp-index \n", name.Data());
	  printf("[%s] ----------------------------------------------------------\n", name.Data());
	  printf("[%s]    seeding      %10.3f   %10.3f   %10.3f   %8i \n", name.Data(),p2.x(), p2.y(), p2.z(), 
		 SeedIndex);
	  printf("[%s]   candidate     %10.3f   %10.3f   %10.3f   %8i \n", name.Data(),p1.x(), p1.y(), p1.z(), 
		 fLastIndex);
	  printf("[%s]  emc cluster    %10.3f   %10.3f   %10.3f   %8i \n", name.Data(),p3.x(), p3.y(), p3.z(), 
		 -1);
	}
	
	printf("[%s] X0 = %10.3f Y0 = %10.3f r = %10.3f dfdz = %5.5f \n", name.Data(), 
	       p0.x(), p0.y(), radius, dfdz);
	printf("[%s]   measured      %10.3f %10.3f %10.3f %8i \n", name.Data(), 
	       shPos.x(), shPos.y(), shPos.z(), i);
	printf("[%s]  predicted      %10.3f %10.3f %10.3f %8i \n", name.Data(), 
	       hePos.x(), hePos.y(), hePos.z(), i);
	printf("[%s] dist-from-prediction = %10.3f  dist-from-seedXY = %5.3f dz-from-seed = %5.3f\n", name.Data(), 
	       dist, distXY, deltaZ);
	printf("[%s] ---------------------------------------------------------\n", name.Data());

      }
//-----------------------------------------------------------------------------
// max_dist: running search window accounts for the finite extrapolation accuracy
//-----------------------------------------------------------------------------
      isStored = false;
      dz_max = _distPatRec + _dfdzErr*deltaZ;
      if ( dist <= dz_max ){
	++CountGoodPoints;

	// 2014-11-12 gianipez:
	//adjust the helix-center coordinates and the radius
	// if the dfdz value has already been evaluated
	if (Mode == 1){
	  p0.setX( sxy.x0());
	  p0.setY( sxy.y0());
	  radius  = sxy.radius();
	}
					            // index of the last good point
	i_last = i; 

	if (fabs(shPos.z()) < dz_z0){
	  dz_z0 = fabs(shPos.z());
	  i_z0  = i;
	}

	//2015-01-27 G. Pezzu and P. Murat
	phi0   = CLHEP::Hep3Vector(shPos - p0).phi();

	//add point to the helixfithack result objet
	sxy.addPoint(shPos.x(),shPos.y(), weight);

	//store the index of the good point found
	markIndexList[i] = 1;
					            // Z-distance from the last point found relying on thr helix
	dzList[i]        = deltaZ;
					            // distance from predition
	distList[i]      = dist;

	if (Mode == 0) {
	  ++mode0GoodPoints;	
	}
	else if ((Mode == 1) && (i<= fLastIndex)){
	  ++mode1GoodPoints;
	}

	for(int j=0; j<np; ++j){
	  if(candidateList[j] == i)
	    isStored = true;
	}
	  
	if(!isStored){
	  if( ( deltaZfromSeed > tollMin ) &&
	      ( deltaZfromSeed < tollMax )){	    
	    if (removeTarget) goodPoint = i;
	  }
	}
      } else {
//-----------------------------------------------------------------------------
// hit outside the search road
//-----------------------------------------------------------------------------
	markIndexList[i] = 0;
	distList[i]      = 0;
	dzList[i]        = 0;
      }

      // 2014-04-23     gianipez fixed a bug
      if ( CountGoodPoints >= 2 &&
	   removeTarget         &&
	   (goodPoint >=0)      &&
	   (goodPoint != fLastIndex) ) {

//recalculate helix parameters using the strawhit candidate "goodPoint"
	p1 = _xyzp[goodPoint]._pos;

	p0.setX( sxy.x0());
	p0.setY( sxy.y0());
	radius  = sxy.radius();

//now calculate more accuratelly the value of dfdz using just the two strawhit positions
	z0 = p2.z();//z coordinate of the seeding point
	z1 = p1.z();
	phi_0 = CLHEP::Hep3Vector(p2 - p0).phi();
	phi_1 = CLHEP::Hep3Vector(p1 - p0).phi();

//2015-01-14 G. Pezzullo added the following condition because in case 
// we have a MPV for dfdz from the procedure findDfDZ we want just to use it
	if (UseMPVDfDz == 0){
	  calculateDfDz(phi_0, phi_1, z0, z1, dfdz);
	} 
	else if (UseMPVDfDz ==1) {
	  dfdz = _hdfdz;                   //_mpDfDz;
	}

	name = banner;
	name += "2strawhitsHelixDef";
	//	printf(" HelixFitHack::_debug: %5i TEEEEST\n",_debug);
	if (_debug > 10) {
	  printf("[%s] strawhit type     X        Y        Z     index\n", name.Data());
	  printf("[%s] ----------------------------------------------------\n", name.Data());
	  printf("[%s]    seeding     %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p2.x(), p2.y(), p2.z(), SeedIndex);
	  printf("[%s]   candidate    %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p1.x(), p1.y(), p1.z(), goodPoint);
	  printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f dfdz = %5.6f chi2 = %5.3f \n", name.Data(),
		 p0.x(), p0.y(), radius, dfdz , sxy.chi2DofCircle());
	
	}
//-----------------------------------------------------------------------------
// what to do if dfdz s negative?
//-----------------------------------------------------------------------------
	if ((dfdz > _maxDfDz) || (dfdz < _minDfDz)) {

	  for(int j=0; j<np; ++j){
	    if(candidateList[j]<0){
	      candidateList[j] = goodPoint;
	      break;
	    }
	  }

	  if (_debug > 10) {
	    printf("[%s] dfdz = %5.5f not in range limits. Continue the search\n",name.Data(),dfdz);
	  }
	  p1 = Hep3Vector(0., 0., 0.);
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
	  dfdz = _mpDfDz;
				// 	    dfdz = (tanLambda/radius)*double(j);
	}else{
	  removeTarget = false;
	  Mode         = 1;
	  fLastIndex   = goodPoint;
	}
	//------------------------------------------------------------//
      }
    NEXT_POINT:;
    }

    if (CountGoodPoints < 3) return;

    //declare a temporary variable for storing info about the dfdz
    //values out of the method 'calculateDfDz(...)'
    double dfdzRes  [3] = {   -1.,    -1.,    -1.};
    double dphi0Res [3] = {-9999., -9999., -9999.};
    double radiusRes[2] = {   -1.,    -1.};

    shPos        = _xyzp[i_z0]._pos;
    phi0         = CLHEP::Hep3Vector(shPos - p0).phi();
    
    if (UseMPVDfDz == 0){
      dfdzRes[0] = dfdz;
    }
    
    dphi0Res [0] = phi0 - shPos.z()*dfdz;
    radiusRes[0] = sxy.radius();
  
    name = banner;
    name += "-results";
   

    Chi2 = sxy.chi2DofCircle();
//-----------------------------------------------------------------------------
// 2015-01-22 G. Pezzullo and P. Murat; update the dfdz value using all hits
//-----------------------------------------------------------------------------
    HelixFitHackResult tmp1HelFitRes(Helix);
    HelixFitHackResult tmp2HelFitRes(Helix);
   
    //2015-01-27 G. Pezzu and P. Murat: initialize only the xy part, z-phi part is not needed here

    tmp1HelFitRes._sxy.init(sxy);
    tmp1HelFitRes._radius = sxy.radius();
    tmp1HelFitRes._center.set(sxy.x0(), sxy.y0(), 0.0);

    radius_end = sxy.radius();

    tmp2HelFitRes._center.set(p0.x(), p0.y(), 0.0);
    tmp2HelFitRes._radius = radius;
    tmp2HelFitRes._dfdz   = dfdz;
    
    int rc0 = refineHelixParameters(tmp1HelFitRes, SeedIndex, markIndexList);
    if ( rc0 >=0){
      tmp2HelFitRes._center.set(tmp1HelFitRes._cw.x(), tmp1HelFitRes._cw.y(), 0.0);
      tmp2HelFitRes._radius = tmp1HelFitRes._rw;
      radius_end            = tmp1HelFitRes._rw;
      sxy.init(tmp1HelFitRes._sxyw);
					// update the Chi2 value
      Chi2 = sxy.chi2DofCircle();
      CountGoodPoints = 0;
                                      //store information for hackdata
      radiusRes[1] = tmp2HelFitRes._radius;

      for (int i=SeedIndex; i<np; ++i){
	if ( markIndexList[i] > 0) ++CountGoodPoints;
      }
    }
    
    //2015-03-26 G. Pezzu added a condition for updating or not the dfdz value
    int rs = findDfDz(tmp2HelFitRes, SeedIndex, markIndexList);
    if (rs ==1 ){
      tmp2HelFitRes._dfdz = _hdfdz;
                                         //fill diag vector
      dfdzRes[1]  = _hdfdz;
      dphi0Res[1] = _hphi0;
    }
    tmp2HelFitRes._fz0    = _hphi0;

//-----------------------------------------------------------------------------
// 2015-01-23 G. Pezzu and P. Murat: when it fails, doLinearFitPhiZ returns negative value 
//                                   in this case, use the previous value for dfdz and phi0
//-----------------------------------------------------------------------------
    bool rcPhiZ = doLinearFitPhiZ(tmp2HelFitRes, SeedIndex, markIndexList);

    if (rcPhiZ) {
      dfdz_end   = tmp2HelFitRes._dfdz;
      phi0_end   = tmp2HelFitRes._fz0;
      srphi.init(tmp2HelFitRes._srphi);
                                         // fill diag vector
      dfdzRes [2] = tmp2HelFitRes._dfdz;
      dphi0Res[2] = tmp2HelFitRes._fz0;

      CountGoodPoints = 0;
      for (int i=SeedIndex; i<np; ++i) {
	if ( markIndexList[i] > 0) 
	  ++CountGoodPoints;
      }
       
    }else {
      dfdz_end = _hdfdz;
      phi0_end = _hphi0;
    }

     if (_debug > 10) {
      printf("[%s] strawhit type     X        Y        Z     index\n", name.Data());
      printf("[%s] ----------------------------------------------------\n", name.Data());
      printf("[%s]    seeding     %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p2.x(), p2.y(), p2.z(), SeedIndex);
      printf("[%s]   candidate    %5.3f  %5.3f  %5.3f   %i  \n", name.Data(),p1.x(), p1.y(), p1.z(), goodPoint);
      printf("[%s]  emc cluster   %5.3f  %5.3f  %5.3f \n", name.Data(),p3.x(), p3.y(), p3.z());
      printf("[%s] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.6fdfdz = %5.6f chi2 = %5.3f \n", name.Data(),
	     sxy.x0(), sxy.y0(), radius_end, phi0_end, dfdz_end , sxy.chi2DofCircle());
      printf("[%s] CountGoodPoints = %i\n", name.Data(), CountGoodPoints);
    }

    if (mode1GoodPoints>0){
      rescuedPoints = mode1GoodPoints - mode0GoodPoints ;
    } else{
      rescuedPoints = -1;
    }

//----------------------------------------------------------------------
    if (CountGoodPoints > 2) {
      //2014-01-29 gianipez added the following line
      Chi2       = sxy.chi2DofCircle();
      name = banner;
      name += "EndFindTrack";
      if (_debug > 10) {
	printf("[%s] Chi2 = %5.3f nGoodPoints = %d dfdz = %5.5f mode = %i\n",
	       name.Data(),
	       sxy.chi2DofCircle(),
	       CountGoodPoints,
	       dfdz_end,
	       Mode);
      }
    }
    
    if ( (Mode == 1) ||
	 UseDefaultDfDz ){
      if ( (CountGoodPoints > _goodPointsTrkCandidate) ||
	   ( (CountGoodPoints == _goodPointsTrkCandidate) && 
	     (Chi2             < _chi2TrkCandidate      ) ) ){
	//      update trackcandidate informations
	fSeedIndex      = SeedIndex;
	if (Mode == 1) {
	  fUseDefaultDfDz = 1;
	  fCandidateIndex = fLastIndex;
	}else if (UseDefaultDfDz){
	  fCandidateIndex = -9999;
	}
//----------------------------------------------------------------------
// 2015 - 01 - 17 G. Pezzu: remove the target center from sxy 
// in order to evaluate more accuratelly the helix parameters
//----------------------------------------------------------------------
	_x0     = sxy.x0();	// p0.x();
	_y0     = sxy.y0();	// p0.y();
	_phi0   = phi0_end;

	_radius = radius_end;//sxy.radius();	// radius_end;
	_dfdz   = dfdz_end;

	_goodPointsTrkCandidate = CountGoodPoints;
	_chi2TrkCandidate       = Chi2;
//-----------------------------------------------------------------------------
// reset the vector holding the informations abot:
// -> hit belonging to the track candidate
// -> distance in the X-Y plane from the prediction
// -> distance from the seeding hit along the z-axes
//-----------------------------------------------------------------------------
	for(int i=0; i<400; ++i){
	  _indicesTrkCandidate[i] = -9999;
	  _distTrkCandidate[i]    = -9999;
	  _dzTrkCandidate[i]      = -9999;
	}
	
	for (int i=SeedIndex; i<np; ++i){
	  _indicesTrkCandidate[i] = markIndexList[i];
	  _distTrkCandidate[i]    = distList[i];
	  _dzTrkCandidate[i]      = dzList[i];
	}

	Helix._center.set(_x0, _y0, 0.0);
	Helix._radius = _radius;
//-----------------------------------------------------------------------------
//now calculate the phi coordinate at z = 0
//-----------------------------------------------------------------------------
	Helix._fz0    = _phi0;// + _dfdz*(0.0 - p2.z()) ;
	Helix._dfdz   = _dfdz;
	Helix._sxy.init(sxy);
	Helix._srphi.init(srphi);
//-----------------------------------------------------------------------------
// fill diagnostics information for histogramming
//-----------------------------------------------------------------------------
	THackData* hack;
	hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
	int loopId(0);
	if (UseDefaultDfDz == 0) {
	  if (UseMPVDfDz) {
	    loopId = 2;
	  }else {
	    loopId = 0;
	  }
	}else {
	  loopId = 1;
	}
	hack->fData[4]  = loopId;
	hack->fData[5]  = _radius;
// 	hack->fData[6]  = phi0_end;
// 	hack->fData[7]  = dfdz_end*_radius;
// 	hack->fData[8]  = dfdz_end;
	hack->fData[9]  = rescuedPoints;

	double dz = p1.z() - p2.z();

	hack->fData[10] = Mode == 1 ? dz : -1.;
	hack->fData[11] = CountGoodPoints;
	hack->fData[12] = sxy.chi2DofCircle();
	hack->fData[13] = srphi.chi2DofLine();

	hack->fData[17] = dfdzRes[0];
	hack->fData[18] = dfdzRes[1];
	hack->fData[19] = dfdzRes[2];

	hack->fData[22] = dphi0Res[0];
	hack->fData[23] = dphi0Res[1];
	hack->fData[24] = dphi0Res[2];


	hack->fData[20] = radiusRes[0];
	hack->fData[21] = radiusRes[1];

	int j=0;
	for (int i=SeedIndex; i<np; ++i){
	  if (_indicesTrkCandidate[i] != 1) continue;
	  hack->fDist[j] = _distTrkCandidate[i];
	  hack->fDz[j]   = _dzTrkCandidate[i];
	  ++j;
	}	
      }
    }
  }

//-----------------------------------------------------------------------------
// 2015-03-15 P.Murat: 'Helix' doesn't seem to be necessary ?
//----------------------------------------------------------------------------- 
  void HelixFitHack::calculateTrackParameters(Hep3Vector&         p0, 
					      double&             radius,
					      double&             phi0, 
					      double&             tanLambda,
					      Hep3Vector          p1, 
					      Hep3Vector          p2,
					      Hep3Vector          p3,
					      HelixFitHackResult& Helix,
					      bool                cleanPattern) 
  {
    p0.setZ(p2.z()); 
   
    double x_m,y_m, x_n, y_n;
    //coordinates of the mean point between p1 and p3
    x_m = (p3.x() + p1.x())/2.;
    y_m = (p3.y() + p1.y())/2.;
    //------------------------------------------------------------//

    //coordianates of the mean point between p2 and p3    
    x_n = (p3.x() + p2.x())/2.;
    y_n = (p3.y() + p2.y())/2.;
    //------------------------------------------------------------//

    //calculate now the term of the line ortoghonal to the mid point of
    //the cord which links p1 and p3
    double m = -1.*(p3.x() - p1.x())/(p3.y() - p1.y());
    double c = y_m - x_m*m;
    //the eq. is: y = x*m + c

    //calculate now the term of the line ortoghonal to the mid point of
    //the cord which links p2 and p3
    double k = -1.*(p3.x() - p2.x())/(p3.y() - p2.y());
    double t = y_n - x_n*k;
    //the eq. is: y = x*k + t
    
    //now we can calculate the x0 and y0 of p0
    double x0 = (t - c)/(m - k);//(c - t) * (k*m)/(m-k);
    p0.setX(x0);
    double y0 = m*x0 + c;//(c - t) * m / (m - k) + t;
    p0.setY(y0);
    
    //now calculate the radius,phi0, tanLambda assuming that the helix
    //crosses the point (0,0,z), which is the center of the capture
    //targets
    double deltaX = p3.x() - x0;
    double deltaY = p3.y() - y0;
    double deltaZ = p3.z() - p2.z();
    //    radius    = std::sqrt(x0*x0 + y0*y0);
    radius    = std::sqrt((x0 - p3.x())*(x0 - p3.x()) + (y0 - p3.y())*(y0 - p3.y()));

    double delta2Y = (p2.y() - y0);
    double delta2X = (p2.x() - x0);
    phi0      = std::acos( delta2X/radius);
    //now compute phi0 using y coordinates so to estabilish which sign phi0 must have
    double phi0_fromY = std::asin( delta2Y/radius);
    //    phi0     = phi0*phi0_fromY/std::fabs(phi0_fromY);
    double phi0_fromXY = std::atan2(  delta2Y, delta2X);//std::atan2(  (p2.y() - y0), (p2.x() - x0));
    //2014 - 02 - 01 gianipez resolved a bug?
    //if( phi0*phi0_fromY < 0.)
    //  phi0     = phi0*phi0_fromY/std::fabs(phi0_fromY);
    if (_debug > 5) {
      printf("[HelixFitHack:calculateTrackParameters] phi0 from X = %5.3f phi0 from Y = %5.3f phi0 from tan = %5.3f p1.z = %10.3f p2.z = %10.3f p3.z = %10.3f\n",
	     phi0, phi0_fromY, phi0_fromXY,p1.z(),p2.z(),p3.z());
    }
    phi0 = phi0_fromXY;

    double deltaPhi_0 = std::atan2(deltaY,deltaX) - phi0;
    double deltaPhi_1 = std::atan(deltaY/deltaX) - phi0;
    if (_debug > 5) {
      printf("[HelixFitHack:calculateTrackParameters] deltaPhi_0 = %5.5f deltaPhi_1 = %5.5f\n",
	     deltaPhi_0,
	     deltaPhi_1);
    }
    
    if(deltaPhi_0 < 0.) 
      deltaPhi_0 += 2.*M_PI;

    tanLambda = (radius/deltaZ)*deltaPhi_0;//(radius/deltaZ)*std::fabs(deltaPhi_0);

    if (_debug > 5) {
      printf("[HelixFitHack:calculateTrackParameters] dfdz = %5.8f \n", 
	     tanLambda/radius);
    }
  }
  
//-----------------------------------------------------------------------------
// 2015-03-15 P.Murat: WHY the input parameters are passed by reference ?
//-----------------------------------------------------------------------------
  void  HelixFitHack::calculateDfDz(double &phi0, double &phi1, 
				    double &z0,   double &z1,
				    double &dfdz) {

    double deltaPhi = TVector2::Phi_mpi_pi(phi1-phi0);

    dfdz = deltaPhi/(z1 - z0);
    if (dfdz > 0.0) {
      if (_debug > 5) {
	printf("[HeliFitHack::calculateDfDZ] df = %5.3f dz = %5.3f dfdz = %5.5f\n",
	       deltaPhi, (z1 - z0), dfdz);
      }
    }
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void HelixFitHack::plotXY(int ISet) {
  //  TFolder* fol;

//   fol = (TFolder*) gROOT->GetRootFolder()->FindObject("Mu2e/CalPatRec");

//   if (fol == NULL) {
//     printf(" fol is not defined, exit\n");
//     return;
//   }

//   Ref* ref = (Ref*) fol->FindObject("Ref");
//   mu2e::HelixFitHack* hfit = ref->fHelixFit;

  std::vector<mu2e::XYZPHack>* xyzp;
  HelixFitHackResult*          helx;

  if ((ISet >=0 ) && (ISet < 6)) {
    xyzp   = &_results[ISet]._xyzp;
    helx   = &_results[ISet]._helix;
  }
  else {
    printf("ISet = %i undefined, return\n",ISet);
    return;
  }

  int nhits = xyzp->size();

  printf("nhits = %10i\n",nhits);

  mu2e::XYZPHack*    hit;
  CLHEP::Hep3Vector* pos;

  double x[1000], y[1000];//, z[1000];
  int    flag[1000];

  for (int i=0; i<nhits; i++) {
    hit  = &xyzp->at(i);
    pos  = &hit->_pos;

    x[i] = hit->_pos.x();
    y[i] = hit->_pos.y();
    //    z[i] = hit->_pos.z();

    flag[i] = *((int*) &hit->_flag);

    printf("i: %3i ind: %5i x: %10.3f y: %10.3f z: %10.3f 0x%08x %5i\n",
	   i,(int) hit->_ind,
	   pos->x(),pos->y(),pos->z(),
	   flag[i], 
	   hit->isOutlier());
  }

  char name_xy[200];

  TCanvas* c;
  TMarker* m;
  int      color;

  sprintf(name_xy,"c_plot_hits_xy_%i",ISet);
  c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name_xy);
  if (c) delete c;

  c = new TCanvas(name_xy,"c_xy",1000,1000);
  //  c->Divide(2,1);
//-----------------------------------------------------------------------------
// plot XY view
//-----------------------------------------------------------------------------
  c->cd(1);

  TH2F* h2_xy = new TH2F("h2_xy",Form("XY View %i",ISet), 140,-700,700,140,-700,700);
  h2_xy->SetStats(0);
  h2_xy->Draw();


  for (int i=0; i<nhits; i++) {
    hit = &xyzp->at(i);
    m   = new TMarker(x[i],y[i],2);

    if (hit->isOutlier()) color = kBlack;
    else                  color = kRed;

    m->SetMarkerSize(0.7);
    m->SetMarkerColor(color);
    m->Draw();
  }

  if (helx->fitIsValid()) {
//-----------------------------------------------------------------------------
// draw unweighted  helix
//-----------------------------------------------------------------------------
    double x0, y0, r;

    x0  = helx->_sxy.x0();
    y0  = helx->_sxy.y0();
    r   = helx->_sxy.radius();

    TEllipse* e = new TEllipse(x0,y0,r);
    e->SetFillStyle(0);
    e->SetLineColor(2);
    e->Draw();
//-----------------------------------------------------------------------------
// draw weighted  helix
//-----------------------------------------------------------------------------
    if (helx->weightedFitIsValid()) {
      x0  = helx->_sxyw.x0();
      y0  = helx->_sxyw.y0();
      r   = helx->_sxyw.radius();
    
      e   = new TEllipse(x0,y0,r);
      e->SetFillStyle(0);
      e->SetLineColor(3);
      e->Draw();
    }
  }

//-----------------------------------------------------------------------------
// plot YZ view
//-----------------------------------------------------------------------------
//   c->cd(2);

//   TH2F* h2_yz = new TH2F("h2_yz",Form("YZ VIEW %i",ISet),1600,-1600,1600,140,-700,700);
//   h2_yz->SetStats(0);
//   h2_yz->Draw();

//   for (int i=0; i<nhits; i++) {
//     hit = &xyzp->at(i);
//     m = new TMarker(z[i],y[i],2);

//     if (hit->isOutlier()) color = kBlack;
//     else                  color = kRed;

//     m->SetMarkerColor(color);
//     m->SetMarkerSize(0.7);
//     m->Draw();
//   }

  helx->print("from HelixFitHack::plotXY");

  printf("All Done\n");
}
//-----------------------------------------------------------------------------
// this routine is supposed to be called interactively from the ROOT prompt
// it has to retrieve a pointer to HelixFitHack called from CalPatRec
//-----------------------------------------------------------------------------
  void HelixFitHack::plotZPhi(int ISet) {

    TCanvas* c;
    int      color;
    TMarker* m;
    char     name[200];
//-----------------------------------------------------------------------------
// retrieve the data points for storing in TGraphs
//-----------------------------------------------------------------------------
    std::vector<mu2e::XYZPHack>* xyzp;
    HelixFitHackResult*          helx;
    
    if ((ISet >=0 ) && (ISet < 6)) {
      xyzp   = &_results[ISet]._xyzp;
      helx   = &_results[ISet]._helix;
    }
    else {
      printf("ISet = %i undefined, return\n",ISet);
      return;
    }

    int nhits = xyzp->size();

    printf("nhits = %10i\n",nhits);

    mu2e::XYZPHack* hit;
    
    CLHEP::Hep3Vector* pos;

    double x[1000], y[1000];
    double z[1000], phi[1000], phi1[1000];

    int    flag[1000];
    
    for (int i=0; i<nhits; i++) {
      hit = & xyzp->at(i);
      pos = & hit->_pos;
      
      x[i]    = hit->_pos.x();
      y[i]    = hit->_pos.y();
      z[i]    = hit->_pos.z();
      flag[i] = *((int*) &hit->_flag);
      phi [i] = hit->_pos.phi();

      phi1[i] = atan2(y[i]-helx->_center.y(),x[i]-helx->_center.x());
      
      printf("i: %3i ind: %5i x: %10.3f y: %10.3f z: %10.3f phi: %10.3f phi1: %10.3f 0x%08x %5i\n",
	     i,(int) hit->_ind,
	     pos->x(),pos->y(),pos->z(),phi[i],phi1[i],
	     flag[i], hit->isOutlier());
    }

//-----------------------------------------------------------------------------
// straight line in RPhi out of the fit
//-----------------------------------------------------------------------------
    double phi0 = helx->_fz0;
    double dfdz = helx->_dfdz;

    sprintf(name,"c_plot_hits_phiz_%i",ISet);
  
    c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);

    if (c != NULL) delete c;

    c = new TCanvas(name,Form("c_phiz %i",ISet),1200,1000);
    c->Divide(1,2);
//-----------------------------------------------------------------------------
// plot PHI-Z picture
//-----------------------------------------------------------------------------
    c->cd(1);
    gPad->SetGrid();

    TH2F* h2_phiz = new TH2F("h2_phiz",Form("phiZ VIEW %i",ISet),2200,-1600,2800,80,-20,20);
    h2_phiz->SetStats(0);
    h2_phiz->Draw();

    for (int i=0; i<nhits; i++) {
      hit = &xyzp->at(i);
//-----------------------------------------------------------------------------
// predict phi at this Z to resolve 2pi ambiguities
//-----------------------------------------------------------------------------
      double pred = phi0+dfdz*z[i];
      double phib = phi1[i];

      while (phib-pred >  M_PI) phib -= 2*M_PI;
      while (phib-pred <- M_PI) phib += 2*M_PI;

      m   = new TMarker(z[i],phib,2);

      if (hit->isOutlier()) color = kBlack;
      else                  color = kRed;

      m->SetMarkerColor(color);
      m->SetMarkerSize(0.7);
      m->Draw();
    }

    TF1 *yf = new TF1("yf","[0]+x*[1]",-1600., 2800.);
    yf->SetParameters(phi0, dfdz);
    yf->SetLineWidth(1);
    yf->Draw("same");
//-----------------------------------------------------------------------------
// plot residuals from the straight line fit
//-----------------------------------------------------------------------------
    c->cd(2);
    gPad->SetGrid();

    TH2F* h22_phiz = new TH2F("h22_phiz",Form("Delta(phi) Vs Z, VIEW %i",ISet),2200,-1600,2800,200,-1,1);
    h22_phiz->SetStats(0);
    h22_phiz->Draw();

    double dphi;

    for (int i=0; i<nhits; i++) {
      hit = &xyzp->at(i);
//-----------------------------------------------------------------------------
// predict phi at this Z to resolve 2pi ambiguities
//-----------------------------------------------------------------------------
      double pred = phi0+dfdz*z[i];
      double phib = phi1[i];

      dphi = TVector2::Phi_mpi_pi(phib-pred);
      m    = new TMarker(z[i],dphi,2);

      if (hit->isOutlier()) color = kBlack;
      else                  color = kRed;

      m->SetMarkerSize(0.7);
      m->SetMarkerColor(color);
      m->Draw();
    }

    TF1 *yf2 = new TF1("yf","[0]+x*[1]",-1600., 2800.);
    yf2->SetLineWidth(1);
    yf2->SetParameters(0, 0);
    yf2->Draw("same");

    helx->print("from HelixFitHack::plotZPhi");
  }

//-----------------------------------------------------------------------------
// 'TrkIndex' = 0 : mark all hits as active 
//-----------------------------------------------------------------------------
  void HelixFitHack::saveResults(XYZPHackVector&      Xyzp    , 
				 HelixFitHackResult&  Helix   , 
				 int                  Index   ) {
    _results[Index]._xyzp  = Xyzp;
    _results[Index]._helix = Helix;
  }

}
