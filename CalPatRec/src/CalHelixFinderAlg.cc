///////////////////////////////////////////////////////////////////////////////
// Helix fit to straw hits
//
// $Id: CalHelixFinderAlg.cc,v 1.13 2014/06/06 21:35:08 murat Exp $
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
#include "BTrk/TrkBase/HelixTraj.hh"
#include "CalPatRec/inc/CalHelixFinderAlg.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
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
#include <array>
#include <string>
#include <algorithm>

#include "CalPatRec/inc/CalHelixFinderAlg.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"

using CLHEP::HepVector;
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;

namespace mu2e {

//   // comparison functor for ordering points
//   struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
//     bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
//   };

  // comparison functor for sorting by z
  struct zcomp : public std::binary_function<CalHelixPoint,CalHelixPoint,bool> {
    bool operator()(CalHelixPoint const& p1, CalHelixPoint const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::defineHelixParams(CalHelixFinderData& Helix) const {

    static const double pi(M_PI), halfpi(pi/2.0);

    HepVector pvec(5,0), perr(5,0);

    // the helix fit introduces a radial bias due to an asymmetry in the detector (more phase space for
    // noise hits outside the circle than inside.  correct for it.
    double radius = Helix._radius; //  + _rbias;
    // omega is the inverse transverse radius of the particle's circular motion.
    // It is signed by the particle angular momentum about the cirle center.
    // This CANNOT be deduced geometrically, so must be supplied as an ad-hoc assumption
    double amsign = copysign(1.0,-Helix._tpart.charge()*bz());
    pvec[HelixTraj::omegaIndex] = amsign/radius;
    // phi0 is the azimuthal angle of the particle velocity vector at the point
    // of closest approach to the origin.  It's sign also depends on the angular
    // momentum.  To translate from the center, we need to reverse coordinates
    float  dx     = amsign*Helix._center.y();
    float  dy     = -amsign*Helix._center.x();
    pvec[HelixTraj::phi0Index] = atan2(dy, dx);//-amsign*Helix._center.x(),amsign*Helix._center.y());
    // d0 describes the distance to the origin at closest approach.
    // It is signed by the particle angular momentum WRT the origin.
    // The Helix fit radial bias is anti-correlated with d0; correct for it here.
    pvec[HelixTraj::d0Index] = amsign*(sqrtf(Helix._center.Perp2()) - Helix._radius); //  - 2*_rbias);
    // the dip angle is measured WRT the perpendicular.  It is signed by the particle Z momentum
    pvec[HelixTraj::tanDipIndex] = amsign/(radius*Helix._dfdz);
    // must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
    // resolve the loop ambiguity such that the POCA is closest to z=0.
    double dphi = deltaPhi(Helix._fz0+amsign*halfpi,pvec[HelixTraj::phi0Index]);
    // choose z0 (which loop) so that f=0 is as close to z=0 as possible
    pvec[HelixTraj::z0Index] = dphi*pvec[HelixTraj::tanDipIndex]/pvec[HelixTraj::omegaIndex];
    // estimated covariance based on average performance.  These should be parameters, FIXME!!!
    perr[HelixTraj::d0Index]     = 34.0;
    perr[HelixTraj::phi0Index]   = 0.02;
    perr[HelixTraj::omegaIndex]  = 0.0002;
    perr[HelixTraj::tanDipIndex] = 0.05;
    perr[HelixTraj::z0Index]     = 15.0;

    HepSymMatrix hcov = vT_times_v(perr);

    HelixTraj ht(pvec,hcov);

    Helix._helix = ht.clone();
  }

//-------------------------------------------------------------------------//
//  2014-12-26 Gianipez added the following method for asking t the helix
// finder if the hit "index" has already been marked as used by a previous
// by a previous search
// 2017-12-27: as _markHitCandidates = 0, currently hits are never marked as used
//-----------------------------------------------------------------------------
  int   CalHelixFinderAlg::isHitUsed(int index) {
    return 0;

    // if (index >= kMaxNHits) {
    //   printf("[CalHelixFinderAlg::isHitUsed] requested index = %i range exceeded the range allowed\n", index);
    //   return 1;
    // }

    // int rc = (_indicesTrkCandidate[index] > 0);
    // return rc;
  }

//-----------------------------------------------------------------------------
  CalHelixFinderAlg::CalHelixFinderAlg(fhicl::ParameterSet const& pset) :
    _diag             (pset.get<int>   ("diagLevel"        )),
    _debug            (pset.get<int>   ("debugLevel"       )),
    _debug2           (pset.get<int>   ("debugLevel2"      )),
    _hsel             (pset.get<vector<string> >("HelixFitSelectionBits"  )),
    _bkgsel           (pset.get<vector<string> >("BackgroundSelectionBits")),
    _maxElectronHitEnergy(pset.get<double>("maxElectronHitEnergy")),
    _minNHits         (pset.get<int>   ("minNHit"          )),
    //    _maxDz            (pset.get<double>("maxdz",35.0)),
    _mpDfDz           (pset.get<double>("mostProbableDfDz")),
    // _minNSt           (pset.get<double>("minNActiveStationPairs")),
    _dzOverHelPitchCut(pset.get<double>("dzOverHelPitchCut")),
    _maxDfDz          (pset.get<double>("maxDfDz",0.01)),
    _minDfDz          (pset.get<double>("minDfDz",5e-04)),
    _sigmaPhi         (pset.get<double>("sigmaPhi")),
    _weightXY         (pset.get<double>("weightXY")),
    _weightZPhi       (pset.get<double>("weightZPhi")),
    _weight3D         (pset.get<double>("weight3D")),
    // _ew               (pset.get<double>("errorAlongWire")),
    _maxXDPhi         (pset.get<double>("maxXDPhi",5.)),
    _maxPanelToHelixDPhi(pset.get<double>("maxPanelToHelixDPhi",1.309)),// 75 degrees
    _distPatRec       (pset.get<double>("distPatRec")),
    // _rhomin           (pset.get<double>("rhomin",350.0)),
    // _rhomax           (pset.get<double>("rhomax",780.0)),
    _mindist          (pset.get<double>("mindist",500.0)),
    // _maxdist          (pset.get<double>("maxdist",500.0)),
    _pmin             (pset.get<double>("minP",50.0)),
    _pmax             (pset.get<double>("maxP",150.0)),
    _tdmin            (pset.get<double>("minAbsTanDip",0.3)),
    _tdmax            (pset.get<double>("maxAbsTanDip",2.0)),
    // _rcmin            (pset.get<double>("rcmin",200.0)),
    // _rcmax            (pset.get<double>("rcmax",350.0)),
    _xyweights        (pset.get<bool>  ("xyWeights",false)),
    _zweights         (pset.get<bool>  ("zWeights",false)),
    _filter           (pset.get<bool>  ("filter",true)),
    _plotall          (pset.get<bool>  ("plotall",false)),
    _usetarget        (pset.get<bool>  ("usetarget",true)),
    _bz               (0.0),
    _nHitsMaxPerPanel      (pset.get<int>("nHitsMaxPerPanel"     )),
    _hitChi2Max            (pset.get<double>("hitChi2Max"         )),
    _chi2xyMax             (pset.get<double>("chi2xyMax")),
    _chi2zphiMax           (pset.get<double>("chi2zphiMax")),
    _chi2hel3DMax          (pset.get<double>("chi2hel3DMax")),
    _dfdzErr               (pset.get<double>("dfdzErr")){

    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    CalHelixPoint::_useflag = StrawHitFlag(bitnames);
  }


//-----------------------------------------------------------------------------
  CalHelixFinderAlg::~CalHelixFinderAlg() {
  }

//-----------------------------------------------------------------------------
  double CalHelixFinderAlg::bz() const {
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
//2016-12-26 gianipez added the following function to make CalHelixFinderAlg compatible with TimeCluster obj
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::findHelix(CalHelixFinderData& Helix, const TimeCluster* TimePeak) {

    fTimeCluster = TimePeak;
    //check presence of a cluster
    const CaloCluster* cl = TimePeak->caloCluster().get();
    if (cl == NULL)   return false;

    //fill the calorimeter cluster info
    Hep3Vector         gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskId(),cl->cog3Vector());
    Hep3Vector         tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
    fCaloTime = cl->time();
    fCaloX    = tpos.x();
    fCaloY    = tpos.y();
    fCaloZ    = tpos.z();
//-----------------------------------------------------------------------------
//  compute the allowed radial range for this fit
//-----------------------------------------------------------------------------
    double pb = fabs((CLHEP::c_light*1e-3)/(bz()*Helix._tpart.charge()));
    _rmin = _pmin/(pb*sqrt(1.0+_tdmax*_tdmax));
    _rmax = _pmax/(pb*sqrt(1.0+_tdmin*_tdmin));
//-----------------------------------------------------------------------------
//  particle charge, field, and direction affect the pitch range
//-----------------------------------------------------------------------------
    _dfdzsign = copysign(1.0,-Helix._tpart.charge()*Helix._fdir.dzdt()*bz());

    if(_dfdzsign > 0.0){
      _smin = 1.0/(_rmax*_tdmax);
      _smax = 1.0/(_rmin*_tdmin);
    } else {
      _smax = -1.0/(_rmax*_tdmax);
      _smin = -1.0/(_rmin*_tdmin);
    }
//-----------------------------------------------------------------------------
// call down
//-----------------------------------------------------------------------------
    bool retval = findHelix(Helix);

    return retval;
  }

//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::findHelix(CalHelixFinderData& Helix, const CalTimePeak* TimePeak) {

    fTimePeak = TimePeak;
					// fill the calorimeter cluster info
    fCaloTime = TimePeak->ClusterT0();
    fCaloX    = TimePeak->ClusterX();
    fCaloY    = TimePeak->ClusterY();
    fCaloZ    = TimePeak->ClusterZ();
//-----------------------------------------------------------------------------
//  compute the allowed radial range for this fit
//-----------------------------------------------------------------------------
    double pb = fabs((CLHEP::c_light*1e-3)/(bz()*Helix._tpart.charge()));
    _rmin = _pmin/(pb*sqrt(1.0+_tdmax*_tdmax));
    _rmax = _pmax/(pb*sqrt(1.0+_tdmin*_tdmin));
//-----------------------------------------------------------------------------
//  particle charge, field, and direction affect the pitch range
//-----------------------------------------------------------------------------
    _dfdzsign = copysign(1.0,-Helix._tpart.charge()*Helix._fdir.dzdt()*bz());

    if(_dfdzsign > 0.0){
      _smin = 1.0/(_rmax*_tdmax);
      _smax = 1.0/(_rmin*_tdmin);
    } else {
      _smax = -1.0/(_rmax*_tdmax);
      _smin = -1.0/(_rmin*_tdmin);
    }
//-----------------------------------------------------------------------------
// call down
//-----------------------------------------------------------------------------
    bool retval = findHelix(Helix);

    return retval;
  }

//-----------------------------------------------------------------------------
// called internally; in the diagnostics mode save several states of _xyzp
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::findHelix(CalHelixFinderData& Helix) {
    bool retval(false);
					// initialize internal array of hits, print if requested
    fillXYZP(Helix);
//-----------------------------------------------------------------------------
// 2014-11-09 gianipez: reset the track candidate parameters if a new time peak is used!
// so the previous candidate should not be compared to the new one at this level
//-----------------------------------------------------------------------------
    resetTrackParamters();
//-----------------------------------------------------------------------------
// save results in the very beginning
//-----------------------------------------------------------------------------
    if (_diag > 0) saveResults(Helix, 0);

    // if (_filter) {
    //   filterDist(Helix);
    //   if (_diag > 0) saveResults(_xyzp,Helix,1);
    // }

    doPatternRecognition(Helix);
//---------------------------------------------------------------------------
// 2014-11-11 gianipez changed the following if() statement to test the
// possibility of spead up the pattern recognition in presence of background
// how the number of good points may be different from the number used in sums ?
//---------------------------------------------------------------------------
    if (_debug != 0) {
      printf("[CalHelixFinderAlg::findHelix] Helix._nXYSh = %i goodPointsTrkCandidate = %i\n",
	     Helix._nXYSh, Helix._nStrawHits);//_goodPointsTrkCandidate);
    }

    if (Helix._nStrawHits < _minNHits ) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,1); // small number of hits
    }
    else if ((Helix._radius < _rmin) || (Helix._radius > _rmax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,2); // initialization failure
    }
    else if ((Helix._nXYSh < _minNHits) || (Helix._sxy.chi2DofCircle() > _chi2xyMax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,3); // xy reconstruction failure
    }
    else if ((Helix._nZPhiSh < _minNHits) || (Helix._szphi.chi2DofLine() > _chi2zphiMax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,4); // phi-z reconstruction failure
    }
    else {
//-----------------------------------------------------------------------------
// success, form output
//-----------------------------------------------------------------------------
      Helix._goodhits.clear();

      PanelZ_t*      panelz(0);
      CalHelixPoint* hit(0);

      for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
	panelz = &Helix._oTracker[p];
	int  nhits          = panelz->fNHits;
	for (int i=0; i<nhits; ++i){   
	  hit = &panelz->fHitData.at(i);
	  int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	  if (Helix._hitsUsed[index] != 1)                continue;
	  Helix._goodhits.push_back(hit->index());
	}
      }

      defineHelixParams(Helix);
      retval = true;
    }

    return retval;
  }

//-----------------------------------------------------------------------------
  int CalHelixFinderAlg::findDfDz(CalHelixFinderData& Helix, 
				  SeedInfo_t          SeedIndex,
 				  // int *IndexVec, 
				  int                 Diag_flag) {
    //    return findDfDz_1(Helix, SeedIndex, IndexVec, Diag_flag);
    return findDfDz_2(Helix, SeedIndex, Diag_flag);
  }

//----------------------------------------------------------------------------------------
// 2015-01-13  calculate track DphiDz using histogrammed distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  int CalHelixFinderAlg::findDfDz_1(CalHelixFinderData& Helix, SeedInfo_t SeedIndex, int  Diag_flag) {

//     double phi, phi_ref(-1e10), z, z_ref, dphi, dz, dzOverHelPitch;

//     CLHEP::Hep3Vector* center = &Helix._center;
//     CLHEP::Hep3Vector pos_ref;

//     _hDfDzRes->Reset();
//     _hPhi0Res->Reset();
// 					// 2015 - 03 -30 G. Pezzu changed the value of tollMax.
// 					// using the initial value of dfdz we can set it more accuratelly:
// 					// tollMax = half-helix-step = Pi / dfdz
//     double tollMin(100.);
// //-----------------------------------------------------------------------------
// // 2017-09-26 gianipez fixed a bug: in case the Helix phi-z fit didn't converge yet, 
// // Helix._dfdz is set to -1e6, so we need to make a check here!
// // this is a tempOrary fix that doesn't take into account the particle helicity. FIX ME!
// //-----------------------------------------------------------------------------
//     double helix_dfdz(_mpDfDz);
//     // 2017-11-14 gianipez: findDfDz shoudl use the dfdz value obtained only from the linearFit
//     if (Helix._szphi.qn() >= 10) helix_dfdz = Helix._szphi.dfdz();
//     //    if (Helix._dfdz > 0) helix_dfdz =  Helix._dfdz;
//     double tollMax = 2.*M_PI / helix_dfdz; 

//     if (_debug > 5) {
//       printf("[CalHelixFinderAlg::findDfDz:BEGIN] x0 = %9.3f y0 = %9.3f Helix._radius = %9.3f",
// 	     center->x(), center->y(), Helix._radius);
//       printf("helix_dfdz = %9.6f Helix._nStrawHits = %3i tollMax = %8.6f\n",helix_dfdz, Helix._nStrawHits, tollMax);
//     }

//     int       nstations, nhits[30];
//     double    phiVec[30], zVec[30], weight(0), weight_cl(0);
//     PanelZ_t* panelz(0);

//     // np        = _xyzp.size();
//     nstations = _tracker->nStations();

//     for (int i=0; i<nstations; i++) {
//       phiVec[i] = 0;
//       zVec  [i] = 0;
//       nhits [i] = 0;
//     }
//     //-----------------------------------------------------------------------------
//     // Part 1: use only contiguous parts of the trajectory
//     //-----------------------------------------------------------------------------
//     for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
//       panelz = &Helix._oTracker[p];
//       int  nhitsPerPanel  = panelz->fNHits;
//       int  seedPanelIndex(0);
//       if (nhitsPerPanel == 0)                                                            continue;
//       if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

//       for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
// 	CalHelixPoint* hit = &panelz->fHitData.at(i);
// 	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
// 	if (Helix._hitsUsed[index] != 1)                     continue;

// 	int ist = hit->_straw->id().getStation();                   // station number
// 	phi     = polyAtan2(hit->y()-center->y(),hit->x()-center->x()); // atan2 returns its result in [-pi,pi], convert to [0,2pi]
// 	if (phi < 0) phi += 2*M_PI;
// 	zVec  [ist] += hit->z();
// 	//-----------------------------------------------------------------------------
// 	// make sure there all hits within the station get close values of phi, although a 
// 	// common 2pi ambiguity is still unresolved
// 	//-----------------------------------------------------------------------------
// 	if (nhits[ist] == 0) phiVec[ist] = phi;
// 	else {
// 	  while (phi-phiVec[ist] >  M_PI) phi -= 2*M_PI;
// 	  while (phi-phiVec[ist] < -M_PI) phi += 2*M_PI;
	
// 	  phiVec[ist] = (phiVec[ist]*nhits[ist]+phi)/(nhits[ist]+1);
// 	}
// 	nhits [ist] += 1;
//       }
//     }

//     for (int i=0; i<nstations; i++) {
//       if (nhits[i] > 0) {
// 	zVec  [i] = zVec  [i]/nhits[i];
//       }
//     }
    
//     if (_debug >5) {
//       printf("[CalHelixFinderAlg::findDfDz] StationID  nhits       z        phi\n");
//       for (int i=0; i<nstations; i++) {
// 	if (nhits[i] > 0) printf("[CalHelixFinderAlg::findDfDz] %5i %6i    %9.3f %8.5f\n", i,nhits[i],zVec[i],phiVec[i]);
//       }
//     }

//     int i0(-1), first_point(1);
// 					// add the cluster phi
//     double zCl   = fCaloZ;
//     double phiCl = polyAtan2(fCaloY-center->y(),fCaloX-center->x());
//     if (phiCl < 0) phiCl += 2*M_PI;

//     for (int i=0; i<nstations; i++) {
//       if (nhits[i] == 0)                                    continue; 
// 				        // i0: fist station with hits
//       if (first_point) {
// 	i0          = i;
// 	first_point = 0;
//       }

//       phi_ref = phiVec[i];
//       z_ref   = zVec  [i];

//       for(int j=i+1; j<nstations; ++j) {
// 	if (nhits[j] == 0)                                  continue;
// 	phi = phiVec[j];
// 	z   = zVec  [j];
// 	dz  = z - z_ref;
	
// 	dzOverHelPitch = dz/tollMax - int(dz/tollMax);
// 	weight         = nhits[i] + nhits[j];

// 	//	if ((phi_ref > -9999) && (dzOverHelPitch < _dzOverHelPitchCut) && (dz > tollMin)) {
// 	if (dz > tollMin) {
// 	  dphi = phi-phi_ref;
//  	  while (dphi >  M_PI) dphi -= 2*M_PI;
//  	  while (dphi < -M_PI) dphi += 2*M_PI;
// //-----------------------------------------------------------------------------
// // add 2*PI to take into account the fact we are in the second loop
// // FIX ME: what to do if we are in the third loop?
// //-----------------------------------------------------------------------------
// 	  if (dz > tollMax) dphi += 2*M_PI*int(dz/tollMax);

// 	  double dphidz = dphi/dz;
// 	  while (dphidz < 0.) {
// 	    dphi    = dphi+2.*M_PI;
// 	    dphidz  = dphi/dz;
// 	  }
// 	  _hDfDzRes->Fill(dphidz, weight);

// 	  double tmpphi0 = phi_ref - dphidz*z_ref;
// 	  tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

// 	  if (_debug > 5) {
// 	    printf("[CalHelixFinderAlg::findDfDz:1] z_ref: %9.3f z: %9.3f dz: %9.3f",z_ref,z,dz);
// 	    printf(" phi_ref: %9.5f phi: %9.5f dphi: %9.5f dz/HelPitch: %10.3f dphi/dz: %9.5f phi0 = %9.6f\n",
// 		   phi_ref,phi,dphi, dzOverHelPitch, dphidz, tmpphi0);
// 	  }
// //-----------------------------------------------------------------------------
// // in case dfdz is out of limits, set tmpphi0 as negative
// //-----------------------------------------------------------------------------
// 	  if ((dphidz < _minDfDz) || (dphidz >  _maxDfDz)) tmpphi0 = -1;
// 	  _hPhi0Res->Fill(tmpphi0, weight);
// 	}
//       }
// //-----------------------------------------------------------------------------
// // include the calorimeter cluster phi
// //-----------------------------------------------------------------------------
//       dz             = zCl - z_ref;
//       dzOverHelPitch = dz/tollMax - int(dz/tollMax);
//       weight_cl      =  nhits[i];

//       //      if ((dzOverHelPitch < _dzOverHelPitchCut) && (dz > tollMin)) {
//       if (dz > tollMin) {
// 	dphi  = phiCl - phi_ref;
// 	dphi  = TVector2::Phi_0_2pi(dphi);
// //-----------------------------------------------------------------------------
// // add 2 pi for taking into account the fact we are in the second loop
// // *FIX ME*: what if we are in the third loop?
// //-----------------------------------------------------------------------------
// 	if (dz > tollMax) dphi += 2*M_PI*int(dz/tollMax);

// 	double dphidz = dphi/dz;
// 	while (dphidz < 0.) {
// 	  dphi   += 2.*M_PI;
// 	  dphidz = dphi/dz;
// 	}

// 	double tmpphi0 = phi_ref - dphidz*z_ref;
// 	tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

// 	if (_debug > 5){
// 	  printf("[CalHelixFinderAlg::findDfDz:2] z_ref: %9.3f z: %9.3f dz: %9.3f",z_ref,zCl,dz);
// 	  printf(" phi_ref: %9.5f phi: %9.5f dphi: %9.5f dz/HelPitch: %10.3f dphi/dz: %9.5f phi0 = %9.6f\n",
// 		 phi_ref,phiCl,dphi, dzOverHelPitch, dphidz, tmpphi0);
// 	}

// 	if (dzOverHelPitch < _dzOverHelPitchCut ) {
// 	  _hDfDzRes->Fill(dphidz, weight_cl);
// 	  if ((dphidz < _minDfDz) || (dphidz >  _maxDfDz)) tmpphi0 = -1;
// 	  _hPhi0Res->Fill(tmpphi0, weight_cl);
// 	}
//       }
//     }
// //-----------------------------------------------------------------------------
// // 2015 - 04- 02 G. Pezzu changed the way the maximum is searched
// // since sometimes a 2pi ambiguity creates two peaks in the histogram
// // we want to use the second, because it is the correct one
// //-----------------------------------------------------------------------------
//     double  maxContent = _hDfDzRes->GetMaximum() - 0.001;
//     int      maxBin    = _hDfDzRes->FindLastBinAbove(maxContent);//GetMaximumBin();
//     _hdfdz             = _hDfDzRes->GetBinCenter(maxBin);//_hDfDzRes->GetMean();
//     double dfdzmean    = _hDfDzRes->GetMean();
//     int    nentries    = _hDfDzRes->GetEntries();
//     int    overflows   = _hDfDzRes->GetBinContent(0)  + _hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1);

//     maxContent         = _hPhi0Res->GetMaximum() - 0.001;
//     maxBin             = _hPhi0Res->FindLastBinAbove(maxContent);//GetMaximumBin();

//     double mpvphi0     = _hPhi0Res->GetBinCenter(maxBin); //_hPhi0Res->GetMean();
//     double menaphi0    = _hPhi0Res->GetMean();
//     int    nentriesphi = _hPhi0Res->GetEntries();

//     _hphi0 = mpvphi0;  // 2018-01-05: *DOUBLE_CHECK*

//     if (_debug > 5) {
//       printf("[CalHelixFinderAlg::findDfDz:DFDZ] nent: %3i mpvDfDz: %9.6f meanDphiDz: %9.6f under: %3.0f over: %3.0f ENTRIES:",
// 	     nentries, _hdfdz, dfdzmean,
// 	     _hDfDzRes->GetBinContent(0),_hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1)
// 	     );
//       for (int i=0; i<_hDfDzRes->GetNbinsX(); i++) {
// 	printf(" %3.0f",_hDfDzRes->GetBinContent(i+1));
//       }
//       printf("\n");

//       printf("[CalHelixFinderAlg::findDfDz:PHI0] nent: %3i mpvPhi0: %9.6f meanPhi0  : %9.6f under: %3.0f over: %3.0f ENTRIES:",
// 	     nentriesphi, mpvphi0,  menaphi0,
// 	     _hPhi0Res->GetBinContent(0),_hPhi0Res->GetBinContent(_hPhi0Res->GetNbinsX()+1)
// 	     );
//       for (int i=0; i<_hPhi0Res->GetNbinsX(); i++) {
// 	printf(" %3.0f",_hPhi0Res->GetBinContent(i+1));
//       }
//       printf("\n");
//     }
// //-----------------------------------------------------------------------------
// // Part 2: try to perform a more accurate estimate - straight line fit
// //-----------------------------------------------------------------------------
//     double z0, phi0, dphidz, pred;

//     z0     = 0.    ;
//     phi0   = _hphi0;
//     dphidz = _hdfdz;
//     //    _sdfdz = -1;

//     if (_debug > 5) {
//       double tmpphi0=phi0+dphidz*z0;
//       printf("[CalHelixFinderAlg::findDfDz:PART2] phi0 = %9.6f dfdz = %9.6f\n", tmpphi0, dphidz);
//     }
// //--------------------------------------------------------------------------------
// // 2015-03-25 G. Pezzu changed the way the 2PI ambiguity is resolved
// //--------------------------------------------------------------------------------
//     LsqSums4 szphi;

//     weight = 1./(_sigmaPhi*_sigmaPhi);

//     double xdphi, zLast(z0), zdist;

//     if (_debug > 5) {
//       printf("[CalHelixFinderAlg::findDfDz:LOOP]  i       z        dz      phiVec     phi    szphi.dfdz    dphi    xdphi     dfdz    \n");
//     }

//     dz   = zCl - z0;
//     dphi = dz*dphidz + phi0 - phiCl;
//     while (dphi > M_PI){
//       phiCl += 2*M_PI;
//       dphi  -= 2*M_PI; // dz*dphidz + phi0 - phiCl;
//     }
//     while (dphi < -M_PI){
//       phiCl -= 2*M_PI;
//       dphi  += 2*M_PI; // = dz*dphidz + phi0 - phiCl;
//     }

//     double errCl    = 2./30.;
//     double weightCl = 1./(errCl*errCl);

//     xdphi = fabs(dphi)/errCl;
// //-----------------------------------------------------------------------------
// // 2015-04-21 Gianipez added the condition (xdphi < 2.*_maxXDPhi) for adding or not
// // the calorimeter point to the fitter. In case the particle scattered in the end of the tracker
// // the calorimeter point is dangerous.
// //-----------------------------------------------------------------------------
//     if (xdphi < 2.*_maxXDPhi){
//       szphi.addPoint(zCl, phiCl, weightCl);
//     }

//     if (_debug > 10) {
//       printf("[CalHelixFinderAlg::findDfDz:LOOP] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
// 	     0, zCl, dz, phiCl, phiCl, szphi.dfdz(), dphi, xdphi, dphidz);
//     }

//     // 2015-07-06 Gianipez added the following line for avoiding infinite loops
//     if ( i0 < 0 ) goto NEXT_STEP;

//     for (int i=i0; i<nstations; i++) {
//       if (nhits[i] > 0) {
// 	z    = zVec[i];
// 	dz   = z-z0;
// 	pred = phi0 + dz*dphidz;
// 	phi  = phiVec[i];
// 	dphi = phi - pred;//pred - phi;

// 	while (dphi > M_PI){
// 	  phi -= 2*M_PI;//+= 2*M_PI;
// 	  dphi = phi - pred;//pred - phi;
// 	}
// 	while (dphi < -M_PI){
// 	  phi += 2*M_PI;//-= 2*M_PI;
// 	  dphi = phi - pred;//pred - phi;
// 	}

// 	xdphi = fabs(dphi)/_sigmaPhi;

// 	if (xdphi < 2.*_maxXDPhi){
// 	  szphi.addPoint(z, phi, weight);

// 	  zdist = z - zLast;

// 	  if ( (szphi.qn() >= 3.) && (zdist > 500.)){
// 	    z0     = 0.;
// 	    phi0   = szphi.phi0();
// 	    dphidz = szphi.dfdz();
// 	  }
// 	}

// 	if (_debug > 10) {
// 	  double tmpDfDz = szphi.dfdz();//, Helix._szphi.chi2DofLine());
// 	  printf("[CalHelixFinderAlg::findDfDz:LOOP] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
// 		 i, z, dz, phiVec[i], phi, tmpDfDz, dphi, xdphi, dphidz);
// 	}
//       }
//     }

//   NEXT_STEP:;

//     if (szphi.qn() >= 3.) {
//       _hdfdz = szphi.dfdz();		// sigxy/sigxx;
//       _hphi0 = szphi.phi0();		// ymean - xmean*sigxy/sigxx;
//       //      _sdfdz = szphi.chi2DofLine();
//     }
//     else {
//       _hphi0 = phi0 + _hdfdz*z0;
//       //      _sdfdz = -1;
//     }

//     int     nActive_stations = nentries - overflows;
    
//     if (Diag_flag > 0){
//       Helix._diag.nStationPairs = nActive_stations;
//     }

//     if (_debug > 5) {
//       printf("[CalHelixFinderAlg::findDfDz] END: _hdfdz = %9.5f _hphi0 = %9.6f chi2 = %9.3f ", _hdfdz,
// 	     _hphi0, -1.);
//       printf(" FIT: szphi.dfdz() = %9.5f szphi.phi0() = %9.6f chi2 = %9.3f qn = %6.0f\n", szphi.dfdz(),
// 	     szphi.phi0(), szphi.chi2DofLine(), szphi.qn());
//     }
// //----------------------------------------------------------------------------- 
// // 2017-11-14 gianipez: in case of less than _minNSt active stations we should 
// //                      probably use the _mpDfDz
// //-----------------------------------------------------------------------------
//     if (nActive_stations < _minNSt) {
//       _hdfdz = _mpDfDz;
//       return 0;
//     }

    return 1;
  }



//----------------------------------------------------------------------------------------
// 2015-01-13  calculate track DphiDz using histogrammed distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  int CalHelixFinderAlg::findDfDz_2(CalHelixFinderData& Helix, SeedInfo_t SeedIndex, int  Diag_flag) {

    double phi, phi_ref(-1e10), z_ref, dphi, dz;

    double hist[20], minX(0), maxX(0.01), stepX(0.0005), nbinsX(20); // make it 20 bins

    XYZVec* center = &Helix._center;
    XYZVec  pos_ref;
//-----------------------------------------------------------------------------
// 2017-09-26 gianipez fixed a bug: in case the Helix phi-z fit didn't converge yet, 
// Helix._dfdz is set to -1e6, so we need to make a check here!
// this is a tempOrary fix that doesn't take into account the particle helicity. FIX ME!
//-----------------------------------------------------------------------------
    if (_debug > 5) {
      printf("[CalHelixFinderAlg::findDfDz:BEGIN] x0 = %9.3f y0 = %9.3f Helix._radius = %9.3f Helix._nStrawHits = %3i",
	     center->x(), center->y(), Helix._radius,Helix._nStrawHits);
    }

    int       nstations, nhits[30], nstations_with_hits(0);
    double    phiVec[30], zVec[30], weight(0);

    // np        = _xyzp.size();
    nstations = _tracker->nStations();

    for (int i=0; i<nstations; i++) {
      phiVec[i] = 0;
      zVec  [i] = 0;
      nhits [i] = 0;
    }

    for (int i=0; i<nbinsX; i++) hist[i] = 0;
//-----------------------------------------------------------------------------
// calorimeter cluster - point number nstations+1
//-----------------------------------------------------------------------------
    double zCl   = fCaloZ;
    double phiCl = polyAtan2(fCaloY-center->y(),fCaloX-center->x());
    if (phiCl < 0) phiCl += 2*M_PI;

    phiVec[nstations] = phiCl;
    zVec  [nstations] = zCl;
    nhits [nstations] = 1;
    //-----------------------------------------------------------------------------
    // Step 1: for each station with track candidate hits, calculate average phi per station
    //-----------------------------------------------------------------------------
    PanelZ_t* panelz(0);

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhitsPerPanel  = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhitsPerPanel == 0)                                      continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

      for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	CalHelixPoint* hit = &panelz->fHitData.at(i);

	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1 )                         continue;

	int ist = hit->strawId().station();//_straw->id().getStation();                   // station number
	phi     = polyAtan2(hit->y()-center->y(),hit->x()-center->x()); // atan2 returns its result in [-pi,pi], convert to [0,2pi]
	if (phi < 0) phi += 2*M_PI;
	zVec  [ist] += hit->z();
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
    }

    for (int i=0; i<nstations; i++) {
      if (nhits[i] > 0) {
	zVec  [i] = zVec  [i]/nhits[i];
      }
    }

    if (_debug >5) {
      printf("[CalHelixFinderAlg::findDfDz] StationID  nhits       z        phi\n");
      for (int i=0; i<nstations+1; i++) {
	if (nhits[i] > 0) printf("[CalHelixFinderAlg::findDfDz] %5i %6i    %9.3f %8.5f\n", i,nhits[i],zVec[i],phiVec[i]);
      }
    }
//-----------------------------------------------------------------------------
// Step 2: determine the most likely value of phi
//-----------------------------------------------------------------------------
    for (int i=0; i<nstations; i++) {
      if (nhits[i] == 0)                                    continue; 

      phi_ref = phiVec[i];
      z_ref   = zVec  [i];

      for(int j=i+1; j<nstations+1; ++j) { // nstations+1 accounts for the cluster
	if (nhits[j] == 0)                                  continue;

	dphi = phiVec[j]-phi_ref;
	dz   = zVec[j] - z_ref;
	double dphidz = dphi/dz;
	
	weight = nhits[i] + nhits[j];
//-----------------------------------------------------------------------------
// calculate N potential choices for 2*PI ambiguity resolution
//-----------------------------------------------------------------------------
	int n(0), nmax(0), nmin(0), nchoices = 0;

	double x = dphidz + n*2*M_PI/dz;

	if (x < _minDfDz) {
//-----------------------------------------------------------------------------
// for n=0, x < _minDfDz
//-----------------------------------------------------------------------------
	  while (x < _maxDfDz) {
	    if (x < _minDfDz) nmin = n+1;
	    nmax = n;
	    if ((x > _minDfDz) && (x < _maxDfDz)) nchoices += 1;
	    n += 1;
	    x += 2*M_PI/dz;
	  }
	}
	else if (x < _maxDfDz) {
//-----------------------------------------------------------------------------
// for n=0,   _xMin <= x < _xMax
//-----------------------------------------------------------------------------
	  while (x < _maxDfDz) {
	    nmax = n;
	    if ((x > _minDfDz) && (x < _maxDfDz)) nchoices += 1;
	    n += 1;
	    x += 2*M_PI/dz;
	  }

	  nmin = 0;
	  x    = dphidz+(nmin-1)*2*M_PI/dz;

	  while (x > _minDfDz) {
	    nchoices += 1;
	    nmin -= 1;
	    x    -= 2*M_PI/dz;
	  }
	}
	else {
//-----------------------------------------------------------------------------
// for n=0, x >= _xMax
//-----------------------------------------------------------------------------
	  while (x > _minDfDz) {
	    if (x > _maxDfDz) nmax = n-1;
	    nmin = n;
	    if ((x > _minDfDz) && (x < _maxDfDz)) nchoices += 1;
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
    if (nstations_with_hits < 2) _hdfdz = _mpDfDz;
    else                         _hdfdz = xmp;
//-----------------------------------------------------------------------------
// last step - determine phi0 = phi(z=0)
//-----------------------------------------------------------------------------
    double phi0(0), sdphi(0);
    int    sn(0);

    for (int i=0; i<nstations+1; i++) {
      if (nhits[i] == 0) continue;

      if (sn == 0) { // first station with hits gives the "2*PI normalization";
	phi0 = phiVec[i]-zVec[i]*_hdfdz;
	sdphi = 0;
	sn    = 1;
      }
      else {
//-----------------------------------------------------------------------------
// for all points different from the first one need to choose the turn number
//-----------------------------------------------------------------------------
	dphi = phiVec[i]-(phi0+zVec[i]*_hdfdz);
	double dphi_min = dphi;

	int n= 0;
	while (1) {
	  n += 1;
	  double dphi = phiVec[i]+2*M_PI*n-(phi0+zVec[i]*_hdfdz);
	  if (fabs(dphi) < fabs(dphi_min)) dphi_min = dphi;
	  else break;
	}

	n=0;
	while (1) {
	  n -= 1;
	  double dphi = phiVec[i]+2*M_PI*n-(phi0+zVec[i]*_hdfdz);
	  if (fabs(dphi) < fabs(dphi_min)) dphi_min = dphi;
	  else break;
	}
	
	sdphi += dphi_min;
	sn += 1;
      }
    }

    _hphi0 = phi0 + sdphi/sn;

    if (Diag_flag > 0){
      Helix._diag.nStationPairs = nstations_with_hits;
    }

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::findDfDz] END: _hdfdz = %9.5f _hphi0 = %9.6f ", _hdfdz, _hphi0);
    }

    return 1;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::doLinearFitPhiZ(CalHelixFinderData& Helix    ,
					  SeedInfo_t          SeedIndex,
					  int                 UseInteligentWeight,
					  int                 DoCleanUp           ) {

    auto              hitsUsed  = Helix._hitsUsed;

    bool              success(false);
    int               nPointsRemoved(0);
 
    double            phi(0.0), z(0.0), weight(0.0);
    double            dphi, err, xdphi;

    XYZVec            helCenter = Helix._center;
    double            radius    = Helix._radius;
    XYZVec            pos;
    PanelZ_t*         panelz(0);
    CalHelixPoint*    hit(0);
//-----------------------------------------------------------------------------
// gianipez: procedure for aligning the phi vector
//-----------------------------------------------------------------------------
    ::LsqSums4       szphi;
    int              count(0), indexWorst;
    double           chi2, chi2min, deltaPhi, dphi_max(0), phi_ref, weightWorst(0);

    SeedInfo_t       iworst(-1,-1);

    if (Helix._sxy.qn() > 0) {
      helCenter = XYZVec( Helix._sxy.x0(), Helix._sxy.y0(), 0);
      radius    = Helix._sxy.radius();
    }
    XYZVec           strawDir;
    const char       banner[200] = "doLinearFitPhiZ";
//--------------------------------------------------------------------------------
// set EMC cluster info and initilize the dfdz for the search
//-----------------------------------------------------------------------------
    double dfdz  = Helix._dfdz;
    double phi0  = Helix._fz0;

    double zCl   = fCaloZ;
    pos          = XYZVec(fCaloX, fCaloY, fCaloZ);
    double dx    = (pos.x() - helCenter.x());
    double dy    = (pos.y() - helCenter.y());
    double phiCl = polyAtan2(dy, dx);//CLHEP::Hep3Vector(pos - helCenter).phi();//center).phi();
    if (phiCl < 0) phiCl = phiCl + 2*M_PI;//        = TVector2::Phi_0_2pi(phiCl);

    deltaPhi = zCl*dfdz + phi0 - phiCl;
    while (deltaPhi > M_PI){
      phiCl   += 2*M_PI;
      deltaPhi = zCl*dfdz + phi0 - phiCl;
    }
    while (deltaPhi < -M_PI){
      phiCl   -= 2*M_PI;
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

    Helix._szphi.clear();
    Helix._nZPhiSh = 0;

    if ( xdphi < 2.*_maxXDPhi ) {
      double weight_cl = 784.;//10.0; 
      Helix._szphi.addPoint(zCl,phiCl,weight_cl);
      Helix._nZPhiSh += 1;
    }

    count = 0;
    double zlast, dz, dx2, dy2;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ:BEGIN] phi0 = %10.6f dfdz = %10.6f chi2N = %10.3f DoCleanup = %i\n",
	     Helix._fz0,  Helix._dfdz, 0.,DoCleanUp);
      printf("[CalHelixFinderAlg::doLinearFitPhiZ]    flag   A   shID   i       z         ");
      printf("    phi         dphi      xdphi      zlast        dz      dphidz  szphidfdz  chi2\n");
      printf("[CalHelixFinderAlg::doLinearFitPhiZ] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
	     0, 1, 0, 0,  zCl, phiCl, deltaPhi, xdphi, 0., 0., dfdz, 0., 0.);
    }

    // int            nLayers(2);
    float          panelHitChi2 [_nHitsMaxPerPanel] = {1e10};
    int            panelHitIndex[_nHitsMaxPerPanel] = {-1};

    zlast = 0;

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      z      = panelz->z;
      int       nhits  = panelz->fNHits;
      int       seedPanelIndex(0);
      if (nhits == 0)                                                                    continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

      // reset the chi2 values calulated in the previous panel
      // for (int l=0; l<nLayers; ++l){
      for (int g=0; g<_nHitsMaxPerPanel; ++g){
	panelHitChi2 [g] = 1e10;
	panelHitIndex[g] = -1;
      }
      // }

      for (int i=seedPanelIndex; i<nhits; ++i){
	CalHelixPoint* hit = &panelz->fHitData.at(i);
	pos      = hit->_pos;
	strawDir = hit->_sdir;

	dx2      = (pos.x() - helCenter.x());
	dy2      = (pos.y() - helCenter.y());
	phi      = polyAtan2(dy2, dx2);//CLHEP::Hep3Vector(pos - helCenter).phi();//center).phi();
	if (phi < 0) phi = phi + 2*M_PI;//      = TVector2::Phi_0_2pi(phi);
	dz       = z - zlast;
                                    
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
	// phi_corrected[i] = phi;
	// _phiCorrected[i] = phi;
	hit->_hphi = phi;
	dphi      = fabs(dphi);
	err       = _sigmaPhi;

	if (UseInteligentWeight == 1){
	  weight           = calculatePhiWeight(*hit,/*pos, strawDir, */ helCenter, radius, 0, banner);
	  err              = 1./sqrt(weight);
	}
	hit->_zphiWeight = weight;

	xdphi = dphi/err;

	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (_debug > 10) {
	  printf("[CalHelixFinderAlg::doLinearFitPhiZ:LOOP] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
		 *((int*) &hit->_flag), Helix._hitsUsed[index],
		 hit->strawId().straw()/*_strawhit->strawIndex().asInt()*/, i,
		 z, phi, dphi,xdphi,zlast,dz,
		 dfdz, Helix._szphi.dfdz(), Helix._szphi.chi2DofLine());
	}
	if (Helix._hitsUsed[index] != 1)                     continue;

	if ( (DoCleanUp == 1) && (xdphi > _maxXDPhi) ) {
	  Helix._hitsUsed[index] = 0;
	  ++nPointsRemoved;
	  continue;
	}

	double delta_min(0);
	int    index_min(-1);
	// int    layer_id(hit->strawId().layer());//_straw->id().getLayer());
	for (int k=0; k<_nHitsMaxPerPanel; ++k){
	  double delta = xdphi - panelHitChi2[k];
	  if (delta < delta_min){
	    delta_min = delta;
	    index_min = k;
	  }
	}
	if (index_min>=0) {
	  panelHitIndex[index_min] = i;
	  panelHitChi2 [index_min] = xdphi;
	}
      }

      //loop over the nHitsMaxPerPanel hits closest to the helix prediction
      double  panel_xdphi(0);
      double  counter_panel_hits(0);
      for (int t=0; t<_nHitsMaxPerPanel; ++t){
	if (panelHitIndex[t] < 0)                          continue;
	hit                = &panelz->fHitData.at(panelHitIndex[t]);
	panel_xdphi        = panel_xdphi + panelHitChi2[t];
	counter_panel_hits = counter_panel_hits + 1;

	Helix._szphi.addPoint(z,hit->_hphi,hit->_zphiWeight);
	Helix._nZPhiSh    += hit->nStrawHits();
	++count;
      }
      panel_xdphi = panel_xdphi/counter_panel_hits;

      if (count == 1) {//FIXME! investigate if it is needed or not
	zlast = z;
	dz    = 0.;
      }
	
	
      if ( (count>=5) &&      //FIXME! should we count and/or the number of panels?
	   // (xdphi < 2.)){
	   (panel_xdphi < 2.)){
	if ( (fabs(dfdz - Helix._szphi.dfdz()) < 8.e-4) ){//  || //require that the new value of dfdz is
	  //close to the starting one. update dfdz only if:
	  if ( (Helix._szphi.dfdz() > 0.) && //{                    // 1. the points browsed are more the half
	       (dz >=_mindist ) ){
	    phi0  = Helix._szphi.phi0();                     // 2. and require dfdz to be positivie! scattered hits or
	    dfdz  = Helix._szphi.dfdz();                     //    delta hits could have moved dfdz to negative value!
	    zlast = z;
	  }
	}
      }
	
    }
    _phiCorrectedDefined = 1;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ:BEFORE_CLEANUP] Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f points removed = %4i\n",
	     Helix._szphi.phi0(),Helix._szphi.dfdz(), Helix._szphi.chi2DofLine(), nPointsRemoved);
    }
    //-----------------------------------------------------------------------------
    // perform a cleanup in RZ
    //-----------------------------------------------------------------------------
    if ( DoCleanUp == 1){
      if ( Helix._szphi.chi2DofLine() > _chi2zphiMax) {
      NEXT_ITERATION:;
	//reset the coordinates of the worst hit
	iworst.Panel         = -1;
	iworst.PanelHitIndex = -1;

	indexWorst  = -1;
	chi2min     = 1e10;
	weightWorst = -1;

	for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
	  panelz = &Helix._oTracker[p];
	  z      = panelz->z;
	  int       nhits  = panelz->fNHits;
	  int       seedPanelIndex(0);
	  if (nhits == 0)                                        continue;
	  if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

	  for (int i=seedPanelIndex; i<nhits; ++i){
	    CalHelixPoint* hit = &panelz->fHitData.at(i);
	    int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	    if (Helix._hitsUsed[index] != 1)                     continue; 

	    szphi.init(Helix._szphi);

	    phi      = hit->_hphi;

	    if (UseInteligentWeight == 1){
	      weight = hit->_zphiWeight;
	    }

	    szphi.removePoint(z, phi, weight);
	    chi2 = szphi.chi2DofLine();
	    //printf("[CalHelixFinderAlg::doLinearFitPhiZ] chi2 = %5.3e chi2min = %5.3e\n", chi2, chi2min);
	    if (chi2 < chi2min) {
	      iworst.Panel         = p;
	      iworst.PanelHitIndex = i;

	      chi2min              = chi2;
	      weightWorst          = weight;
	    }
	  }
	}
	if ((iworst.Panel >= 0) && (Helix._nZPhiSh > _minNHits)) {
	  panelz  = &Helix._oTracker[iworst.Panel];
	  hit     = &panelz->fHitData.at(iworst.PanelHitIndex);
	  int index = iworst.Panel*CalHelixFinderData::kNMaxHitsPerPanel + iworst.PanelHitIndex;
	  Helix._hitsUsed[index] = 0;

	  z   = panelz->z;
	  phi = hit->_hphi;

	  Helix._szphi.removePoint(z, phi, hit->_zphiWeight);
	  Helix._nZPhiSh -= hit->nStrawHits();
	  
	  chi2min = Helix._szphi.chi2DofLine();
	  if (_debug > 5) {
	    printf("[CalHelixFinderAlg::doLinearFitPhiZ_removed:LOOP2] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);//FIXME! remove indexworst
	  }
	}

      CHECK_RESIDUALS:;
	dphi_max    = _maxXDPhi;
	//reset the coordinates of the worst hit
	iworst.Panel         = -1;
	iworst.PanelHitIndex = -1;

	weightWorst = -1;

	for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
	  panelz = &Helix._oTracker[p];
	  z      = panelz->z;
	  int       nhits  = panelz->fNHits;
	  int       seedPanelIndex(0);
	  if (nhits == 0)                                        continue;
	  if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

	  for (int i=seedPanelIndex; i<nhits; ++i){
	    CalHelixPoint* hit = &panelz->fHitData.at(i);
	    int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	    if (Helix._hitsUsed[index] != 1)                     continue;

	    pos      = hit->_pos;
	    strawDir = hit->_sdir;
	    phi      = z*Helix._szphi.dfdz() + Helix._szphi.phi0();
	    dphi     = fabs(hit->_hphi - phi);
	    err      = _sigmaPhi;

	    if (UseInteligentWeight == 1){
	      weight = calculatePhiWeight(*hit, /*pos, strawDir, */ helCenter, radius, 0, banner);//hit->_zphiWeight;//
	      err    = 1./sqrt(weight);
	    }

	    xdphi = dphi/err;

	    if ( xdphi > dphi_max) {
	      iworst.Panel         = p;
	      iworst.PanelHitIndex = i;

	      dphi_max             = xdphi;
	      weightWorst          = weight;
	    }
	  }
	}
	//remove the point
	if(iworst.Panel>=0 && Helix._nZPhiSh > _minNHits){
	  panelz  = &Helix._oTracker[iworst.Panel];
	  hit     = &panelz->fHitData.at(iworst.PanelHitIndex);

	  int index = iworst.Panel*CalHelixFinderData::kNMaxHitsPerPanel + iworst.PanelHitIndex;
	  Helix._hitsUsed[index] = 0;

	  z           = panelz->z;
	  phi         = hit->_hphi;
	  weightWorst = hit->_zphiWeight;

	  Helix._szphi.removePoint(z, phi, weightWorst);
	  Helix._nZPhiSh -= hit->nStrawHits();

	  chi2min = Helix._szphi.chi2DofLine();
	  if (_debug > 5) {
	    printf("[CalHelixFinderAlg::doLinearFitPhiZ:REMOVED] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);
	  }
	  goto CHECK_RESIDUALS;
	}

	if(Helix._nZPhiSh <= _minNHits) chi2min = Helix._szphi.chi2DofLine();
	
	if ( (chi2min >= _chi2zphiMax) ||
	     (iworst.Panel>=0 )) {

	  //--------------------------------------------------------------------------------
	  // 2016-04-27 gianipez: why should I not check the chi2 f I have 10 hits?
	  //--------------------------------------------------------------------------------
	  if (Helix._nZPhiSh > _minNHits) {
	    goto NEXT_ITERATION;
	  }
	}
      }
    }
    //-----------------------------------------------------------------------------
    // 2015-04-21 Gianipez changed the threshold from 3 to _minNHits. there is no reason
    // which should allow to keep a result which selects a number of points lower than the threshold!
    //-----------------------------------------------------------------------------
    if ( (Helix._nZPhiSh >= _minNHits) && (Helix._szphi.chi2DofLine() < _chi2zphiMax) ){
      success = true;
    }
    //----------------------------------------------------------------------//
    if (Helix._szphi.dfdz() < 0.) { // *FIXME* : negative helicity handling
      success = false;
    }
    else if (success) {                               // update helix results
      Helix._fz0  = Helix._szphi.phi0();
      Helix._dfdz = Helix._szphi.dfdz();
    }
    
    if ((SeedIndex.Panel == 0) && (SeedIndex.PanelHitIndex == 0) && (_diag > 0)) {
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
      Helix._diag.phi0_6           = Helix._szphi.phi0();
      Helix._diag.rdfdz_7          = Helix._szphi.dfdz()* Helix._diag.n_rescued_points_9;
      Helix._diag.dfdz_8           = Helix._szphi.dfdz();
      Helix._diag.chi2_dof_line_13 = Helix._szphi.chi2DofLine();

      if (success) {
	int h=0;

	for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
	  panelz = &Helix._oTracker[p];
	  int  nhitsPerPanel  = panelz->fNHits;
	  int  seedPanelIndex(0);
	  if (nhitsPerPanel == 0)                                                            continue;
	  if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

	  for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	    hit = &panelz->fHitData.at(i);
	    int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	    if (Helix._hitsUsed[index] != 1)                      continue; 
	    z        = panelz->z;//pos.z();
	    phi      = z* Helix._dfdz + Helix._fz0;
	    deltaPhi = hit->_hphi - phi;

	    if (h < Helix.maxIndex()) {
	      Helix._diag.resid[h] = deltaPhi;
	      ++h;
	    }
	    else {
	      printf (" ERROR: too many hits. Ignore \n");
	    }
	  }
	}
      }
    }

    if (_debug > 5) {

      printf("[CalHelixFinderAlg::doLinearFitPhiZ:END] retval = %d Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f szphi: phi_0 = %5.3f dfdz = %5.5f\n",
	     success ? 1:0, Helix._szphi.phi0(),Helix._szphi.dfdz(), Helix._szphi.chi2DofLine(), szphi.phi0(), szphi.dfdz() );

      if (_debug > 10) {
	printf("[CalHelixFinderAlg::doLinearFitPhiZ:END2]    flag   A   shID       z             phi      phi-dfdz*z-phi0\n");

	int lastPanel = SeedIndex.Panel;

	for (int p=CalHelixFinderData::kNTotalPanels-1; p>=lastPanel;  --p){
	  panelz = &Helix._oTracker[p];
	  int  nhitsPerPanel  = panelz->fNHits;
	  if (nhitsPerPanel == 0)                                                            continue;
	  if (p==SeedIndex.Panel) nhitsPerPanel = SeedIndex.PanelHitIndex;  

	  for (int  i=nhitsPerPanel-1;i>=0; --i){   
	    hit      = &panelz->fHitData.at(i);
	    z        = panelz->z;//pos.z();
	    phi      = z* Helix._dfdz + Helix._fz0;
	    deltaPhi = hit->_hphi - phi;

	    int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	    printf("[CalHelixFinderAlg::doLinearFitPhiZ:END2] %08x %2i %6i %12.5f %12.5f %12.5f\n",
		   *((int*) &hit->_flag), Helix._hitsUsed[index],
		   hit->strawId().straw()/*_strawhit->strawIndex().asInt()*/,  z, hit->_hphi, deltaPhi);
	  }
	}
      }
    }

    if (!success) {
      Helix._hitsUsed = hitsUsed;    
    }

    return success;
  }

//  void  CalHelixFinderAlg::fillHitLayer(CalHelixFinderData& Helix) {
    
// //--------------------------------------------------------------------------------
// // fill some geometrical info
// //--------------------------------------------------------------------------------
//     for (int ist=0; ist<_tracker->nStations(); ist++) {
//       const Station* st = &_tracker->getStation(ist);
// 
//       for (int ipl=0; ipl<st->nPlanes(); ipl++) {
// 	const Plane* pln = &st->getPlane(ipl);
// 	for (int ipn=0; ipn<pln->nPanels(); ipn++) {
// 	  const Panel* panel = &pln->getPanel(ipn);
// 	  for (int il=0; il<panel->nLayers(); il++) {
// 	    LayerZ_t* lz = &_hitLayer[il];
// 	    lz->fHitData.clear();
// 	    lz->fPanel = panel;
// //-----------------------------------------------------------------------------
// // panel caches phi of its center and the z
// //-----------------------------------------------------------------------------
// 	    lz->wx  = panel->straw0Direction().x();
// 	    lz->wy  = panel->straw0Direction().y();
// 	    lz->phi = panel->straw0MidPoint().phi();
// 	    lz->z   = (panel->getLayer(0).straw0MidPoint().z()+panel->getLayer(1).straw0MidPoint().z())/2.;
// 	  }
// 	}	
//       }
//     }
// 
// 
//     const vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();
// 
//     int size = Helix._timeCluster->nhits();
// //--------------------------------------------------------------------------------
//     if (Helix.shpos() != 0) {
//       int loc;
//       StrawHitFlag flag;
//       for (int i=0; i<size; ++i) {
// 	loc                = shIndices[i];	 // index in shcol of i-th timecluster hit
// 	flag               = Helix.shfcol()->at(loc);
// 	int good_hit = flag.hasAllProperties(_hsel  );
// 	int bkg_hit  = flag.hasAnyProperty  (_bkgsel);
// 	int used_hit = flag.hasAnyProperty  (StrawHitFlag::calosel);
// 
// 	if (good_hit && (! bkg_hit) && (! used_hit)) {
// 
// 	  const StrawHit*         sh    = &Helix.shcol()->at(loc);
// 	  const Straw*            straw = &_tracker->getStraw(sh->strawId());
// 	  const StrawHitPosition* shp   = &Helix.shpos()->at(loc);
//       
// 	  if (sh->energyDep() > _maxElectronHitEnergy)         continue;
// 
// 	  int       layerId  = straw->id().getLayer();
// 	  float     sigw     = shp->posRes(StrawHitPosition::wire);
// 	  
// 	  LayerZ_t* lz       = &_hitLayer[layerId];
// 
// 	  lz->fHitData.push_back(HitData_t(sh, shp, straw, sigw));
// 	}	
//       }
//     }
//     

//-----------------------------------------------------------------------------
// 12-09-2013 gianipez modified this procedure to avoid the doubling of the
// same stereohitposition
// points in filled array are ordered in Z coordinate
//-------------------------------------------------------------------------
  void CalHelixFinderAlg::fillXYZP(CalHelixFinderData& Helix) {

    static const double pi(M_PI);
    static const double twopi(2*pi);

    double clPhi(-9999.);

    if (fCaloTime > 0) clPhi = polyAtan2(fCaloY,fCaloX);

    const vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();
    ChannelID cx, co;

    int     size           = Helix._timeCluster->nhits();
    int     nFiltPoints(0), nFiltStrawHits(0);
    // int nTotalStations = _tracker->nStations();
    //--------------------------------------------------------------------------------
    // if (Helix.shpos() != 0) {
    int loc;
    StrawHitFlag flag;
    for (int i=0; i<size; ++i) {
      loc          = shIndices[i];	 // index in chcol of i-th timecluster hit
      flag         = Helix.shfcol()->at(loc);
      //-----------------------------------------------------------------------------
      // select hits: don't reuse straw hits
      //-----------------------------------------------------------------------------
      int good_hit = flag.hasAllProperties(_hsel  );
      int bkg_hit  = flag.hasAnyProperty  (_bkgsel);
      int used_hit = flag.hasAnyProperty  (StrawHitFlag::calosel);

      if (good_hit && (! bkg_hit) && (! used_hit)) {

	const ComboHit& ch          = Helix.chcol()->at(loc);

	if (ch.energyDep() > _maxElectronHitEnergy)         continue;

	//skip the hit if it doesn't rely on the semi-plane where the calo-lcuster is
	if (_filter) {
	  double chPhi = polyAtan2(ch.pos().y(), ch.pos().x());
	  double dphi  = chPhi - clPhi;

	  if (dphi >  pi) dphi -= twopi;
	  if (dphi < -pi) dphi += twopi;
	    
	  if (fabs(dphi) > pi/2)                            continue;
	}

	cx.Station                 = ch.strawId().station();//straw.id().getStation();
	cx.Plane                   = ch.strawId().plane() % 2;//straw.id().getPlane() % 2;
	cx.Face                    = ch.strawId().face();
	cx.Panel                   = ch.strawId().panel();//straw.id().getPanel();

	// get Z-ordered location
	Helix.orderID(&cx, &co);
     
	int os       = co.Station; 
	int of       = co.Face;
	int op       = co.Panel;

	int       stationId = os;
	int       faceId    = of + stationId*CalHelixFinderData::kNFaces;
	int       panelId   = op + faceId*CalHelixFinderData::kNPanelsPerFace;
	PanelZ_t* pz        = &Helix._oTracker[panelId];

	if ((os < 0) || (os >= CalHelixFinderData::kNStations     )) printf(" >>> ERROR: wrong station number: %i\n",os);
	if ((of < 0) || (of >= CalHelixFinderData::kNFaces        )) printf(" >>> ERROR: wrong face    number: %i\n",of);
	if ((op < 0) || (op >= CalHelixFinderData::kNPanelsPerFace)) printf(" >>> ERROR: wrong panel   number: %i\n",op);

	// pz->fHitData.push_back(CalHelixPoint(loc,sh,shp,straw,flag));
	pz->fHitData.push_back(CalHelixPoint(loc,ch,flag));
	pz->fNHits  = pz->fNHits + 1;
	if (pz->fNHits > CalHelixFinderData::kNMaxHitsPerPanel) printf("[CalHelixFinderAlg::fillXYZP] number of hits with the panel exceed the limit: NHits =  %i MaxNHits = %i\n", pz->fNHits, CalHelixFinderData::kNMaxHitsPerPanel);
	++nFiltPoints;
	nFiltStrawHits += ch.nStrawHits();
      }
    }
    // }
    
    Helix._nFiltPoints    = nFiltPoints;     //ComboHit counter
    Helix._nFiltStrawHits = nFiltStrawHits;  //StrawHit counter

    if (_debug > 0) printXYZP(Helix);
  }




//----------------------------------------------------------------------------
//2015-01-17 G. Pezzullo: the following procedure looks the hit with
// z-coordinate smaller then the seeding one and calculates distance from
// prediction in order to check if they are good or outliers
//----------------------------------------------------------------------------
  void CalHelixFinderAlg::rescueHitsBeforeSeed(CalHelixFinderData& Helix){
    PanelZ_t*   panelz = &Helix._oTracker[Helix._seedIndex.Panel];

    double      weight(-1), radius, phi0, dfdz, x0, y0;
    dfdz        = Helix._dfdz;
    //    phi0        = Helix._fz0 + dfdz*(_xyzp[Helix._seedIndex]._pos.z());
    phi0        = Helix._fz0 + dfdz*(panelz->z);
    x0          = Helix._center.x();
    y0          = Helix._center.y();
    radius      = Helix._radius;

    double      dx,dy,phi,max_dist;
    XYZVec      shPos, hePos, strawDir, helCenter(x0, y0, 0);

    double      deltaZ(0.); // , deltaX(0.), deltaY(0.);
    double      distXY(0.);
    double      dist(0.), dist2(0.); // help parameter for storing strawhit position residual
    int         rescuedStrawHits(0), rescuedPoints(0);

    char banner[]="CalHelixFinderAlg::rescueHitsBeforeSeed";

    if (_debug > 0) {
      printf("[%s:BEGIN] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.5f dfdz = %5.6f chi2 = %5.3f \n", banner,
	     x0, y0, radius, phi0, dfdz , Helix._sxy.chi2DofCircle());
      printf("[%s] SeedIndex = %i N-points = %5.3f\n",  banner, Helix._seedIndex.Panel, Helix._sxy.qn()-1);//FIXME!
      if (Helix._seedIndex.Panel >= 0) {//FIXME!
	printf("[%s] index      Z        xi      yi       xp       yp       X0        Y0         R        dfdZ  dXY(pred) dXY(seed) dZ(seed)     wt\n",banner);
	printf("[%s]-------------------------------------------------------------------------------------------------------------------------------\n",banner);
      }
    }
//-----------------------------------------------------------------------------
// given a helix candidate, move upstream and pick up points with Z < _xyzp[fSeedIndex].z
//-----------------------------------------------------------------------------
    PanelZ_t*      seedPanelz = &Helix._oTracker[Helix._seedIndex.Panel];
    PanelZ_t*      lastPanelz = &Helix._oTracker[Helix._seedIndex.Panel];
    CalHelixPoint* hit(0);
    // int            nLayers(2);
    float          panelHitChi2 [_nHitsMaxPerPanel] = {1e10};
    int            panelHitIndex[_nHitsMaxPerPanel] = {-1};

    for (int p=Helix._seedIndex.Panel; p>=0;  --p){
      panelz = &Helix._oTracker[p];
      int  nhitsPerPanel  = panelz->fNHits;

      if (nhitsPerPanel == 0)                                              continue;
      if (p==Helix._seedIndex.Panel) nhitsPerPanel = Helix._seedIndex.PanelHitIndex;//the seedHit is already clusterized!      

      //-----------------------------------------------------------------------------
      // dfdz = tanLambda/radius; phi0 is the last found hit phi
      //-----------------------------------------------------------------------------
      deltaZ    = panelz->z - lastPanelz->z;
      phi       = phi0 + (deltaZ)*dfdz;     
      //evaluate the helix prediction using the z coordinate of the panel
      hePos     = XYZVec(x0 + radius*std::cos(phi),
			  y0 + radius*std::sin(phi),
			  panelz->z);
      //check the Panel-phi wrt to the DS center
      // double  hePosPhi = TVector2::Phi_0_2pi(polyAtan2(hePos.y(), hePos.x()));
      double  hePosPhi = polyAtan2(hePos.y(), hePos.x());
      if (hePosPhi < 0) hePosPhi = hePosPhi + 2*M_PI;
      double  deltaPhi = hePosPhi - panelz->phi;
      if ( deltaPhi > M_PI ) deltaPhi -= 2*M_PI;
      if ( deltaPhi < -M_PI) deltaPhi += 2*M_PI;
      if ( fabs(deltaPhi) > _maxPanelToHelixDPhi)                             continue;
	
      // reset the chi2 values calulated in the previous panel
      // for (int l=0; l<nLayers; ++l){
      for (int g=0; g<_nHitsMaxPerPanel; ++g){
	panelHitChi2 [g] = 1e10;
	panelHitIndex[g] = -1;
      }
      // }
  
      for (int  i=nhitsPerPanel-1;i>=0; --i){   
	hit       = &panelz->fHitData.at(i);
	shPos     = hit->_pos;
	strawDir  = hit->_sdir;

	dx        = hePos.x() - shPos.x();
	dy        = hePos.y() - shPos.y();
	dist2     = dx*dx + dy*dy;
	dist      = std::sqrt(dist2);

	if (_debug > 10) {
	  printf("[%s:LOOP] %5i %9.3f %8.3f %8.3f %8.3f %8.3f %9.3f %9.3f %9.3f %8.5f %8.3f %8.3f %8.3f %8.3f",
		 banner,i,shPos.z(),shPos.x(),shPos.y(),hePos.x(),hePos.y(),
		 x0,y0,radius,dfdz,dist,distXY,deltaZ,weight);
	}

	max_dist = _distPatRec + _dfdzErr*fabs(deltaZ);
	if (dist <= max_dist) {
	  double delta_min(0);
	  int    index_min(-1);
	  // int    layer_id(hit->strawId().layer());
	  for (int k=0; k<_nHitsMaxPerPanel; ++k){
	    double delta = dist - panelHitChi2[k];
	    if (delta < delta_min){
	      delta_min = delta;
	      index_min = k;
	    }
	  }
	  if (index_min>=0) {
	    panelHitIndex[index_min] = i;
	    panelHitChi2 [index_min] = dist;
	  }
	}else {
	  if (_debug > 10) {
	    printf(" missed\n");
	  }
	}
      }//end loop pver the hits within the panel

      //loop over the nHitsMaxPerPanel hits closest to the helix prediction
      for (int t=0; t<_nHitsMaxPerPanel; ++t){
	int hitIndex = panelHitIndex[t];
	if (hitIndex < 0)                        continue;
	hit        = &panelz->fHitData.at(hitIndex);
	  
	lastPanelz = &Helix._oTracker[p];

	// add point to the helixfithack result objet
	weight     = calculateWeight(*hit, helCenter, radius);
	Helix._sxy.addPoint(hit->x(),hit->y(),weight);
	Helix._nXYSh += hit->nStrawHits();

	// update helix parameters
	x0      = Helix._sxy.x0();
	y0      = Helix._sxy.y0();
	radius  = Helix._sxy.radius();

	helCenter.SetX(x0);
	helCenter.SetY(y0);
	double dx  = (hit->x() - helCenter.x());
	double dy  = (hit->y() - helCenter.y());
	phi0       =  polyAtan2(dy,dx);

	//update hit info
	hit->_xyWeight   = weight;
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + hitIndex;
	Helix._hitsUsed[index] = 1;

	hit->_hphi        = phi0;
	hit->_dzFromSeed = deltaZ;
	hit->_drFromPred = panelHitChi2[t];
	double dzFromSeed = panelz->z - seedPanelz->z;         // expected to be negative
	hit->_dzFromSeed = dzFromSeed;
	hit->_drFromPred = panelHitChi2[t];  

	rescuedStrawHits += hit->nStrawHits();
	++rescuedPoints;

	if (_debug > 0) {
	  printf("rescued %08x %2i %12.5f %12.5f %12.5f\n",
		 *((int*) &hit->_flag), Helix._hitsUsed[index],
		 hit->x(), hit->y(), hit->z());
	}
      }
    }
    //-----------------------------------------------------------------------------
    // update Helix info
    //-----------------------------------------------------------------------------
    Helix._center.SetXYZ(x0, y0, 0.0);
    Helix._radius      = radius;
    Helix._nComboHits += rescuedPoints;
    Helix._nStrawHits += rescuedStrawHits;

    if (_debug > 5) {
      printf("[%s:END] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.3f dfdz = %5.6f chi2 = %5.3f",
	     banner,
	     x0, y0, radius,  Helix._fz0, dfdz , Helix._sxy.chi2DofCircle());
      printf(" SeedIndex: %i N(rescued points): %i\n",Helix._seedIndex.Panel,rescuedPoints);//FIXME!
    }

    if (_diag) Helix._diag.n_rescued_points_16 = rescuedPoints;
  }

//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::printInfo(CalHelixFinderData&  Helix){
    const char banner [] = "CalHelixFinderAlg::printInfo";
    double dr(0), dx, dy;

    if (_debug > 0) {
      printf("[%s] N(points): %3i x0: %12.5f y0: %12.5f r: %12.5f chi2c: %12.5f phi0: %5.5f dfdz: %5.6f chi2l: %5.3f\n",
	     banner,
	     Helix._nStrawHits,
	     Helix._sxy.x0(),Helix._sxy.y0(),Helix._sxy.radius(),Helix._sxy.chi2DofCircle(),
	     Helix._fz0, Helix._dfdz , Helix._szphi.chi2DofLine());

      PanelZ_t*      panelz(0);
      CalHelixPoint* hit(0);

      for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
	panelz = &Helix._oTracker[p];
	int  nhits          = panelz->fNHits;
	for (int i=0; i<nhits; ++i){   
	  hit = &panelz->fHitData.at(i);
	  dx = hit->x() - Helix._sxy.x0();
	  dy = hit->y() - Helix._sxy.y0();
	  dr = sqrt(dx*dx+dy*dy) - Helix._sxy.radius();
	  printf("[%s] %08x %6i %3i %6i %12.5f %12.5f %12.5f %10.3f\n",banner,
		 *((int*) &hit->_flag),  int(hit->index()), i, hit->strawId().straw(),
		 hit->x(), hit->y(), hit->z(), dr
		 );//FIXME!
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// this routine simply checks '_indicesTrkCandidate' array and for negative
// indices sets the 'outlier' flag to the corresponding 'xyzp'
// no actual check of residuals is performed
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::filterUsingPatternRecognition(CalHelixFinderData& Helix) {

    if (Helix._seedIndex.Panel < 0) return;

    int            nActive(0), nActive_hel(0);
    int            nSh = Helix._nFiltStrawHits;

    double         straw_mean_radius(0), chi2_global_helix(0), total_weight(0);
    double         x_center(Helix._center.x()), y_center(Helix._center.y()), radius(Helix._radius);
    double         fz0(Helix._fz0), lambda(1./Helix._dfdz);
    XYZVec         hel_pred(0., 0., 0.);

    PanelZ_t*      panelz(0);
    CalHelixPoint* hit(0); bool isFirst(true);

    for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhits     = panelz->fNHits;
      if (nhits == 0)             continue;
	
      for (int i=0; i<nhits; ++i){   
	hit = &panelz->fHitData.at(i);
	if (_debug > 0) {
	  if (isFirst) {
	    isFirst = false;
	    printf("[CalHelixFinderAlg::filterUsingPatternRecognition]  filterUsingPatternRecognition() will set asOutlier the following hits using helix parameters\n");
	    printf("[CalHelixFinderAlg::filterUsingPatternRecognition] X0 = %5.3f Y0 = %5.3f r = %5.3f chi2N = %5.5f phi0 = %5.5f dfdz = %5.5f chi2N = %5.5f straw-hits = %i\n",
		   Helix._sxy.x0(), Helix._sxy.y0(), Helix._radius, Helix._sxy.chi2DofCircle(),
		   Helix._szphi.phi0(), Helix._szphi.dfdz(), Helix._szphi.chi2DofLine(),
		   Helix._nStrawHits);//goodPointsTrkCandidate );// +1 for counting also the seeding hit
	    printf("[CalHelixFinderAlg::filterUsingPatternRecognition] index  shID type           X        Y         Z        dist      Dz\n");
	    printf("[CalHelixFinderAlg::filterUsingPatternRecognition] -------------------------------------------------------------------\n");
	  }
	}

	double dist = hit->_drFromPred;
	double dz   = hit->_dzFromSeed;

	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1)        hit->setOutlier(); // *FIXME* ? ? do we want to call this an outlier ?
	else  {
	  // ++nActive;
	  nActive          += hit->nStrawHits();
	  double  x         = hit->x();
	  double  y         = hit->y();
	  double  z         = panelz->z;
	  double  phi_pred  = fz0 + z/lambda;
	  double  x_pred    = x_center + radius*cos(phi_pred);
	  double  y_pred    = y_center + radius*sin(phi_pred);
	  hel_pred.SetX(x_pred);
	  hel_pred.SetY(y_pred);
	  double  weight    = calculateWeight(*hit, hel_pred, radius)*_weight3D;
	  double  x_resid2  = (x - x_pred)*(x - x_pred);
	  double  y_resid2  = (y - y_pred)*(y - y_pred);
	  double  hitResi2  = (x_resid2 + y_resid2)*weight;
	  if (hitResi2 > _chi2hel3DMax) {
	    hit->setOutlier();  // *FIXME* ? ? do we want to call this an outlier ?
	  } else {
	    // ++nActive_hel;
	    nActive_hel      += hit->nStrawHits();
	    chi2_global_helix = chi2_global_helix + hitResi2;
	    straw_mean_radius = straw_mean_radius + sqrt(x*x + y*y)*weight;
	    total_weight      = total_weight + weight;	
	  }
	}

	if (_debug > 10) {
	  XYZVec*     shPos = &hit->_pos;
	  int         is_outlier    = hit->isOutlier();
	  string      type;
	  if      ((p == Helix._seedIndex.Panel) && (i == Helix._seedIndex.PanelHitIndex)) type = "seed";
	  else if ((p == Helix._candIndex.Panel) && (i == Helix._candIndex.PanelHitIndex)) type = "cand";

	  printf("[CalHelixFinderAlg::filterUsingPatternRecognition] %5i %5i %4i %4s  %8.3f %8.3f %9.3f %8.3f %8.3f\n",
		 i,hit->strawId().straw(),is_outlier,type.data(),shPos->x(),shPos->y(),shPos->z(),dist,dz);
	}
      }
    }
    if (_debug > 5) {
      printf("[CalHelixFinderAlg::filterUsingPatternRecognition:END] N( active StrawHit) = %i over N-StrawHits = %i\n", nActive, nSh);
    }

    Helix._nStrawHits = nActive_hel; 
 
    if (_diag){
      Helix._diag.n_active_11       = nActive_hel;
      Helix._diag.straw_mean_radius = (nActive > 0) ? straw_mean_radius/total_weight : -1;
      Helix._diag.chi2d_helix       = (nActive > 0) ? chi2_global_helix/double(nActive_hel - 5.): -1;
    }
  }

//-----------------------------------------------------------------------------
  double CalHelixFinderAlg::deltaPhi(double phi1, double phi2){
    static const double pi(M_PI);
    static const double twopi(2*pi);
    double dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }


//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::searchBestTriplet   (CalHelixFinderData& Helix, CalHelixFinderData& TmpHelix, int UseMPVdfdz){
    int       nSh = Helix._nFiltStrawHits;
    int       nHitsTested(0);

    PanelZ_t* panelz;
    
    for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int       nhits  = panelz->fNHits;
      for (int i=0; i<nhits; ++i){   
	if (Helix._nStrawHits > (nSh - nHitsTested))   continue;	  
	if ((nSh - nHitsTested) < _minNHits     )   continue;	  
	//clear the info of the tmp object used to test the triplet
	TmpHelix.clearResults();

	SeedInfo_t          seed(p,i);
	findTrack(seed,TmpHelix,UseMPVdfdz);
	// ++nHitsTested;
	nHitsTested += panelz->fHitData.at(i).nStrawHits();

	//compare tripletHelix with bestTripletHelix
	if (( TmpHelix._nStrawHits >  Helix._nStrawHits) ||
	    ((TmpHelix._nStrawHits == Helix._nStrawHits) && (Helix._helixChi2 < Helix._helixChi2))) {
	  Helix = TmpHelix;
	}
	if (_debug > 5) {
	  printf("[CalHelixFinderAlg::doPatternRecognition]: calling findTrack(i=%i,Helix,useDefaltDfDz=FALSE,useMPVdfdz=%i)",i,UseMPVdfdz);
	  printf(" : np=%3i _goodPointsTrkCandidate=%3i\n",nSh,Helix._nStrawHits);
	}
      }
    }
    
  }


//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::doPatternRecognition(CalHelixFinderData& Helix) {

    int    useMPVdfdz(1), useIntelligentWeight(1);//, nHitsTested(0);

    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition:BEGIN] fUseDefaultDfDz = %i\n",fUseDefaultDfDz);

    if (_debug2 == 0){
      _debug2 = _debug;
      _debug  = 0;
    }

    CalHelixFinderData tripletHelix(Helix);
    tripletHelix._helix = NULL;//FIXME!

    _findTrackLoopIndex = 1; 		// debugging
    searchBestTriplet(Helix, tripletHelix);
    //-----------------------------------------------------------------------------
    // 2014-11-09 gianipez: if no track was found requiring the recalculation of dfdz
    // look for a track candidate using the default value of dfdz and the target center
    //-----------------------------------------------------------------------------
    _findTrackLoopIndex = 2; 		// *DEBUGGING*
    if (fUseDefaultDfDz == 0) {
      searchBestTriplet(Helix, tripletHelix, useMPVdfdz);
   }

    if (_debug == 0){
      _debug  = _debug2;
      _debug2 = 0;
    }

    char banner[200];
    bool rc;
    int  rs, usePhiResid;


    if ((Helix._seedIndex.Panel < 0) || (Helix._nXYSh < _minNHits) ) goto  PATTERN_RECOGNITION_END;

    // 2015-01-17 G. Pezzullo: rescue points with z-coordinate less than the seed hit
    if (_debug != 0) {
      printf("[CalHelixFinderAlg::doPatternRecognition]: calling rescueHitsBeforeSeed\n");
      printInfo(Helix);
    }

    rc = doLinearFitPhiZ(Helix, SeedInfo_t(0,0), useIntelligentWeight);

    //2017-10-05 Gianipez added the following line to make some tests
    if (Helix._szphi.qn() == 0.)                                 goto  PATTERN_RECOGNITION_END;

    rescueHitsBeforeSeed(Helix);
//-----------------------------------------------------------------------------
// finally, assume that the found helix is already close enough and refine
// the helix parameters accounting for different weights
//-----------------------------------------------------------------------------
    if (_debug != 0)  printInfo(Helix);
//--------------------------------------------------------------------------------
// 2017-09-25 gianipez
// ::refineHelixParameters uses _sxy, so we need to update it using the  result
// from the inteligent-weight fit _sxyw
//--------------------------------------------------------------------------------
    strcpy(banner,"refineHelixParameters");

    refineHelixParameters(Helix, SeedInfo_t(0,0), banner, _debug);

     if (_debug != 0)  printInfo(Helix);
//---------------------------------------------------------------------------------------
// use the results of the helix search to see if points along the track can be rescued
//---------------------------------------------------------------------------------------
    if ((Helix._nZPhiSh < _minNHits) || (!rc)) usePhiResid = 0;
    else                                       usePhiResid = 1;

    rescueHits(Helix, SeedInfo_t(0,0), usePhiResid);
    
    if ((Helix._nXYSh - 1) != Helix._nZPhiSh) rc = doLinearFitPhiZ(Helix, SeedInfo_t(0,0), useIntelligentWeight);//the factor "-1" takes into account that the XY fit includes the target center

    if (_debug != 0)  printInfo(Helix);
//--------------------------------------------------------------------------------------------------------------
// 2015-03-25 G. Pezzu added the following call to findDfDz(...) in order to help the fitter on finding
// the more reliable value of dfdz which is needed for resolving the 2pi ambiguity.
// Since in the previous step we could have rescued few points, that would give us an help!
// re-evaluate the df/dz and phi0 including rescued hits and new XY parameters
//--------------------------------------------------------------------------------------------------------------
    if (Helix._nZPhiSh < _minNHits || (!rc)){
      rs = findDfDz(Helix, SeedInfo_t(0,0));
      
      if (rs == 1) {			// update Helix Z-phi part
	Helix._dfdz = _hdfdz;
	Helix._fz0  = _hphi0;
      }
    }

    rc = doLinearFitPhiZ(Helix, SeedInfo_t(0,0), useIntelligentWeight);

    if (rc) {
      usePhiResid = 1;
      rescueHits(Helix, SeedInfo_t(0,0), usePhiResid);
      if ((Helix._nXYSh - 1) != Helix._nZPhiSh) rc = doLinearFitPhiZ(Helix, SeedInfo_t(0,0), useIntelligentWeight);  //the factor "-1" takes into account that the XY fit includes the target center

      if (_debug != 0)  printInfo(Helix);
      strcpy(banner,"refineHelixParameters-after-doLinearFitPhiZ");
      refineHelixParameters(Helix,SeedInfo_t(0,0),banner,_debug);
      if (_debug != 0)  printInfo(Helix);
    }
//-----------------------------------------------------------------------------
// more diagnostic data, evaluate helix radius without using the target center
//-----------------------------------------------------------------------------
    if (_diag > 0) {
      Helix._diag.radius_14          = Helix._radius;
      Helix._diag.chi2_dof_circle_15 = Helix._sxy.chi2DofCircle();
      Helix._diag.z0_6               = Helix._fz0;
      Helix._diag.rdfdz_7            = Helix._dfdz*Helix._radius;
      Helix._diag.dfdz_8             = Helix._dfdz;

      Helix._sxy.removePoint(0., 0., 1./900.);
      Helix._diag.dr                 = Helix._radius - Helix._sxy.radius();
      Helix._sxy.addPoint(0., 0., 1./900.);
    }
//-----------------------------------------------------------------------------
// 2014-11-09 gianipez changed the cleanup process. now it is faster and cleaner
//-----------------------------------------------------------------------------
    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition]: calling filterUsingPatternRecognition\n");
    filterUsingPatternRecognition(Helix);

  PATTERN_RECOGNITION_END:;
//-----------------------------------------------------------------------------
// if running in the diagnostics mode, save state of the Xyzp (this is a deep copy)
//-----------------------------------------------------------------------------
    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition:END]\n");
  }

//--------------------------------------------------------------------------------
// a straw man attempt to account for significantly different resolutions 
// along the wire and in the drift direction
//--------------------------------------------------------------------------------
  double  CalHelixFinderAlg::calculateWeight(const CalHelixPoint& Hit,
					     // const Hep3Vector& HitPos   ,
					     // const Hep3Vector& StrawDir ,
					     const XYZVec&        HelCenter,
					     double               Radius   ) {

    // double    rs(2.5);   // straw radius, mm
    double    transErr = 5./sqrt(12.);
    //scale the error based on the number of the strawHits that are within teh ComboHit
    if (Hit.nStrawHits() > 1) transErr *= 1.5;
    double    transErr2 = transErr*transErr;

    double x   = Hit.x();
    double y   = Hit.y();
    double dx  = x-HelCenter.x();
    double dy  = y-HelCenter.y();
    double dxn = dx*Hit._sdir.x()+dy*Hit._sdir.y();

    double costh2 = dxn*dxn/(dx*dx+dy*dy);
    double sinth2 = 1-costh2;

    // double e2     = _ew*_ew*sinth2+rs*rs*costh2;
    double e2     = Hit.wireErr2()*sinth2+transErr2*costh2;
    double wt     = 1./e2;
                                                    // scale the weight for having chi2/ndof distribution peaking at 1
    wt *= _weightXY;

    return wt;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  double  CalHelixFinderAlg::calculatePhiWeight(const CalHelixPoint& Hit,
						// const XYZVec&  HitPos   ,
						// const XYZVec&  StrawDir ,
						const XYZVec&  HelCenter,
						double             Radius   ,
						int                Print    ,
						const char*        Banner   ) {
    // double    rs(2.5);  // straw radius, mm
    double    transErr = 5./sqrt(12.);
    //scale the error based on the number of the strawHits that are within teh ComboHit
    if (Hit.nStrawHits() > 1) transErr *= 1.5;
    double    transErr2 = transErr*transErr;

    double x  = Hit.x();
    double y  = Hit.y();
    double dx = x-HelCenter.x();
    double dy = y-HelCenter.y();

    double dxn = dx*Hit._sdir.x()+dy*Hit._sdir.y();

    double costh2  = dxn*dxn/(dx*dx+dy*dy);
    double sinth2 = 1-costh2;

    // double e2     = _ew*_ew*costh2+rs*rs*sinth2;
    double e2     = Hit.wireErr2()*costh2+transErr2*sinth2;
    double wt     = Radius*Radius/e2;
    wt           *= _weightZPhi;

    if (Print > 0) {
      double dr = calculateRadialDist(Hit.pos(),HelCenter,Radius);
      printf("[CalHelixFinderAlg::%s] %9.3f %9.3f %10.5f %10.5f %10.5f %10.5f %12.5e %10.3f\n", 
	                       Banner, x, y, dx, dy, costh2, sinth2, e2, dr);
    }

    return wt;
  }

//--------------------------------------------------------------------------------
// calculate the radial distance of a straw hit from the helix prediction
//--------------------------------------------------------------------------------
  double  CalHelixFinderAlg::calculateRadialDist (const XYZVec& HitPos   ,
						  const XYZVec& HelCenter,
						  double            Radius   ) {
    double dx = HitPos.x()-HelCenter.x();
    double dy = HitPos.y()-HelCenter.y();
    double dr = sqrt(dx*dx+dy*dy)-Radius;

    return dr;
  }


//-----------------------------------------------------------------------------
  void   CalHelixFinderAlg::doWeightedCircleFit (CalHelixFinderData& Helix,
						 SeedInfo_t          SeedIndex,
						 XYZVec&             HelCenter,
						 double&             Radius   ,
						 int                 Print    ,
						 const char*         Banner   ) {
    double     wt;
//-----------------------------------------------------------------------------
// add calorimeter cluster with a position error of 10 mm => wt = 1/100
//-----------------------------------------------------------------------------
    Helix._sxy.clear();
    Helix._nXYSh = 0;
    Helix._nComboHits = 0;

    Helix._sxy.addPoint(fCaloX,fCaloY,1./100.);
    Helix._nXYSh += 1;
    Helix._nComboHits += 1;
//-------------------------------------------------------------------------------
// add stopping target center with a position error of 100 mm/sqrt(12) ~ 30mm => wt = 1/900
//-------------------------------------------------------------------------------
    Helix._sxy.addPoint(0.,0.,1./900.);
    Helix._nXYSh += 1;
    Helix._nComboHits += 1;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doWeightedCircleFit] BEGIN: x0 = %8.3f y0 = %8.3f radius = %8.3f chi2dof = %8.3f\n",
	     HelCenter.x(),HelCenter.y(),Radius,Helix._sxy.chi2DofCircle());
      if (_debug > 10) {
	printf("[CalHelixFinderAlg::doWeightedCircleFit:LOOP] Index      X          Y         Z          wt        wireNx     wireNy\n");
      }
    }

    PanelZ_t*      panelz(0);
    CalHelixPoint* hit   (0);

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhits          = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhits == 0)                                                                       continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;
	
      for (int i=seedPanelIndex; i<nhits; ++i){   
	hit = &panelz->fHitData.at(i);
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1)                     continue;

	wt             = calculateWeight(*hit,HelCenter,Radius);
	hit->_xyWeight = wt;

	Helix._sxy.addPoint(hit->_pos.x(),hit->_pos.y(),wt);
	Helix._nXYSh      += hit->nStrawHits();
	Helix._nComboHits += 1;
    
	if (_debug > 10) {
	  printf("[CalHelixFinderAlg::doWeightedCircleFit:LOOP] %4i %10.3f %10.3f %10.3f %10.3e %10.4f %10.4f\n",
		 (int)hit->index(), hit->_pos.x(), hit->_pos.y(), hit->_pos.z(), wt, hit->_sdir.x(), hit->_sdir.y());
	}
      }
    }
					// update helix info
    Radius  = Helix._sxy.radius();
    HelCenter.SetX(Helix._sxy.x0());
    HelCenter.SetY(Helix._sxy.y0());

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doWeightedCircleFit:END] : npt = %3.0f  chi2dof = %8.3f x0 = %8.3f y0 = %8.3f radius = %8.3f\n",
	     Helix._sxy.qn(),Helix._sxy.chi2DofCircle(),HelCenter.x(),HelCenter.y(),Radius);
    }
  }


//-----------------------------------------------------------------------------
// this is a rather "primitive" definition of the worst hit, should do for now
//-----------------------------------------------------------------------------
  void    CalHelixFinderAlg::searchWorstHitWeightedCircleFit(CalHelixFinderData& Helix,
							     SeedInfo_t          SeedIndex,
							     const XYZVec&       HelCenter,
							     double&             Radius,
							     SeedInfo_t&         Iworst,
							     double&             HitChi2Worst)
  {
    HitChi2Worst         = _hitChi2Max;
    Iworst.Panel         = -1;
    Iworst.PanelHitIndex = -1;

    double     dr, hitChi2;
  
    CalHelixPoint* hit(0);
    PanelZ_t*      panelz(0);
    PanelZ_t*      seedPanelz = &Helix._oTracker[Helix._seedIndex.Panel];

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhitsPerPanel  = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhitsPerPanel == 0)                                                            continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

      for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	hit = &panelz->fHitData.at(i);
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1)                    continue;

	dr      = calculateRadialDist(hit->_pos,HelCenter,Radius);
	hitChi2 = dr*dr*hit->_xyWeight;

	// store info out the radial residual
	if ((SeedIndex.Panel == 0) && (SeedIndex.PanelHitIndex == 0)) {
	  hit->_drFromPred  = fabs(dr);//hitChi2;
	  double dzFromSeed = panelz->z - seedPanelz->z;              // expected to be positive (non-negative)
	  hit->_dzFromSeed  = dzFromSeed;
	}

	if (hitChi2 > HitChi2Worst) {
	  HitChi2Worst         = hitChi2;
	  Iworst.Panel         = p;
	  Iworst.PanelHitIndex = i;
	}
      }
    }
  }

//--------------------------------------------------------------------------------
// IWorst is always defined
// returns the index of the hit which provides the highest contribute to the chi2
//--------------------------------------------------------------------------------
  void    CalHelixFinderAlg::cleanUpWeightedCircleFit(CalHelixFinderData& Helix,
						      SeedInfo_t          SeedIndex,
						      SeedInfo_t&         IWorst)
  {
    LsqSums4   sxy;
    double     chi2, chi2_min (-1.), x, y;

    //reset the coordinates of the worst hit found previousl
    IWorst.Panel         = -1;
    IWorst.PanelHitIndex = -1;

    CalHelixPoint* hit(0);
    PanelZ_t*      panelz(0);
 
    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhitsPerPanel  = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhitsPerPanel == 0)                                                             continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

      for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	hit = &panelz->fHitData.at(i);
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1)                    continue;

	sxy.init(Helix._sxy);

	x  = hit->_pos.x();
	y  = hit->_pos.y();

	sxy.removePoint(x, y, hit->_xyWeight);

	chi2  = sxy.chi2DofCircle();

	if ((chi2 < chi2_min) || ( (i == SeedIndex.PanelHitIndex) && (p == SeedIndex.Panel))) {
	  chi2_min             = chi2;
	  IWorst.Panel         = p;
	  IWorst.PanelHitIndex = i;	    
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// use hits only, at this point the cluster is no longer needed 
//-----------------------------------------------------------------------------
  int CalHelixFinderAlg::refineHelixParameters(CalHelixFinderData& Trk,
					       SeedInfo_t          SeedIndex,
					       const char*         Banner,
					       int                 Print  ) {
    auto           hitsUsed = Trk._hitsUsed;
    double         x, y, r, r_start;
    double         hitChi2Worst;

    // ::LsqSums4     sxyw;
    int            pointsRemoved(0);

    SeedInfo_t     iworst(-1, -1);
    double         wtWorst;
    double         chi2, chi2_min;

    XYZVec         hitPos, strawDir, helCenter, helCenter_start;
    CalHelixPoint* hit(0);
    PanelZ_t*      panelz(0);
    int            hitUsedIndex(-1);

    int            rc(0);               // success-oriented initialization :)
					// initialize helix 
    r  = Trk._sxy.radius();
    helCenter.SetX( Trk._sxy.x0());
    helCenter.SetY( Trk._sxy.y0());

    helCenter_start = helCenter;
    r_start         = r;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::refineHelixParameters] BEGIN               x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	     Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle());
      printf("[CalHelixFinderAlg::refineHelixParameters] i       X        Y        dx        dy         costh        sinth2         e2     radial-dist\n");
    }

    doWeightedCircleFit (Trk,SeedIndex,helCenter,r,Print,Banner);
//-----------------------------------------------------------------------------
// recalcute weights using most recent helix parameters
//-----------------------------------------------------------------------------
    doWeightedCircleFit (Trk,SeedIndex,helCenter,r,Print,Banner);

    //now initialize the LsqSum4 variable
    // sxyw.init(Trk._sxy);

    searchWorstHitWeightedCircleFit(Trk,SeedIndex,helCenter,r,iworst,hitChi2Worst);

    chi2     = Trk._sxy.chi2DofCircle();
    chi2_min = chi2;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::refineHelixParameters] npt = %3.0f x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f iworst=%3i chi2Worst = %8.3f\n",
	     Trk._sxy.qn(),Trk._sxy.x0(),Trk._sxy.y0(),Trk._sxy.radius(),Trk._sxy.chi2DofCircle(),iworst.Panel,hitChi2Worst);
    }

    if ((chi2 <= _chi2xyMax) && (hitChi2Worst <= _hitChi2Max)) goto F_END;
//-----------------------------------------------------------------------------
// one of the chi2's is above the threshold, cleanup is needed
//-----------------------------------------------------------------------------
    if (_debug > 5) printf("[CalHelixFinderAlg::refineHelixParameters] : START CLEANUP\n");  
  NEXT_ITERATION:;

    cleanUpWeightedCircleFit(Trk,SeedIndex,iworst);

    if (iworst.Panel >= 0) {
      panelz  = &Trk._oTracker[iworst.Panel];
      hit     = &panelz->fHitData.at(iworst.PanelHitIndex);
      x       = hit->_pos.x();
      y       = hit->_pos.y();
      wtWorst = hit->_xyWeight;//weights[iworst];
      
					// remove point from the track, this is why need to return weights
      // sxyw.removePoint(x, y, wtWorst);
      Trk._sxy.removePoint(x, y, wtWorst);
      Trk._nXYSh -=hit->nStrawHits();

      hitUsedIndex = iworst.Panel*CalHelixFinderData::kNMaxHitsPerPanel + iworst.PanelHitIndex;
      Trk._hitsUsed[hitUsedIndex] = 0;

      // Trk._nPoints -= 1;
      Trk._nComboHits -= 1;
      Trk._nStrawHits -= hit->nStrawHits();
      ++pointsRemoved;

      if (_debug > 5) {
	printf("[CalHelixFinderAlg::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
	       Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(), _chi2xyMax, iworst.Panel);//FIXME!
      }
					// update helix parameters and refit
      r  = Trk._sxy.radius();
      helCenter.SetX(Trk._sxy.x0());
      helCenter.SetY(Trk._sxy.y0());

      doWeightedCircleFit (Trk,SeedIndex,helCenter,r,0,Banner);

					// update the chi2 value
      chi2_min = Trk._sxy.chi2DofCircle();
    }
//-----------------------------------------------------------------------------
// recalculate the worst radial residual
//-----------------------------------------------------------------------------
  CHECK_RESIDUALS: ;
    searchWorstHitWeightedCircleFit(Trk,SeedIndex,helCenter,r,iworst,hitChi2Worst);
//-----------------------------------------------------------------------------
// if a hit contributes chi2 > _hitCHi2Max, remove it and go back looking for the next such hit
//-----------------------------------------------------------------------------
    if (iworst.Panel >= 0) {
      panelz  = &Trk._oTracker[iworst.Panel];
      hit     = &panelz->fHitData.at(iworst.PanelHitIndex);
      x       = hit->_pos.x();
      y       = hit->_pos.y();
      wtWorst = hit->_xyWeight;
					// remove point from the track and mark it
      // sxyw.removePoint(x, y, wtWorst);
      Trk._sxy.removePoint(x, y, wtWorst);
      Trk._nXYSh -=hit->nStrawHits();

      hitUsedIndex = iworst.Panel*CalHelixFinderData::kNMaxHitsPerPanel + iworst.PanelHitIndex;
      Trk._hitsUsed[hitUsedIndex] = 0;

      // Trk._nPoints -= 1;
      Trk._nComboHits -= 1;
      Trk._nStrawHits -= hit->nStrawHits();
      ++pointsRemoved;

      if (_debug > 5) {
	printf("[CalHelixFinderAlg::refineHelixParameters:REMOVE] iworst=%3i (x0,y0,R) = (%8.3f, %8.3f, %8.3f) chi2 = %8.3f chi2Maxxy = %8.3f\n",
	       iworst.Panel, Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(), _chi2xyMax);//FIXME!
      }
					// update helix info
      r  = Trk._sxy.radius();
      helCenter.SetX(Trk._sxy.x0());
      helCenter.SetY(Trk._sxy.y0());
					// refit helix and update the chi2 value

      doWeightedCircleFit (Trk,SeedIndex,helCenter,r,0,Banner);
      chi2_min = Trk._sxy.chi2DofCircle();
                                                            goto CHECK_RESIDUALS;
    }

    if ((chi2_min >= _chi2xyMax) && (iworst.Panel >= 0)) {
//-----------------------------------------------------------------------------
// still bad chi2, repeat the cleanup cycle
//-----------------------------------------------------------------------------
      if (Trk._nStrawHits > _minNHits)                         goto NEXT_ITERATION;
      else  {
//-----------------------------------------------------------------------------
// number of points too small, bail out
//-----------------------------------------------------------------------------
	rc = -1;
      }
    }
    else if (chi2_min < _chi2xyMax) {
//-----------------------------------------------------------------------------
// 2018-01-03 PM: not sure what exactly this is supposed to mean, but...
//-----------------------------------------------------------------------------
      rc = 1;
    }

  F_END:;
    if (_debug > 5) {
      printf("[CalHelixFinderAlg::refineHelixParameters] END: RC=%1i", rc);
      printf(" npt:%3.0f x0: %8.3f y0: %8.3f R: %8.3f chi2: %9.3f N(removed):%3i\n",
	     Trk._sxy.qn(), Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(), pointsRemoved);
    }

    if (rc >= 0) {
//-----------------------------------------------------------------------------
// update circle parameters
//-----------------------------------------------------------------------------
      // Trk._sxyw.init(sxyw);
      Trk._center.SetX(Trk._sxy.x0());
      Trk._center.SetY(Trk._sxy.y0());
      Trk._radius = Trk._sxy.radius();
      Trk._chi2   = Trk._sxy.chi2DofCircle();

    }else {
      Trk._hitsUsed = hitsUsed;   //restore the info of the used-hits that was originally passed to the procedure      
      doWeightedCircleFit (Trk,SeedIndex,helCenter_start,r_start,0,Banner);
    }

    return rc;
  }


//-----------------------------------------------------------------------------
// The previous helix search may have  thrown away point which instead have
// small radial residual, so this function is devoted for rescueing these
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::rescueHits(CalHelixFinderData& Helix          ,
				     SeedInfo_t          SeedIndex      ,
				     int                 UsePhiResiduals) {

    const char  banner[] = "rescueHits";
    double      wt, e2, x, y, r;
    double      phiwt(-9999.);

    XYZVec  hitPos, strawDir, helCenter, hel_pred(0.,0.,0.);

    double      dfdz, phi0, dphi, dphiChi2(0.0), phi_pred;
    
    ::LsqSums4  sxy;
    int         n_added_points(0);
    SeedInfo_t  ibest(-1,-1);
    double      wtBest, phiwtBest;
    double      chi2, chi2_min, dr, hitChi2, drChi2;

    double      x_pred(0), y_pred(0), weight_hel(0), x_resid2(0), y_resid2(0);
 
    PanelZ_t*      panelz(0);
    PanelZ_t*      seedPanelz = &Helix._oTracker[Helix._seedIndex.Panel];
    CalHelixPoint* hit(0);

   //set  dfdz and phi0
    dfdz = Helix._dfdz;
    phi0 = Helix._fz0;

    //update helix info
    r  = Helix._sxy.radius();
    helCenter.SetX( Helix._sxy.x0());
    helCenter.SetY( Helix._sxy.y0());

    doWeightedCircleFit (Helix, SeedIndex, helCenter, r);

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::%s:BEGIN] x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f phi0 = %9.6f dfdz = %9.6f chi2 = %8.3f\n",
	     banner,
	     Helix._sxy.x0(), Helix._sxy.y0(), Helix._sxy.radius(), Helix._sxy.chi2DofCircle(),
	     phi0, dfdz, Helix._szphi.chi2DofLine());
    }

    // perform some clean up if needed

    chi2 = Helix._sxy.chi2DofCircle();
    if (chi2 >= _chi2xyMax) {
      if (_debug > 5) {
	printf("[CalHelixFinderAlg::%s] chi2 = %8.3f already at limit! no point can be added\n",banner,chi2);
      }
      goto F_END;
    }
//-----------------------------------------------------------------------------
// now add points
//-----------------------------------------------------------------------------
    if (_debug > 5) {
      printf("[CalHelixFinderAlg::%s] x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	     banner,Helix._sxy.x0(), Helix._sxy.y0(), Helix._sxy.radius(), Helix._sxy.chi2DofCircle());
      if (_debug > 10) {
	printf("[CalHelixFinderAlg::%s:LOOP] i       X        Y        dx        dy         costh        sinth2         e2     radial-dist\n",
	       banner);
      }
    }
  NEXT_ITERATION:;
    if (UsePhiResiduals == 1) {
      chi2_min    = _hitChi2Max;   // 2.*_hitChi2Max;
    } else{
      chi2_min    = _hitChi2Max;
    }
    dphiChi2    = 0.;

    ibest.Panel         = -1;
    ibest.PanelHitIndex = -1;
    
    wtBest      = -1;
    phiwtBest   = -1;

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhits          = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhits == 0)                                       continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;

      for (int i=seedPanelIndex; i<nhits; ++i){   
	hit = &panelz->fHitData.at(i);
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] >= 1)                    continue;

	hitPos    = hit->_pos;
	strawDir  = hit->_sdir;

	dr = calculateRadialDist(hitPos,helCenter,r);
	wt = calculateWeight    (*hit,helCenter,r);

	drChi2  = (dr*dr)*wt;
	  
	if ((UsePhiResiduals == 1) && (_phiCorrectedDefined)) {
	  phi_pred = panelz->z*dfdz + phi0;//hitPos.z()*dfdz + phi0;
	  dphi     = phi_pred - hit->_hphi;//_phiCorrected[i];
	  phiwt    = calculatePhiWeight(*hit,/*hitPos, strawDir,*/ helCenter, r, 0, banner);
	  dphiChi2 = dphi*dphi*phiwt;
	  // calculate distance from predicted point
	  x         = hitPos.x();
	  y         = hitPos.y();
	  x_pred    = helCenter.x() + r*cos(phi_pred);
	  y_pred    = helCenter.y() + r*sin(phi_pred);
	  hel_pred.SetX(x_pred);
	  hel_pred.SetY(y_pred);
	  weight_hel = calculateWeight(*hit, hel_pred, r)*_weight3D;
	  x_resid2  = (x - x_pred)*(x - x_pred);
	  y_resid2  = (y - y_pred)*(y - y_pred);
	  hitChi2   = (x_resid2 + y_resid2)*weight_hel;
	} 
	else {
	  hitChi2  = drChi2;
	}

	if (_debug > 10) {
	  printf("[CalHelixFinderAlg::%s:LOOP] sigmaphi2 = %10.3f drChi2 = %5.3f dphiChi2 = %5.3f chi2 = %5.3f wt = %8.3f\n",
		 banner, 1./phiwt, drChi2, dphiChi2, hitChi2,wt);
	}
	//-----------------------------------------------------------------------------
	// require chi2 < _hitChi2Max, identify the closest point
	//-----------------------------------------------------------------------------
	if ( (hitChi2 < chi2_min) && (drChi2 < _hitChi2Max) && (dphiChi2 < _hitChi2Max)) {
	  //-----------------------------------------------------------------------------
	  // check if XY-chi2 and ZPhi-chi2 are less than chi2xyMax and chi2zphiMax
	  //-----------------------------------------------------------------------------
	  x = hitPos.x();
	  y = hitPos.y();
	  Helix._sxy.addPoint  (x, y, wt);
	  Helix._szphi.addPoint(panelz->z, hit->_hphi, phiwt);

	  //update the number of strawHits
	  int     nHitSh = hit->nStrawHits();
	  Helix._nXYSh   += nHitSh;
	  Helix._nZPhiSh += nHitSh;

	  //update hit info
	  hit->_xyWeight   = wt;
	  hit->_zphiWeight = phiwt;

	  if (Helix._sxy.chi2DofCircle() < _chi2xyMax){
	    if (UsePhiResiduals == 1){
	      if (Helix._szphi.chi2DofLine() < _chi2zphiMax){
		chi2_min    = hitChi2;

		ibest.Panel         = p;
		ibest.PanelHitIndex = i;

		wtBest      = wt;
		phiwtBest   = phiwt;
	      }
	    }else {
	      chi2_min    = hitChi2;

	      ibest.Panel         = p;
	      ibest.PanelHitIndex = i;

	      wtBest      = wt;
	    }
	  }

	  Helix._sxy.removePoint  (x, y, wt);
	  Helix._szphi.removePoint(panelz->z, hit->_hphi, phiwt);

	}
      }
    }

    if (ibest.Panel >= 0){
      
      //evaluate the number of hits already in use in the same panel
      int    nHitsUsed(0), idSamePanel(-1);//, nLayers(2);
      for (int l=0; l<CalHelixFinderData::kNMaxHitsPerPanel; ++l){
	int   id = ibest.Panel*CalHelixFinderData::kNMaxHitsPerPanel + l;
	if (Helix._hitsUsed[id] == 1) {
	  ++nHitsUsed;
	  idSamePanel = l;
	}
      }
      
      panelz = &Helix._oTracker[ibest.Panel];
      hit    = &panelz->fHitData.at(ibest.PanelHitIndex);
      // int    layer_id = hit->strawId().layer();//_straw->id().getLayer();
      int    index    = ibest.Panel*CalHelixFinderData::kNMaxHitsPerPanel + ibest.PanelHitIndex;
      
      // now check two condition we want to skip:
      //  1) the hit we want to add belongs to a panel where we already reached the maximum number of hits allowed
      //  2) the hit belongs to the same layer where  we already reached the maximum number of hits allowed
      if (nHitsUsed >= _nHitsMaxPerPanel) {
	Helix._hitsUsed[index] = 10;
                                     goto NEXT_ITERATION;
      }
      if ( (_nHitsMaxPerPanel == 1) && (idSamePanel>=0) ){
	// int    pre_layer_id = panelz->fHitData.at(idSamePanel).strawId().layer();//_straw->id().getLayer();
	// if (pre_layer_id  == layer_id) {
	Helix._hitsUsed[index] = 10;
                                     goto NEXT_ITERATION;
	// }
      }

      x      = hit->x();
      y      = hit->y();
                                       //add point from the track
      Helix._sxy.addPoint(x, y, wtBest);
      int    nHitSh = hit->nStrawHits();
      Helix._nXYSh += nHitSh;

      if (UsePhiResiduals == 1){
      	Helix._szphi.addPoint(panelz->z, hit->_hphi, phiwtBest);
	Helix._nZPhiSh += nHitSh;
	dfdz  = Helix._szphi.dfdz(); 
	phi0  = Helix._szphi.phi0();
      }

      if (_debug > 5) {
	printf("[CalHelixFinderAlg::%s:PT2] x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %6.3f chi2Maxxy = %6.3f index point added = %i straw-id = %6i hitChi2 = %6.3f x = %8.3f y = %8.3f z = %9.3f\n",
	       banner,
	       Helix._sxy.x0(), Helix._sxy.y0(), Helix._sxy.radius(), Helix._sxy.chi2DofCircle(), _chi2xyMax, ibest.Panel,
	       hit->strawId().straw()/*_strawhit->strawIndex().asInt()*/, chi2_min,
	       x, y, panelz->z);//FIXME!
      }
					// mark point as active
      Helix._hitsUsed[index] = 1;

      r  = Helix._sxy.radius();
      helCenter.SetX( Helix._sxy.x0());
      helCenter.SetY( Helix._sxy.y0());

                                        // update helix      
      doWeightedCircleFit (Helix, SeedIndex, helCenter,  r);
      
      ++n_added_points;
	                              goto NEXT_ITERATION;
    }
//----------------------------------------------------------------------
//now update information about the radial residual of the hits
//----------------------------------------------------------------------
    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhits          = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhits == 0)                                                                    continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;

      for (int i=seedPanelIndex; i<nhits; ++i){   
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1)               continue;

	hit     = &panelz->fHitData.at(i);
	wt      = hit->_xyWeight;
	e2      = 1./wt;
	hitPos  = hit->_pos;
	dr      = calculateRadialDist(hitPos,helCenter,r);
	hitChi2 = dr*dr/e2;
	// store residual
	if ( (SeedIndex.Panel == 0) && (SeedIndex.PanelHitIndex == 0)){
	  hit->_drFromPred = fabs(dr);
	  double dzFromSeed = panelz->z - seedPanelz->z;              // expected to be positive (non-negative)
	  hit->_dzFromSeed  = dzFromSeed;
	}
      }
    }
//-----------------------------------------------------------------------------
// update circle parameters
//-----------------------------------------------------------------------------
    Helix._center.SetX(Helix._sxy.x0());
    Helix._center.SetY(Helix._sxy.y0());
    Helix._radius  = Helix._sxy.radius();
    Helix._chi2    = Helix._sxy.chi2DofCircle();

  F_END:;
    if (_debug > 5 ) {
      printf("[CalHelixFinderAlg::%s:END] N(added) = %i chi2 = %5.5f\n",banner,n_added_points,Helix._chi2);
    }

    if ((SeedIndex.Panel == 0) && (SeedIndex.PanelHitIndex == 0)){
      Helix._diag.n_rescued_points_16 = n_added_points;
    }
  }



//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::findTrack(SeedInfo_t&         SeedIndex     ,
				    CalHelixFinderData& Helix         ,
				    int                 UseMPVDfDz    ) {
//-----------------------------------------------------------------------------
// points used in the pattern recongition:
// ----------------------------------------------------------------
// center : center of the helix, center.z = z(hit with SeedIndex)
// p1     : center of the stopping target
// p2     : point in the position SeedIndex on the vector Xyzp
// p3     : postion of the EMC cluster
//-----------------------------------------------------------------------------
    double         radius, phi0, dx, dy, phi, Chi2;
    XYZVec         center, shPos, hePos;
//-----------------------------------------------------------------------------
// mode0GoodPoints: number of points belonging to a trajectory when dfdz is not
//                  re-calculated using the function calculateDfDz()
// mode1GoodPoints: number of points belonging to a trajectory when dfdz is
//                  re-computed using calculateDfDz()
//-----------------------------------------------------------------------------
    int            Mode(0), mode0GoodPoints(0), mode1GoodPoints(0), rescuedPoints(0);
//--------------------------------------------------------------------------------
// reset the temporary flags used to test the helix-search using the given triplet
//--------------------------------------------------------------------------------
                                                 //mark the seedHit as used
    PanelZ_t*      seedPanelz = &Helix._oTracker[SeedIndex.Panel];
    CalHelixPoint* seedHit    = &seedPanelz->fHitData.at(SeedIndex.PanelHitIndex);

    //mark the seed-hit as used
    int            index      = SeedIndex.Panel*CalHelixFinderData::kNMaxHitsPerPanel + SeedIndex.PanelHitIndex;
    Helix._hitsUsed[index] = 1;

    PanelZ_t*      panelz(0);

//---------------------------------------------------------------------
// define constrains on the z coordinate of the strawhit candidate for re-calculating dfdz
// If the candidate and the seeding strawhit are too close, the dfdz calculated could be 
// affected by several effects which lead to a wrong estimate
// We are asking that the candidate straw must be at a distance along
// the z axes greater than tollMin and less than tollMax.
// These parameters still need to be optimized
//-----------------------------------------------------------------------------
    double tollMin(100.) ; // , tollMax(500.);

					// parameters used to calculate the strawhit position residuals
    double weight(1.), wtarget(0.1);
    double deltaZ(0.), dist(0.), dist2(0.);
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
    double dfdz_end(-1e10), phi0_end(-1e10), radius_end(-1e10);

					// two flags are needed:
    bool removeTarget(true);            // avoid the recalculation of dfdz
					// and helix parameters in case when
                                        // others strawhit candidates are found
    double dfdz = _mpDfDz;		// tanLambda/radius (set to most probable);
//----------------------------------------------------------------------
// calculate helix paramters using the center of the stopping target,
// the EMC cluster which seeded the CalTimePeak and the seeding strawhit.
// The z coordinate of the target center is set to 0 because in the formula
// inside calculateTrackParameters(...) the z coordinate is not used
//-----------------------------------------------------------------------------
    //CalHelixPoint* lastHit = seedHit;
    PanelZ_t*      lastPanel = &Helix._oTracker[SeedIndex.Panel];

    float          panelHitChi2 [_nHitsMaxPerPanel] = {1e10};
    int            panelHitIndex[_nHitsMaxPerPanel] = {-1};

    XYZVec p1(0.,0.,0.);	       // target, z(ST) = 5971. - 10200. is not used
    XYZVec p2(seedHit->_pos);          // seed hit
    XYZVec p3(fCaloX,fCaloY,fCaloZ);   // cluster
    
    calculateTrackParameters(p1,p2,p3,center,radius,phi0,dfdz);
    
    double     tollMax = 2.*M_PI/dfdz;
//------------------------------------------------------------------------------
// helix parameters, in particular, phi0, are defined at Z=p2.z()
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE 
//------------------------------------------------------------------------------
    if (UseMPVDfDz ==1 ) {
      dfdz    = _hdfdz;			// _mpDfDz;
      tollMax = 2.*M_PI/dfdz;
    }

    SeedInfo_t lastIndex(-1,-1);

    ::LsqSums4 sxy;
    ::LsqSums4 szphi;

    sxy.addPoint(p2.x(), p2.y(), 1.     );  // seed hit
    sxy.addPoint(p3.x(), p3.y(), 1.     );  // EMC cluster position
    sxy.addPoint(    0.,     0., wtarget);  // Target center in the transverse plane, with small weight

    std::string name("CalHelixFinderAlg::findTrack");

    int  NPoints = seedHit->nStrawHits();     // nhits, associated with the track, sxy has NPoints+2 or NPoints+1 points
    int  NComboHits(1);

    double    z_phi0 = p2.z();

    CalHelixPoint* hit(0);

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhits          = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhits == 0)                                                                        continue;
	
      if( _debug > 10){
	if( p==SeedIndex.Panel ) {
	  printf("[%s:LOOP]  findTrack() starts with helix parameters derived from these points \n",name.data());
	  printf("[%s:LOOP]   point  type      X         Y         Z       xyzp-index \n", name.data());
	  printf("[%s:LOOP] ----------------------------------------------------------\n", name.data());
	  printf("[%s:LOOP] seeding        %9.3f %9.3f %9.3f %5i\n",name.data(),p2.x(),p2.y(),p2.z(),SeedIndex.Panel );
	  printf("[%s:LOOP] candidate      %9.3f %9.3f %9.3f %5i\n",name.data(),p1.x(),p1.y(),p1.z(),lastIndex.Panel );
	  printf("[%s:LOOP] emc cluster    %9.3f %9.3f %9.3f %5i\n",name.data(),p3.x(),p3.y(),p3.z(),             -1);
	  printf("[%s:LOOP]----------------------------------------------------------------------------------------------------------------------------------------\n",name.data());
	  printf("[%s:LOOP]  P     Z        xi       yi       xp       yp    dXYpred  dXYseed   dZseed    X0       Y0        R        phi      dfdz    chi2    \n",name.data());
	  printf("[%s:LOOP]----------------------------------------------------------------------------------------------------------------------------------------\n",name.data());
	}
      }
	
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex + 1;  

      //-----------------------------------------------------------------------------
      // dfdz = tanLambda/radius; phi0 is the last found hit phi
      //-----------------------------------------------------------------------------
      deltaZ = panelz->z - lastPanel->z;
      phi    = phi0 + deltaZ*dfdz;
      //evaluate the helix prediction using the z coordinate of the panel
      hePos.SetXYZ(center.x()+radius*cos(phi),center.y()+radius*sin(phi),panelz->z);

      //check the Panel-phi wrt to the DS center
      // double  hePosPhi = TVector2::Phi_0_2pi(polyAtan2(hePos.y(), hePos.x()));
      double  hePosPhi = polyAtan2(hePos.y(), hePos.x());
      if (hePosPhi < 0) hePosPhi = hePosPhi + 2*M_PI;
      double  deltaPhi = hePosPhi - panelz->phi;
      if ( deltaPhi > M_PI ) deltaPhi -= 2*M_PI;
      if ( deltaPhi < -M_PI) deltaPhi += 2*M_PI;
      if ( fabs(deltaPhi) > _maxPanelToHelixDPhi)                    continue;

      int  goodPoint(-1);               // index of the strawhit candidate for dfdz and helix parameters recalculation
	
      // reset the chi2 values calulated in the previous panel
      // for (int l=0; l<nLayers; ++l){
      for (int g=0; g<_nHitsMaxPerPanel; ++g){
	panelHitChi2 [g] = 1e10;
	panelHitIndex[g] = -1;
      }
      // }
      if (nhits > CalHelixFinderData::kNMaxHitsPerPanel) printf("[CalHelixFinderAlg::findTrack]ERROR!! more than 10 hits within the same panel!!!");
      for (int i=seedPanelIndex; i<nhits; ++i){   
	hit = &panelz->fHitData.at(i);

	hit->_dzFromSeed = 0;
	hit->_drFromPred = 0;

	shPos  = hit->_pos;

	// residuals in XY
	dx              = hePos.x() - hit->x();
	dy              = hePos.y() - hit->y();
	dist2           = dx*dx + dy*dy;
	dist            = std::sqrt(dist2);
	   
	if( _debug > 10){
	  // dist betw the straw hit and the seed hit in the transverse plane
	  double dx  = std::fabs(seedHit->x() - hit->x());
	  double dy  = std::fabs(seedHit->y() - hit->y());
	  double dxy = std::sqrt(dx*dx+dy*dy);
	  double chi2   = sxy.chi2DofCircle();

	  printf("[%s:LOOP] %3i %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f %8.5f %8.3f\n",
		 name.data(),p,hit->z(),hit->x(),hit->y(),hePos.x(),hePos.y(),dist,dxy,deltaZ,center.x(),center.y(),radius,phi,dfdz,chi2) ;
	}
	//-----------------------------------------------------------------------------
	// dxy_max: running search window accounts for the finite extrapolation accuracy
	//-----------------------------------------------------------------------------
	double dxy_max = _distPatRec + _dfdzErr*deltaZ;
	if (dist <= dxy_max) {
	  double delta_min(0);
	  int    index_min(-1);
	  // int    layer_id(hit->strawId().layer());//_straw->id().getLayer());
	  for (int k=0; k<_nHitsMaxPerPanel; ++k){
	    double delta = dist - panelHitChi2[k];
	    if (delta < delta_min){
	      delta_min = delta;
	      index_min = k;
	    }
	  }
	  if (index_min>=0) {
	    panelHitIndex[index_min] = i;
	    panelHitChi2 [index_min] = dist;
	  }
	}
      }

      double dist_min(1e10);
      //loop over the nHitsMaxPerPanel hits closest to the helix prediction
      // for (int l=0; l<nLayers; ++l){
      for (int t=0; t<_nHitsMaxPerPanel; ++t){
	int hitIndex = panelHitIndex[t];
	if (hitIndex < 0)                        continue;

	hit = &panelz->fHitData.at(hitIndex);

	//mark the hit as used
	index = p*CalHelixFinderData::kNMaxHitsPerPanel + hitIndex;
	Helix._hitsUsed[index] = 1;
	// NPoints++;
	++NComboHits;
	NPoints += hit->nStrawHits();

	//-----------------------------------------------------------------------------
	// Mode = 0: helix parameters evaluated using stopping_target+fitst_hit(seed)+cluster
	//           dphi/dz is fixed and so far set to the most proable value for CE
	//-----------------------------------------------------------------------------
	sxy.addPoint(hit->x(),hit->y(),weight);
	if (Mode == 1) {
	  // dfdz has already been evaluated, update XY part of the helix
	  center.SetX(sxy.x0());
	  center.SetY(sxy.y0());
	  radius = sxy.radius();
	}

	phi0      = polyAtan2(hit->y()-center.y(),hit->x()-center.x());  // *DOUBLE_CHECK*
	z_phi0    = panelz->z;//hitz;			             // *DOUBLE_CHECK*
	lastPanel = &Helix._oTracker[p];
	  
	if      ( Mode == 0 ) ++mode0GoodPoints;
	else if ((Mode == 1) && (panelHitIndex[t]<= lastIndex.PanelHitIndex)) ++mode1GoodPoints;//FIXME!

	double dzFromSeed = panelz->z - seedPanelz->z;              // expected to be positive (non-negative)
	hit->_dzFromSeed = dzFromSeed;
	hit->_drFromPred = panelHitChi2[t];  

	if ((panelHitChi2[t] < dist_min) && (dzFromSeed > tollMin) && (dzFromSeed < tollMax)) {       // FIXME! should we use the panel-z?
	  //-----------------------------------------------------------------------------
	  // goodPoint - index of the first hit separated in Z from the seed by > 10 cm
	  // after such found, the target center is removed and the circle parameters 
	  // recalculated using the cluster, the seed hit and the 'goodPoint' hit
	  // an additional requirement is that at the recalculation time there are 3 or more 
	  // hits found in total 
	  //-----------------------------------------------------------------------------
	  if (removeTarget) {
	    goodPoint = panelHitIndex[t];//i;
	    dist_min  = panelHitChi2 [t];
	  }
	}
      }
      // }

      if ((goodPoint >= 0) && (NComboHits >= 2) && (UseMPVDfDz == 0)) {
	//-----------------------------------------------------------------------------
	// the first point separated from the seed one by more than 10 cm has been found
	// recalculate helix parameters: for XY part use accumulated sxy sums
	// replace stopping target with the hit
	//-----------------------------------------------------------------------------
	//	sxy.removePoint(0.,0.,wtarget);
	hit = &panelz->fHitData.at(goodPoint);
	p1  = hit->_pos;

	center.SetX(sxy.x0());
	center.SetY(sxy.y0());
	radius = sxy.radius();
	//-----------------------------------------------------------------------------
	// now calculate more accuratelly the value of dfdz using just the two strawhit positions
	// change in the circle parameterization changes the phi0 value
	//-----------------------------------------------------------------------------
	phi0 = polyAtan2(hit->y()-center.y(),hit->x()-center.x());

	if (UseMPVDfDz == 0) {
	  // calculateDphiDz_2(Helix,SeedIndex,NPoints,center.x(),center.y(),dfdz);
	  calculateDphiDz_2(Helix,SeedIndex,NComboHits,center.x(),center.y(),dfdz);
	}
	else if (UseMPVDfDz ==1) {
	  dfdz = _hdfdz;
	}

	if (_debug > 10) {
	  printf("[%s:DEF2]  P      Z       X        Y   type\n", name.data());
	  printf("[%s:DEF2] ----------------------------------------------------\n", name.data());
	  printf("[%s:DEF2] %3i %9.3f %8.3f %8.3f seed \n", name.data(),SeedIndex.Panel,seedHit->z(),seedHit->x(),seedHit->y());//FIXME!
	  printf("[%s:DEF2] %3i %9.3f %8.3f %8.3f seed \n", name.data(), p,hit->z()    ,hit->x()    ,hit->y()    );
	}
	//-----------------------------------------------------------------------------
	// what to do if dfdz is negative? - the case of negative helicity is not covered yet
	//-----------------------------------------------------------------------------
	if ((dfdz > _maxDfDz) || (dfdz < _minDfDz)) {
	  //-----------------------------------------------------------------------------
	  // 2014-11-05 dPhi/Dz doesn't make sense, back to ground zero: 
	  // gianipez set dfdz equal to the most probable value for CE, 
	  //-----------------------------------------------------------------------------
	  if (_debug > 10) printf("[%s:DEF3] dfdz = %8.5f outside the limits. Continue the search\n",name.data(),dfdz);
	  p1.SetXYZ(0.,0.,0.);
	  dfdz = _mpDfDz;
	}
	else {
	  //-----------------------------------------------------------------------------
	  // dPhi/Dz makes sense, exclude the stopping target and change the search mode
	  //-----------------------------------------------------------------------------
	  removeTarget = false;
	  Mode         = 1;

	  lastIndex.Panel         = p; 
	  lastIndex.PanelHitIndex = goodPoint;
	  //	  goodPoint    = -1;
	}
      }
	  
    }

    if (NPoints < _minNHits) return;
    //-----------------------------------------------------------------------------
    // 3 or more points have been found on a helix candidate, update Chi2
    //-----------------------------------------------------------------------------
    Chi2 = sxy.chi2DofCircle();
    //-----------------------------------------------------------------------------
    // temporary variables to store dfdz values out of the method 'calculateDfDz(...)'
    //-----------------------------------------------------------------------------
    double dfdzRes  [3] = {   -1.,    -1.,    -1.};
    double dphi0Res [3] = {-9999., -9999., -9999.};
    double radiusRes[2] = {   -1.,    -1.};

    if (_diag > 0) {
      if (UseMPVDfDz == 0) dfdzRes[0] = dfdz;
      dphi0Res [0] = phi0 - dfdz*z_phi0;
      radiusRes[0] = sxy.radius();
    }
//-----------------------------------------------------------------------------
// initialize only the xy part, z-phi part is not needed here
//-----------------------------------------------------------------------------
    Helix._sxy.init(sxy);
    Helix._nXYSh      = NPoints;
    Helix._radius     = sxy.radius();
    Helix._center.SetXYZ(sxy.x0(), sxy.y0(), 0.0);
    Helix._nStrawHits = NPoints;
    Helix._nComboHits = NComboHits;
    Helix._dfdz       = dfdz;
    Helix._fz0        = phi0 - dfdz*z_phi0; // *DOUBLE_CHECK*
 
    radius_end = Helix._radius;

    int rc = refineHelixParameters(Helix, SeedIndex);
//-----------------------------------------------------------------------------
// if weighted XY fit didn't converge, there is nothing else one can do, return
//-----------------------------------------------------------------------------
    if (rc < 0) return;

    // Helix._center.SetXYZ(Helix._cw.x(), Helix._cw.y(), 0.0);
    // Helix._radius  = Helix._rw;
    radius_end     = Helix._radius;

    // Helix._sxy.init(Helix._sxyw); 
					// doWeightedCircleFit still adds the ST and the cluster
    Chi2    = Helix._sxy.chi2DofCircle();
    NPoints = Helix._nStrawHits;   //  *FIXME*  in principle, the fit can remove ST as well as the cluster
					// diagnostics
    radiusRes[1] = Helix._radius;
//-----------------------------------------------------------------------------
// 2015-01-22 G. Pezzullo and P. Murat; update the dfdz value using all hits
//-----------------------------------------------------------------------------
    int rs = findDfDz(Helix, SeedIndex);
    if (rs ==1 ) {
      Helix._dfdz = _hdfdz;
      Helix._fz0  = _hphi0;
					// fill diag vector
      dfdzRes[1]  = _hdfdz;
      dphi0Res[1] = _hphi0;
    }
//-----------------------------------------------------------------------------
// 2015-01-23 G. Pezzu and P. Murat: when it fails, doLinearFitPhiZ returns negative value
//                                   in this case, use the previous value for dfdz and phi0
//-----------------------------------------------------------------------------
    bool rcPhiZ = doLinearFitPhiZ(Helix, SeedIndex);

    if (rcPhiZ) {
      dfdz_end   = Helix._dfdz;
      phi0_end   = Helix._fz0;
					// fill diagnostic vector
      dfdzRes [2] = Helix._dfdz;
      dphi0Res[2] = Helix._fz0;

      NPoints    = 0;
      NComboHits = 0;

      //FIXME! implement a function
      for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
	panelz = &Helix._oTracker[p];
	int  nhitsPerPanel  = panelz->fNHits;
	int  seedPanelIndex(0);
	if (nhitsPerPanel == 0)                                 continue;
	if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

	for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	  hit   = &panelz->fHitData.at(i);
	  index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	  if (Helix._hitsUsed[index] > 0 )  {
	    ++NComboHits;
	    NPoints += hit->nStrawHits();
	  }
	}
      }
    }
    else {
      dfdz_end = _hdfdz;
      phi0_end = _hphi0;
    }

    if (_debug > 10) {
      printf("[%s] strawhit type     X        Y        Z     index     \n", name.data());
      printf("[%s] ----------------------------------------------------\n", name.data());
      printf("[%s]    seeding   %9.3f %9.3f %9.3f   %i  \n", name.data(),p2.x(), p2.y(), p2.z(), SeedIndex.Panel);//FIXME!
      printf("[%s]   candidate  %9.3f %9.3f %9.3f   %i  \n", name.data(),p1.x(), p1.y(), p1.z(), lastIndex.Panel);//FIXME!
      printf("[%s]  emc cluster %9.3f %9.3f %9.3f       \n", name.data(),p3.x(), p3.y(), p3.z());
      printf("[%s] NPoints = %i x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.6fdfdz = %5.6f chi2 = %5.3f \n", 
	     name.data(),NPoints,Helix._sxy.x0(), Helix._sxy.y0(), radius_end, phi0_end, dfdz_end , Helix._sxy.chi2DofCircle());
    }

    if (mode1GoodPoints > 0) rescuedPoints = mode1GoodPoints - mode0GoodPoints ;
    else                     rescuedPoints = -1;
//-----------------------------------------------------------------------------
// all points in the time window checked
//----------------------------------------------------------------------
    if (NComboHits > 2) Chi2 = Helix._sxy.chi2DofCircle();
    else                Chi2 = -1;

    if (_debug > 5) {
      printf("[%s:END] SeedIndex-Panel:%3i Mode= %i UseMPVDfDz:%1i NPoints= %3i Chi2=%5.3f dfdz= %5.5f \n",
	     name.data(),SeedIndex.Panel,Mode,UseMPVDfDz,NPoints,Chi2,dfdz_end);//FIXME!
    }

    // can execution really come here with Mode == 0 ? - YES! not sure, why
    //      if (( NPoints >  _goodPointsTrkCandidate) ||
    // 	  ((NPoints == _goodPointsTrkCandidate) && (Chi2< _chi2TrkCandidate))) {
    // //-----------------------------------------------------------------------------
    // found candidate is better, than the best previous one
    //-----------------------------------------------------------------------------
    // _goodPointsTrkCandidate = NPoints;
    // _chi2TrkCandidate       = Chi2;
    //-----------------------------------------------------------------------------
    // reset the vector holding the informations about:
    // -> hit belonging to the track candidate
    // -> distance in the X-Y plane from the prediction
    // -> distance from the seeding hit along the z-axes
    //-----------------------------------------------------------------------------

    Helix._nStrawHits = NPoints;
    Helix._nComboHits = NComboHits;

    Helix._helixChi2 = Chi2;

    Helix._seedIndex = SeedIndex;
    Helix._candIndex = lastIndex;
    //-----------------------------------------------------------------------------
    // helix parameters and LSQ sums, phi is defined at z=0
    //-----------------------------------------------------------------------------
    Helix._center.SetXYZ(Helix._sxy.x0(),Helix._sxy.y0(), 0.0);
    Helix._radius = radius_end;
    Helix._fz0    = phi0_end;
    Helix._dfdz   = dfdz_end;
	
    if (_diag > 0){
      Helix._diag.loopId_4           = _findTrackLoopIndex;
      Helix._diag.radius_5           = Helix._radius;
      Helix._diag.n_rescued_points_9 = rescuedPoints;

      double dz = p1.z() - p2.z();

      Helix._diag.dz_10              = (Mode == 1) ? dz : -1.;
      Helix._diag.n_active_11        = NPoints;
      Helix._diag.chi2_dof_circle_12 = sxy.chi2DofCircle();
      Helix._diag.chi2_dof_line_13   = szphi.chi2DofLine();

      Helix._diag.dfdzres_17 = dfdzRes[0];
      Helix._diag.dfdzres_18 = dfdzRes[1];
      Helix._diag.dfdzres_19 = dfdzRes[2];

      Helix._diag.dr_20 = radiusRes[0];
      Helix._diag.dr_21 = radiusRes[1];

      Helix._diag.dphi0res_22 = dphi0Res[0];
      Helix._diag.dphi0res_23 = dphi0Res[1];
      Helix._diag.dphi0res_24 = dphi0Res[2];

      int j=0;
      for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
	panelz = &Helix._oTracker[p];
	int  nhitsPerPanel  = panelz->fNHits;
	int  seedPanelIndex(0);
	if (nhitsPerPanel == 0)                                                            continue;
	if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

	for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	  hit = &panelz->fHitData.at(i);
	  // if (hit->_used != 1)              continue;
	  int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	  if (Helix._hitsUsed[index] != 1)               continue; 

	  if (j < Helix.maxIndex()) {
	    Helix._diag.dist[j] = hit->_drFromPred;
	    Helix._diag.dz  [j] = hit->_dzFromSeed;
	    ++j;
	  }
	  else {
	    printf("ERROR in CalHelixFinderAlg::findTrack : index out limits. IGNORE; \n");
	  }
	}
      }
    }

  }

  //-----------------------------------------------------------------------------
  // helix parameters are defined at Z=p2.z, Phi0 corresponds to p2
  //-----------------------------------------------------------------------------
  void CalHelixFinderAlg::calculateTrackParameters(const XYZVec&   p1       ,
						   const XYZVec&   p2       ,
						   const XYZVec&   p3       ,
						   XYZVec&         Center   ,
						   double&             Radius   ,
						   double&             Phi0     ,
						   double&             DfDz23) {
    Center.SetZ(p2.z());

    double x_m, y_m, x_n, y_n;
    // coordinates of the mean point between p1 and p3
    x_m = (p3.x() + p1.x())/2.;
    y_m = (p3.y() + p1.y())/2.;
    // coordinates of the mean point between p2 and p3
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

    //check we are not in a degenerate case where m is close to k, which rapresents two almost parallel lines
    double limit = 0.8;//FIXME!
    
    if ( (m/k>0) && ( (m/k) - int(m/k) > limit) ) {//invert p3 with p1 and recalculate: x_n, y_n, k, t
      x_n = (p1.x() + p2.x())/2.;
      y_n = (p1.y() + p2.y())/2.;
      k   = -1.*(p1.x() - p2.x())/(p1.y() - p2.y());
      t   = y_n - x_n*k;
     }

    // calculate Center.x and Center.y
    double x0 = (t - c)/(m - k);//(c - t) * (k*m)/(m-k);
    Center.SetX(x0);
    double y0 = m*x0 + c;   //(c - t) * m / (m - k) + t;
    Center.SetY(y0);
//-----------------------------------------------------------------------------
// calculate the radius,phi0, tanLambda assuming that the helix also crosses 
// the point (0,0). Note that the Z-position of the stopping target is not used
//-----------------------------------------------------------------------------
    double dx3  = p3.x() - x0;
    double dy3  = p3.y() - y0;
    double dz32 = p3.z() - p2.z();

    Radius      = std::sqrt(dx3*dx3+dy3*dy3);

    double dx2  = (p2.x() - x0);
    double dy2  = (p2.y() - y0);

    Phi0        = polyAtan2(dy2,dx2);
//-----------------------------------------------------------------------------
// this assumes that the helix is right-handed, *FIXME*
// make sure that we are lookign for a particle which makes the number of turns
// close to the expected 
//-----------------------------------------------------------------------------
    double dphi32 = polyAtan2(dy3,dx3) - Phi0;
    if (dphi32 < 0.) dphi32 += 2.*M_PI;

    //    double exp_dphi = _mpDfDz*dz32;

    DfDz23 = dphi32/dz32; 

    double   diff      = fabs(DfDz23 - _mpDfDz);
    double   diff_plus = fabs( (dphi32 + 2.*M_PI)/dz32 -_mpDfDz );
    while ( diff_plus < diff ){
      dphi32  = dphi32 + 2.*M_PI;
      DfDz23      = dphi32/dz32;
      diff      = fabs(DfDz23 - _mpDfDz);
      diff_plus = fabs( (dphi32 + 2.*M_PI)/dz32 -_mpDfDz );
    }
    
    double   diff_minus = fabs( (dphi32 - 2.*M_PI)/dz32 -_mpDfDz );
    while ( diff_minus < diff ){
      dphi32   = dphi32 - 2.*M_PI;
      DfDz23       = dphi32/dz32;
      diff       = fabs(DfDz23 - _mpDfDz);
      diff_minus = fabs( (dphi32 - 2.*M_PI)/dz32 -_mpDfDz );
    }

    //check id DfDz is within the range 
    if ( (DfDz23 < _minDfDz) || (DfDz23 > _maxDfDz)) DfDz23 = _mpDfDz;

    if (_debug > 5) {
//-----------------------------------------------------------------------------
// in debug mode also want to print the helix parameters, calculate them
//-----------------------------------------------------------------------------
      double d0     = sqrt(x0*x0+y0*y0)-Radius;
      double phi00  = polyAtan2(y0,x0)+M_PI/2;   // for negatively charged particle
      double tandip = DfDz23*Radius;
      double dphi   = phi00-Phi0-M_PI/2;
      if (dphi < 0) dphi += 2*M_PI;           // *FIXME* right-handed ellipse

      double z0     = p2.z()-dphi*dz32/dphi32;

      printf("[CalHelixFinderAlg:calculateTrackParameters] X0: %9.3f Y0: %9.3f phi0: %8.5f p1.z = %9.3f p2.z = %9.3f p3.z = %9.3f dphi32 = %8.5f dfdz = %8.5f\n",
	     Center.x(),Center.y(),Phi0,p1.z(),p2.z(),p3.z(),dphi32,DfDz23);
      printf("[CalHelixFinderAlg:calculateTrackParameters] z0 = %9.3f d0 = %8.4f  phi00 = %8.5f omega = %8.5f tandip = %8.4f\n",z0,d0,phi00,1/Radius,tandip);
    }
  }

//-----------------------------------------------------------------------------
// the function is currently not called
//-----------------------------------------------------------------------------
  void  CalHelixFinderAlg::calculateDfDz(double phi0, double phi1, double z0, double z1, double& DfDz) {
    double   deltaPhi  = TVector2::Phi_mpi_pi(phi1-phi0);
    DfDz               = deltaPhi/(z1-z0);

    // 2018-01-02: don't do that!
    // double   diff      = fabs(DfDz - _mpDfDz);
    // double   diff_plus = fabs((deltaPhi + 2.*M_PI)/(z1-z0) -_mpDfDz);
    // while (diff_plus < diff) {
    //   deltaPhi  = deltaPhi + 2.*M_PI;
    //   DfDz      = deltaPhi/(z1-z0);
    //   diff      = fabs(DfDz - _mpDfDz);
    //   diff_plus = fabs( (deltaPhi + 2.*M_PI)/(z1-z0) -_mpDfDz );
    // }
    
    // double   diff_minus = fabs((deltaPhi - 2.*M_PI)/(z1-z0) -_mpDfDz);
    // while (diff_minus < diff) {
    //   deltaPhi   = deltaPhi - 2.*M_PI;
    //   DfDz       = deltaPhi/(z1-z0);
    //   diff       = fabs(DfDz - _mpDfDz);
    //   diff_minus = fabs( (deltaPhi - 2.*M_PI)/(z1-z0) -_mpDfDz );
    // }
  }

//-----------------------------------------------------------------------------
// assume particle makes less than a full turn
//-----------------------------------------------------------------------------
  void  CalHelixFinderAlg::calculateDphiDz_2(CalHelixFinderData& Helix, SeedInfo_t SeedIndex, 
					     int NHits, double X0, double Y0, double& DphiDz) {
    LsqSums2       sphiz;

    double         phi(0), phi0(0);//, phiCl(0);
    
    CalHelixPoint* hit(0);
    PanelZ_t*      panelz(0);
    XYZVec*        pos(0);
    
    int            counter(0);

    bool           isFirst(true);

    for (int p=SeedIndex.Panel; p<CalHelixFinderData::kNTotalPanels; ++p){
      if (counter > NHits  )                               break;
      panelz = &Helix._oTracker[p];
      int  nhitsPerPanel  = panelz->fNHits;
      int  seedPanelIndex(0);
      if (nhitsPerPanel == 0)                           continue;
      if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;  

      for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){   
	hit = &panelz->fHitData.at(i);
	int index = p*CalHelixFinderData::kNMaxHitsPerPanel + i;
	if (Helix._hitsUsed[index] != 1)                     continue;

	pos = &(hit->_pos);
	phi = polyAtan2(pos->y()-Y0,pos->x()-X0);
	if (isFirst) { phi0 = phi; isFirst = false;}

	if (phi-phi0 >  M_PI) phi -= 2*M_PI;
	if (phi-phi0 < -M_PI) phi += 2*M_PI;

	sphiz.addPoint(pos->z(),phi);
	++counter;

	if (_debug > 10) printf("[CalHelixFinderAlg::calculateDphiDz_2:LOOP] panel,id,phi,z=%3i %3i %8.5f %9.3f\n",p,i,pos->z(),phi);
      }
    }
//-----------------------------------------------------------------------------
// define straight line phi = phi0+dPhi/Dz*z , where phi0 = phi(z=0)
//-----------------------------------------------------------------------------
    phi0   = sphiz.yMean();
    DphiDz = sphiz.dydx();

    if (_debug > 5) printf("[CalHelixFinderAlg::calculateDphiDz_2:END] phi0,DphiDz = %9.5f %9.5f \n",phi0,DphiDz);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void CalHelixFinderAlg::plotXY(int ISet) {

  std::vector<mu2e::CalHelixPoint>* xyzp;
  CalHelixFinderData*          helx;

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

  mu2e::CalHelixPoint*    hit;
  XYZVec* pos;

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
	   i,(int) hit->index(),
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
      x0  = helx->_sxy.x0();
      y0  = helx->_sxy.y0();
      r   = helx->_sxy.radius();

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

  helx->print("from CalHelixFinderAlg::plotXY");

  printf("All Done\n");
}
//-----------------------------------------------------------------------------
// this routine is supposed to be called interactively from the ROOT prompt
// it has to retrieve a pointer to CalHelixFinderAlg called from CalPatRec
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::plotZPhi(int ISet) {

    TCanvas* c;
    int      color;
    TMarker* m;
    char     name[200];
//-----------------------------------------------------------------------------
// retrieve the data points for storing in TGraphs
//-----------------------------------------------------------------------------
    std::vector<mu2e::CalHelixPoint>* xyzp;
    CalHelixFinderData*          helx;

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

    mu2e::CalHelixPoint* hit;

    XYZVec* pos;

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
      phi [i] = polyAtan2(hit->_pos.y(), hit->_pos.x());//hit->_pos.phi();

      phi1[i] = polyAtan2(y[i]-helx->_center.y(),x[i]-helx->_center.x());

      printf("i: %3i ind: %5i x: %10.3f y: %10.3f z: %10.3f phi: %10.3f phi1: %10.3f 0x%08x %5i\n",
	     i,(int) hit->index(),
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
    
    helx->print("from CalHelixFinderAlg::plotZPhi");
  }

//-----------------------------------------------------------------------------
// 'TrkIndex' = 0 : mark all hits as active
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::saveResults(CalHelixFinderData& Helix, int Index) {
    //    _results[Index]._xyzp  = Xyzp;
    _results[Index]._helix = Helix;
  }


//-----------------------------------------------------------------------------
// 2017-01-25 P.Murat: print
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::printXYZP(CalHelixFinderData& Helix) {
    int n = Helix._nFiltPoints;//_xyzp.size();

    printf("[CalHelixFinderAlg::printXYZP]-----------------------------------------------------------------------------------------\n");
    printf("[CalHelixFinderAlg::printXYZP]     i index strawID   flag    Used     X         Y         Z      _debug: %5i nhits: %5i\n",_debug,n);
    printf("[CalHelixFinderAlg::printXYZP]-----------------------------------------------------------------------------------------\n");
    
    const vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();

    PanelZ_t* panelz(0);
    
    for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhitsPerPanel  = panelz->fNHits;


      for (int i=0; i<nhitsPerPanel; ++i){   
	CalHelixPoint* pt = &panelz->fHitData.at(i);

	int loc    = shIndices[i];	 // index in chcol of i-th timecluster hit
    
	// const StrawHit& sh          = Helix.chcol()->at(loc);

	printf("[CalHelixFinderAlg::printXYZP] %5i %5i %5i   %08x   %2i %9.3f %9.3f %9.3f \n",
	       i, loc,  pt->strawId().straw()/*sh.strawIndex().asInt()*/, *((int*) &pt->_flag), pt->use(), pt->_pos.x(), pt->_pos.y(), pt->_pos.z());
      }
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::resolve2PiAmbiguity(CalHelixFinderData& Helix,const XYZVec& Center, double DfDz, double Phi0){
  
    const XYZVec*  pos;
    double         z, phi, phi_ref, dphi, dx, dy;

    PanelZ_t*      panelz(0);
    CalHelixPoint* hit(0);

    for (int p=0; p<CalHelixFinderData::kNTotalPanels; ++p){
      panelz = &Helix._oTracker[p];
      int  nhits          = panelz->fNHits;
      for (int i=0; i<nhits; ++i){   
	hit = &panelz->fHitData.at(i);
	pos = &hit->_pos;
	z   = pos->z();

	dx  = (pos->x() - Center.x());
	dy  = (pos->y() - Center.y());
	phi = polyAtan2(dy, dx);
	if (phi < 0) phi = phi + 2*M_PI;//phi = TVector2::Phi_0_2pi(phi);
	// predicted value of phi
	phi_ref = z*DfDz + Phi0;
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
	// _phiCorrected[i] = phi;
      }
			
    }
		// don't know
    _phiCorrectedDefined = 0;
  }



//---------------------------------------------------------------------------
// reset track paramters
// indices on the xyzp vector of: the straw hit seeding the search,
// the second strawhit used for recalculating the dfdz value
//---------------------------------------------------------------------------
  void CalHelixFinderAlg::resetTrackParamters() {

    //    fSeedIndex      = -9999;
    //    fCandIndex      = -9999;
    //    fLastIndex      = -9999;
    fUseDefaultDfDz = 0;

    // _x0             = -9999.;
    // _y0             = -9999.;
    //_phi0           = -9999.;
    //    _radius         = -9999.;
    _hphi0                  = -9999.;
//-----------------------------------------------------------------------------
// quality paramters used for doing comparison between several track candidates
//-----------------------------------------------------------------------------
    // _goodPointsTrkCandidate = 0;
    // _chi2TrkCandidate       = 1e10;
    _hdfdz                  = _mpDfDz;
  }

}


