///////////////////////////////////////////////////////////////////////////////
// helix fit to straw hits
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
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "art/Framework/Services/Optional/TFileService.h"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
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
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "boost_fix/accumulators/statistics.hpp"
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

#include "CalPatRec/inc/CalHelixFinderAlg.hh"

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

//--------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------
  // double     CalHelixFinderAlg::evalApproximatePhi(double& Y, double& x){


  // }


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
    pvec[HelixTraj::phi0Index] = atan2(-amsign*Helix._center.x(),amsign*Helix._center.y());
    // d0 describes the distance to the origin at closest approach.
    // It is signed by the particle angular momentum WRT the origin.
    // The Helix fit radial bias is anti-correlated with d0; correct for it here.
    pvec[HelixTraj::d0Index] = amsign*(Helix._center.perp() - Helix._radius); //  - 2*_rbias);
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
    //    if ((_goodPointsTrkCandidate < _minPointsTrkCandidate) || (_markCandidateHits == 0)) return 0;
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
    _smartTag         (pset.get<int>   ("smartTag"         )),
    _hsel             (pset.get<vector<string> >("HelixFitSelectionBits"  )),
    _bkgsel           (pset.get<vector<string> >("BackgroundSelectionBits")),
    _maxElectronHitEnergy(pset.get<double>("maxElectronHitEnergy")),
    _minNHits         (pset.get<int>   ("minNHit"          )),
    //    _maxDz            (pset.get<double>("maxdz",35.0)),
    _mpDfDz           (pset.get<double>("mostProbableDfDz")),
    _minNSt           (pset.get<double>("minNActiveStationPairs")),
    _dzOverHelPitchCut(pset.get<double>("dzOverHelPitchCut")),
    _maxDfDz          (pset.get<double>("maxDfDz",0.01)),
    _minDfDz          (pset.get<double>("minDfDz",5e-04)),
    _sigmaPhi         (pset.get<double>("sigmaPhi")),
    _weightXY         (pset.get<double>("weightXY")),
    _weightZPhi       (pset.get<double>("weightZPhi")),
    _weight3D         (pset.get<double>("weight3D")),
    _ew               (pset.get<double>("errorAlongWire")),
    _maxXDPhi         (pset.get<double>("maxXDPhi",5.)),
    _distPatRec       (pset.get<double>("distPatRec")),
    _rhomin           (pset.get<double>("rhomin",350.0)),
    _rhomax           (pset.get<double>("rhomax",780.0)),
    _mindist          (pset.get<double>("mindist",50.0)),
    _maxdist          (pset.get<double>("maxdist",500.0)),
    _pmin             (pset.get<double>("minP",50.0)),
    _pmax             (pset.get<double>("maxP",150.0)),
    _tdmin            (pset.get<double>("minAbsTanDip",0.3)),
    _tdmax            (pset.get<double>("maxAbsTanDip",2.0)),
    _rcmin            (pset.get<double>("rcmin",200.0)),
    _rcmax            (pset.get<double>("rcmax",350.0)),
    _xyweights        (pset.get<bool>  ("xyWeights",false)),
    _zweights         (pset.get<bool>  ("zWeights",false)),
    _filter           (pset.get<bool>  ("filter",true)),
    _plotall          (pset.get<bool>  ("plotall",false)),
    _usetarget        (pset.get<bool>  ("usetarget",true)),
    _bz               (0.0),
    // _x0               (-9999.),
    // _y0               (-9999.),
    // _phi0             (-9999.),
    //    _radius           (-9999.),
    //    _dfdz             (-9999.),
    _goodPointsTrkCandidate(-9999),
    _minPointsTrkCandidate (pset.get<int>("minPointsTrkCandidate")),
    _chi2TrkCandidate      (1e10),
    _maxChi2TrkCandidate   (pset.get<double>("maxChi2TrkCandidate")),
    _hitChi2Max            (pset.get<double>("hitChi2Max"         )),
    _chi2xyMax             (pset.get<double>("chi2xyMax")),
    _chi2zphiMax           (pset.get<double>("chi2zphiMax")),
    _chi2hel3DMax          (pset.get<double>("chi2hel3DMax")),
    _dfdzErr               (pset.get<double>("dfdzErr")){

    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    CalHelixPoint::_useflag = StrawHitFlag(bitnames);

    for (int i=0; i<kMaxNHits; ++i){
      _indicesTrkCandidate[i] = -9999;
      _distTrkCandidate   [i] = -9999;
      _dzTrkCandidate     [i] = -9999;
    }

    _chi2nFindZ    = 0.0,
    _eventToLook   = -1;

    _hDfDzRes = new TH1F("hDfDzRes","dfdz residuals" , 20, _minDfDz, _maxDfDz);
    _hPhi0Res = new TH1F("hPhi0Res", "phi0 residuals", 20, 0., 2.*M_PI);

  }


//-----------------------------------------------------------------------------
  CalHelixFinderAlg::~CalHelixFinderAlg() {
    delete _hDfDzRes;
    delete _hPhi0Res;
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
      printf("[CalHelixFinderAlg::findHelix] Helix._sxy.qn() = %5.0f goodPointsTrkCandidate = %i\n",
	     Helix._sxy.qn(), _goodPointsTrkCandidate);
    }

    if (_goodPointsTrkCandidate < _minNHits ) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,1); // small number of hits
    }
    else if ((Helix._radius < _rmin) || (Helix._radius > _rmax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,2); // initialization failure
    }
    else if ((Helix._sxy.qn() < _minNHits) || (Helix._sxy.chi2DofCircle() > _chi2xyMax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,3); // xy reconstruction failure
    }
    else if ((Helix._srphi.qn() < _minNHits) || (Helix._srphi.chi2DofLine() > _chi2zphiMax)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,4); // phi-z reconstruction failure
    }
    else {
//-----------------------------------------------------------------------------
// success, form output
//-----------------------------------------------------------------------------
      Helix._goodhits.clear();

      int n = _xyzp.size();

      for (int i=0; i<n; ++i) {
	if (_xyzp[i].isOutlier()) continue;
	Helix._goodhits.push_back(_xyzp[i]._ind);
      }

      defineHelixParams(Helix);
      retval = true;
    }

    return retval;
  }

//----------------------------------------------------------------------------------------
// 2015-01-13  calculate track DphiDz using histogrammed distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  int CalHelixFinderAlg::findDfDz(CalHelixFinderData& Helix, int SeedIndex, int *IndexVec, int  Diag_flag) {

    double phi, phi_ref(-1e10), z, z_ref, dphi, dz, dzOverHelPitch;

    CLHEP::Hep3Vector* center = &Helix._center;
    CLHEP::Hep3Vector pos_ref;

    _hDfDzRes->Reset();
    _hPhi0Res->Reset();
					// 2015 - 03 -30 G. Pezzu changed the value of tollMax.
					// using the initial value of dfdz we can set it more accuratelly:
					// tollMax = half-helix-step = Pi / dfdz
    double tollMin(100.);
//-----------------------------------------------------------------------------
// 2017-09-26 gianipez fixed a bug: in case the Helix phi-z fit didn't converge yet, 
// Helix._dfdz is set to -1e6, so we need to make a check here!
// this is a tempOrary fix that doesn't take into account the particle helicity. FIX ME!
//-----------------------------------------------------------------------------
    double helix_dfdz(_mpDfDz);
    // 2017-11-14 gianipez: findDfDz shoudl use the dfdz value obtained only from the linearFit
    if (Helix._srphi.qn() >= 10) helix_dfdz = Helix._srphi.dfdz();
    //    if (Helix._dfdz > 0) helix_dfdz =  Helix._dfdz;
    double tollMax = 2.*M_PI / helix_dfdz; 

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::findDfDz:BEGIN] x0 = %9.3f y0 = %9.3f radius = %9.3f",
	     center->x(), center->y(), Helix._radius);
      printf("dfdz = %9.6f straw-hits = %9.5f dzPitch = %8.6f\n",
	     helix_dfdz, (Helix._sxy.qn() - 1), tollMax);
    }

    int       np, nstations, nhits[30];
    double    phiVec[30], zVec[30], weight(0), weight_cl(0);

    np        = _xyzp.size();
    nstations = _tracker->nStations();

    for (int i=0; i<nstations; i++) {
      phiVec[i] = 0;
      zVec  [i] = 0;
      nhits [i] = 0;
    }
//-----------------------------------------------------------------------------
// Part 1: use only contiguous parts of the trajectory
//-----------------------------------------------------------------------------
    for (int i=SeedIndex; i<np; i++) {
      if ((IndexVec[i] > 0) && (!_xyzp[i].isOutlier())) {
	int ist = _xyzp[i]._straw->id().getStation();
	CLHEP::Hep3Vector* pos = &_xyzp[i]._pos;
	phi = atan2(pos->y()-center->y(),pos->x()-center->x()); // atan2 returns its result in [-pi,pi], convert to [0,2pi]
	if (phi < 0) phi += 2*M_PI;
	zVec  [ist] += pos->z();
	if (nhits[ist] == 0) phiVec[ist] = phi;
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
      for (int i=0; i<nstations; i++) {
	if (nhits[i] > 0) printf("[CalHelixFinderAlg::findDfDz] %5i %6i    %9.3f %8.5f\n", i,nhits[i],zVec[i],phiVec[i]);
      }
    }

    int i0(-1), first_point(1);
					// add the cluster phi
    double zCl   = fCaloZ;
    double phiCl = atan2(fCaloY-center->y(),fCaloX-center->x());
    if (phiCl < 0) phiCl += 2*M_PI;

    for (int i=0; i<nstations; i++) {
      if (nhits[i] == 0)                                    continue; 
				        // find station corresponding to the first point
      if (first_point) {
	i0          = i;
	first_point = 0;
      }

      phi_ref = phiVec[i];
      z_ref   = zVec  [i];

      for(int j=i+1; j<nstations; ++j){
	if (nhits[j] == 0)                                  continue;
	phi = phiVec[j];
	z   = zVec  [j];
	dz  = z - z_ref;
	
	dzOverHelPitch = dz/tollMax - int(dz/tollMax);
	weight         = nhits[i] + nhits[j];

	if ((phi_ref > -9999) && (dzOverHelPitch < _dzOverHelPitchCut) && (dz > tollMin)){
	  dphi = phi-phi_ref;
 	  while (dphi >  M_PI) dphi -= 2*M_PI;
 	  while (dphi < -M_PI) dphi += 2*M_PI;
//-----------------------------------------------------------------------------
// add 2 pi for taking into account the fact we are in the second loop
// FIX ME: what to do if we are in the third loop?
//-----------------------------------------------------------------------------
	  if (dz > tollMax) dphi += 2*M_PI*int(dz/tollMax);

	  double dphidz = dphi/dz;
	  while (dphidz < 0.) {
	    dphi    = dphi+2.*M_PI;
	    dphidz  = dphi/dz;
	  }
	  _hDfDzRes->Fill(dphidz, weight);

	  double tmpphi0 = phi_ref - dphidz*z_ref;
	  tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

	  if (_debug >5) {
	    printf("[CalHelixFinderAlg::findDfDz] z_ref = %9.3f z = %9.3f phi_ref = %9.5f",z_ref,z,phi_ref);
	    printf(" dphi-phi_ref = %9.5f dz = %10.3f dz/HelPitch = %10.3f dphi/dz = %9.5f phi0 = %9.6f\n",
		   dphi-phi_ref, dz, dzOverHelPitch, dphi/dz, tmpphi0);
	  }
//-----------------------------------------------------------------------------
// in case dfdz is out of limits set tmpphi0 as negative
//-----------------------------------------------------------------------------
	  if ((dphidz < _minDfDz) || (dphidz >  _maxDfDz)) tmpphi0 = -1;
	  _hPhi0Res->Fill(tmpphi0, weight);
	}
      }
					// use the calorimeter cluster phi
      dz             = zCl - z_ref;
      dzOverHelPitch = dz/tollMax - int(dz/tollMax);
      weight_cl      =  nhits[i];

      if ((phi_ref > -9999 ) && (dzOverHelPitch < _dzOverHelPitchCut) && (dz > tollMin)) {
	dphi  = phiCl - phi_ref;
	dphi  = TVector2::Phi_0_2pi(dphi);
//-----------------------------------------------------------------------------
// add 2 pi for taking into account the fact we are in the second loop
// FIX ME: what to do if we are in the third loop?
//-----------------------------------------------------------------------------
	if (dz > tollMax) dphi += 2*M_PI*int(dz/tollMax);

	double dphidz = dphi/dz;
	while (dphidz < 0.) {
	  dphi   += 2.*M_PI;
	  dphidz = dphi/dz;
	}

	double tmpphi0 = phi_ref - dphidz*z_ref;
	tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

	if (_debug >5){
	  printf("[CalHelixFinderAlg::findDfDz] z_ref = %9.3f z = %9.3f phi_ref = %9.5f phi = %9.5f dz = %10.3f dz/HelPitch = %10.3f df/dz = %9.5f phi0 = %9.6f\n",
		 z_ref, zCl, phi_ref, dphi-phi_ref, dz, dzOverHelPitch, dphidz, tmpphi0);
	}

	if (dzOverHelPitch < _dzOverHelPitchCut ) {
	  _hDfDzRes->Fill(dphidz, weight_cl);
	  if ((dphidz < _minDfDz) || (dphidz >  _maxDfDz)) tmpphi0 = -1;
	  _hPhi0Res->Fill(tmpphi0, weight_cl);
	}
      }
    }
//-----------------------------------------------------------------------------
// 2015 - 04- 02 G. Pezzu changed the way the maximum is searched
// since sometimes a 2pi ambiguity creates two peaks in the histogram
// we want to use the second, because it is the correct one
//-----------------------------------------------------------------------------
    double  maxContent = _hDfDzRes->GetMaximum() - 0.001;
    int      maxBin    = _hDfDzRes->FindLastBinAbove(maxContent);//GetMaximumBin();
    _hdfdz             = _hDfDzRes->GetBinCenter(maxBin);//_hDfDzRes->GetMean();
    double dfdzmean    = _hDfDzRes->GetMean();
    int    nentries    = _hDfDzRes->GetEntries();
    int    overflows   = _hDfDzRes->GetBinContent(0)  + _hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1);

    maxContent         = _hPhi0Res->GetMaximum() - 0.001;
    maxBin             = _hPhi0Res->FindLastBinAbove(maxContent);//GetMaximumBin();

    double mpvphi0     = _hPhi0Res->GetBinCenter(maxBin); //_hPhi0Res->GetMean();
    double menaphi0    = _hPhi0Res->GetMean();
    int    nentriesphi = _hPhi0Res->GetEntries();

    _hphi0 = mpvphi0;  // 2018-01-05: *DOUBLE_CHECK*

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::findDfDz:DFDZ] nent: %3i mpvDfDz: %9.6f meanDphiDz: %9.6f under: %3.0f over: %3.0f ENTRIES:",
	     nentries, _hdfdz, dfdzmean,
	     _hDfDzRes->GetBinContent(0),_hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1)
	     );
      for (int i=0; i<_hDfDzRes->GetNbinsX(); i++) {
	printf(" %3.0f",_hDfDzRes->GetBinContent(i+1));
      }
      printf("\n");

      printf("[CalHelixFinderAlg::findDfDz:PHI0] nent: %3i mpvPhi0: %9.6f meanPhi0  : %9.6f under: %3.0f over: %3.0f ENTRIES:",
	     nentriesphi, mpvphi0,  menaphi0,
	     _hPhi0Res->GetBinContent(0),_hPhi0Res->GetBinContent(_hPhi0Res->GetNbinsX()+1)
	     );
      for (int i=0; i<_hPhi0Res->GetNbinsX(); i++) {
	printf(" %3.0f",_hPhi0Res->GetBinContent(i+1));
      }
      printf("\n");
    }
//-----------------------------------------------------------------------------
// Part 2: try to perform a more accurate estimate - straight line fit
//-----------------------------------------------------------------------------
    double z0, phi0, dphidz, pred;

    z0     = 0.    ;
    phi0   = _hphi0;
    dphidz = _hdfdz;
    //    _sdfdz = -1;

    if (_debug > 5) {
      double tmpphi0=phi0+dphidz*z0;
      printf("[CalHelixFinderAlg::findDfDz:PART2] phi0 = %9.6f dfdz = %9.6f\n", tmpphi0, dphidz);
    }
//--------------------------------------------------------------------------------
// 2015-03-25 G. Pezzu changed theway the 2Phi ambiguity is resolved
//--------------------------------------------------------------------------------
    LsqSums4 srphi;

    weight = 1./(_sigmaPhi*_sigmaPhi);

    double xdphi, zLast(z0), zdist;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::findDfDz:LOOP]  i       z        dz      phiVec     phi    srphi.dfdz    dphi    xdphi     dfdz    \n");
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

    double errCl    = 2./30.;
    double weightCl = 1./(errCl*errCl);

    xdphi = fabs(dphi)/errCl;
//-----------------------------------------------------------------------------
// 2015-04-21 Gianipez added the condition (xdphi < 2.*_maxXDPhi) for adding or not
// the calorimeter point to the fitter. In case the particle scattered in the end of the tracker
// the calorimeter point is dangerous.
//-----------------------------------------------------------------------------
    if (xdphi < 2.*_maxXDPhi){
      srphi.addPoint(zCl, phiCl, weightCl);
    }

    if (_debug > 10) {
      printf("[CalHelixFinderAlg::findDfDz:LOOP] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
	     0, zCl, dz, phiCl, phiCl, srphi.dfdz(), dphi, xdphi, dphidz);
    }

    // 2015-07-06 Gianipez added the following line for avoiding infinite loops
    if ( i0 < 0 ) goto NEXT_STEP;

    for (int i=i0; i<nstations; i++) {
      if (nhits[i] > 0) {
	z    = zVec[i];
	dz   = z-z0;
	pred = phi0 + dz*dphidz;
	phi  = phiVec[i];
	dphi = phi - pred;//pred - phi;

	while (dphi > M_PI){
	  phi -= 2*M_PI;//+= 2*M_PI;
	  dphi = phi - pred;//pred - phi;
	}
	while (dphi < -M_PI){
	  phi += 2*M_PI;//-= 2*M_PI;
	  dphi = phi - pred;//pred - phi;
	}

	xdphi = fabs(dphi)/_sigmaPhi;

	if (xdphi < 2.*_maxXDPhi){
	  srphi.addPoint(z, phi, weight);

	  zdist = z - zLast;

	  if ( (srphi.qn() >= 3.) && (zdist > 500.)){
	    z0     = 0.;
	    phi0   = srphi.phi0();
	    dphidz = srphi.dfdz();
	  }
	}

	if (_debug > 10) {
	  double tmpDfDz = srphi.dfdz();//, Helix._srphi.chi2DofLine());
	  printf("[CalHelixFinderAlg::findDfDz:LOOP] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
		 i, z, dz, phiVec[i], phi, tmpDfDz, dphi, xdphi, dphidz);
	}
      }
    }

  NEXT_STEP:;

    if (srphi.qn() >= 3.) {
      _hdfdz = srphi.dfdz();		// sigxy/sigxx;
      _hphi0 = srphi.phi0();		// ymean - xmean*sigxy/sigxx;
      //      _sdfdz = srphi.chi2DofLine();
    }
    else {
      _hphi0 = phi0 + _hdfdz*z0;
      //      _sdfdz = -1;
    }

    int     nActive_stations = nentries - overflows;
    
    if (Diag_flag > 0){
      Helix._diag.nStationPairs = nActive_stations;
    }

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::findDfDz] END: _hdfdz = %9.5f _hphi0 = %9.6f chi2 = %9.3f ", _hdfdz,
	     _hphi0, -1.);
      printf(" FIT: srphi.dfdz() = %9.5f srphi.phi0() = %9.6f chi2 = %9.3f qn = %6.0f\n", srphi.dfdz(),
	     srphi.phi0(), srphi.chi2DofLine(), srphi.qn());
    }
//----------------------------------------------------------------------------- 
// 2017-11-14 gianipez: in case of less than _minNSt active stations we should 
//                      probably use the _mpDfDz
//-----------------------------------------------------------------------------
    if (nActive_stations < _minNSt) {
      _hdfdz = _mpDfDz;
      return 0;
    }

    return 1;
  }




//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::doLinearFitPhiZ(CalHelixFinderData& Helix    ,
					  int                 SeedIndex,
					  int                *IndexVec ,
					  int                 UseInteligentWeight,
					  int                 DoCleanUp           ) {
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
    double z_previous(-99999), weight_previous(0), xdphi_previous(-9999);
    int    i_previous(-1);
    double dphi, err, xdphi;

    CLHEP::Hep3Vector helCenter = Helix._center;
    double            radius    = Helix._radius;

    CLHEP::Hep3Vector pos;

    // calculate the mean pitch assuming conversion electron
    // hypothesis (so assuming a pitch angle of ~ 0.65)

    ///    double pitch = twopi*std::pow(Helix._radius,2.)*Helix._dfdz; //*tan(0.65);


//-----------------------------------------------------------------------------
// gianipez: procedure for aligning the phi vector
//-----------------------------------------------------------------------------
    ::LsqSums4 srphi;
    int        iworst, count(0), indexWorst;
    double     chi2, chi2min, deltaPhi, dphi_max(0), phi_ref, weightWorst(0);


    if (Helix._sxy.qn() > 0) {
      helCenter = CLHEP::Hep3Vector( Helix._sxy.x0(), Helix._sxy.y0(), 0);
      radius    = Helix._sxy.radius();
    }
    CLHEP::Hep3Vector strawDir;
    const char        banner[200] = "doLinearFitPhiZ";
//--------------------------------------------------------------------------------
// set EMC cluster info and initilize the dfdz for the search
//-----------------------------------------------------------------------------
    double dfdz  = Helix._dfdz;
    double phi0  = Helix._fz0;

    double zCl   = fCaloZ;
    pos          =  Hep3Vector(fCaloX, fCaloY, fCaloZ);
    double phiCl = CLHEP::Hep3Vector(pos - helCenter).phi();//center).phi();
    phiCl        = TVector2::Phi_0_2pi(phiCl);

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

    Helix._srphi.clear();
    if ( xdphi < 2.*_maxXDPhi ) {
      double weight_cl = 10.0;
      Helix._srphi.addPoint(zCl,phiCl,weight_cl);
    }

    count = 0;
    double zlast, dz;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ:BEGIN] phi0 = %10.6f dfdz = %10.6f chi2N = %10.3f DoCleanup = %i\n",
	     Helix._fz0,  Helix._dfdz, 0.,DoCleanUp);
      printf("[CalHelixFinderAlg::doLinearFitPhiZ]    flag   A   shID   i       z         ");
      printf("    phi         dphi      xdphi      zlast        dz      dphidz  szphidfdz  chi2\n");
      printf("[CalHelixFinderAlg::doLinearFitPhiZ] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
	     0, 1, 0, 0,  zCl, phiCl, deltaPhi, xdphi, 0., 0., dfdz, 0., 0.);
    }

    zlast = 0;
    for (int i=SeedIndex; i<N; ++i) {
      pos                  = _xyzp[i]._pos;
      z                    = pos.z();
      strawDir             = _xyzp[i]._sdir;

      phi      = CLHEP::Hep3Vector(pos - helCenter).phi();//center).phi();
      phi      = TVector2::Phi_0_2pi(phi);
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
      phi_corrected[i] = phi;
      _phiCorrected[i] = phi;

      dphi             = fabs(dphi);
      err              = _sigmaPhi;

      if (UseInteligentWeight == 1){
	weight = calculatePhiWeight(pos, strawDir, helCenter, radius, 0, banner);
	err    = 1./sqrt(weight);
      }

      xdphi = dphi/err;

      if (_debug > 10) {
	printf("[CalHelixFinderAlg::doLinearFitPhiZ:LOOP] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
	       *((int*) &_xyzp[i]._flag), idVec[i] < 0 ? 0 : idVec[i],
	       _xyzp[i]._strawhit->strawIndex().asInt()/*int(_xyzp[i]._ind)*/, i,
	       z, phi_corrected[i], dphi,xdphi,zlast,dz,
	       dfdz, Helix._srphi.dfdz(), Helix._srphi.chi2DofLine());
      }

      if (_xyzp[i].isOutlier())                        continue;
      if (  idVec[i] < 1     )                         continue;

      if ( (DoCleanUp == 1) && (xdphi > _maxXDPhi) ) {
	idVec[i] = 0;
	++nPointsRemoved;
	                                               continue;
      }

      if ( fabs(z-z_previous) > 1e-3){
	Helix._srphi.addPoint(z, phi_corrected[i], weight);
	++count;

	//set the values
	i_previous      = i;
	z_previous      = z;
	weight_previous = weight;
	xdphi_previous  = xdphi;
      }else {//the point is in the same layer of the previous one: make a choice
	if (xdphi < xdphi_previous){
	  //remove the point previously added
	  Helix._srphi.removePoint(z_previous, phi_corrected[i_previous], weight_previous);
	  idVec[i_previous] = 0;

	  //add the new point 
	  Helix._srphi.addPoint   (z         , phi_corrected[i]         , weight         );

	  //set the values
	  i_previous      = i;
	  z_previous      = z;
	  weight_previous = weight;
	  xdphi_previous  = xdphi;
	}else {//a point in the same layer closer to predictio was previously added
	  idVec[i] = 0;
	}
      }

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

    _phiCorrectedDefined = 1;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ:BEFORE_CLEANUP] Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f points removed = %4i\n",
	     Helix._srphi.phi0(),Helix._srphi.dfdz(), Helix._srphi.chi2DofLine(), nPointsRemoved);
    }
//-----------------------------------------------------------------------------
// perform a cleanup in RZ
//-----------------------------------------------------------------------------
    if ( DoCleanUp == 1){
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
	  //printf("[CalHelixFinderAlg::doLinearFitPhiZ] chi2 = %5.3e chi2min = %5.3e\n", chi2, chi2min);
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
	    printf("[CalHelixFinderAlg::doLinearFitPhiZ_removed:LOOP2] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);
	  }
	}

      CHECK_RESIDUALS:;
	dphi_max    = _maxXDPhi;
	iworst      = -1;
	weightWorst = -1;

	for(int i=SeedIndex; i < N; ++i) {
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
	    printf("[CalHelixFinderAlg::doLinearFitPhiZ:REMOVED] %6i %5.3f     %5.3f chi2 = %5.3f  \n", indexWorst, z, phi, chi2min);
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
    }
//-----------------------------------------------------------------------------
// 2015-04-21 Gianipez changed the threshold from 3 to _minNHits. there is no reason
// which should allow to keep a result which selects a number of points lower than the threshold!
//-----------------------------------------------------------------------------
    if ( (Helix._srphi.qn() >= _minNHits) && (Helix._srphi.chi2DofLine() < _chi2zphiMax) ){
      success = true;
    }
    //----------------------------------------------------------------------//
    if (Helix._srphi.dfdz() < 0.) { // *FIXME* : negative helicity handling
      success = false;
    }
    else if (success) {                               // update helix results
      Helix._fz0  = Helix._srphi.phi0();
      Helix._dfdz = Helix._srphi.dfdz();
    }

    if (SeedIndex == 0) {
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
      Helix._diag.phi0_6           = Helix._srphi.phi0();
      Helix._diag.rdfdz_7          = Helix._srphi.dfdz()* Helix._diag.n_rescued_points_9;
      Helix._diag.dfdz_8           = Helix._srphi.dfdz();
      Helix._diag.chi2_dof_line_13 = Helix._srphi.chi2DofLine();

      if (success) {
	int h=0;
	for (int i=SeedIndex; i<N; ++i){
	  if (_xyzp[i].isOutlier())     continue;
	  if (  idVec[i] < 1  )         continue;
	  pos      = _xyzp[i]._pos;
	  z        = pos.z();
	  phi      = z* Helix._dfdz + Helix._fz0;
	  deltaPhi = phi_corrected[i] - phi;

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

    _chi2nFindZ = Helix._srphi.chi2DofLine();
    if(_chi2nFindZ < 0.0)  _eventToLook = 0; // FIXME!!

    if (_debug > 5) {

      printf("[CalHelixFinderAlg::doLinearFitPhiZ:END] retval = %d Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f srphi: phi_0 = %5.3f dfdz = %5.5f\n",
	     success ? 1:0, Helix._srphi.phi0(),Helix._srphi.dfdz(), Helix._srphi.chi2DofLine(), srphi.phi0(), srphi.dfdz() );

      if (_debug > 10) {
	printf("[CalHelixFinderAlg::doLinearFitPhiZ:END2]    flag   A   shID       z             phi      phi-dfdz*z-phi0\n");

	for (int i=N-1; i>=SeedIndex; --i) {
	  pos      = _xyzp[i]._pos;
	  z        = pos.z();
	  phi      = z* Helix._dfdz + Helix._fz0;
	  deltaPhi = phi_corrected[i] - phi;

	  printf("[CalHelixFinderAlg::doLinearFitPhiZ:END2] %08x %2i %6i %12.5f %12.5f %12.5f\n",
		 *((int*) &_xyzp[i]._flag), idVec[i] < 0 ? 0 : idVec[i],
		 _xyzp[i]._strawhit->strawIndex().asInt(),  z, phi_corrected[i], deltaPhi);
	}
      }
    }

    if (success) {
      for (int i=0; i<N; ++i) {
	IndexVec[i] = idVec[i];
      }
    }
    return success;
  }


  void  CalHelixFinderAlg::fillHitLayer(CalHelixFinderData& Helix) {
    
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
// 	  const Straw*            straw = &_tracker->getStraw(sh->strawIndex());
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
  }

//-----------------------------------------------------------------------------
// 12-09-2013 gianipez modified this procedure to avoid the doubling of the
// same stereohitposition
// points in filled array are ordered in Z coordinate
//-------------------------------------------------------------------------
  void CalHelixFinderAlg::fillXYZP(CalHelixFinderData& Helix) {

    _xyzp.clear();

    const vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();

    int size = Helix._timeCluster->nhits();
//--------------------------------------------------------------------------------
    if (Helix.shpos() != 0) {
      int loc;
      StrawHitFlag flag;
      for (int i=0; i<size; ++i) {
	loc                = shIndices[i];	 // index in shcol of i-th timecluster hit
	flag               = Helix.shfcol()->at(loc);
//-----------------------------------------------------------------------------
// select hits: don't reuse straw hits
//-----------------------------------------------------------------------------
	int good_hit = flag.hasAllProperties(_hsel  );
	int bkg_hit  = flag.hasAnyProperty  (_bkgsel);
	int used_hit = flag.hasAnyProperty  (StrawHitFlag::calosel);

	if (good_hit && (! bkg_hit) && (! used_hit)) {

	  const StrawHit& sh          = Helix.shcol()->at(loc);
	  const Straw& straw          = _tracker->getStraw(sh.strawIndex());
	  const StrawHitPosition& shp = Helix.shpos()->at(loc);

	  if (sh.energyDep() > _maxElectronHitEnergy)         continue;

	  CalHelixPoint pos(loc,sh,shp,straw,flag);
	  _xyzp.push_back(pos);
	}
      }
    }
//----------------------------------------------------------------------
// 2014-11-06 gianipez added the following line for ordering the xyzp
// strawhits along their z coordinate
//----------------------------------------------------------------------
    std::sort(_xyzp.begin(), _xyzp.end(), [ ]( const CalHelixPoint& lhs,   // 2018-01-04 PM *FIXME* sort hits in advance
					       const CalHelixPoint& rhs )
	      {
		return lhs._pos.z() < rhs._pos.z();
	      } );

    if (_debug > 0) printXYZP(Helix);
  }




//----------------------------------------------------------------------------
//2015-01-17 G. Pezzullo: the following procedure looks the hit with
// z-coordinate smaller then the seeding one and calculates distance from
// prediction in order to check if they are good or outliers
//----------------------------------------------------------------------------
  void CalHelixFinderAlg::rescueHitsBeforeSeed(CalHelixFinderData& Helix){

    double      weight, radius, phi0, dfdz, x0, y0;
    dfdz        = Helix._dfdz;
    //    phi0        = CLHEP::Hep3Vector(_xyzp[Helix._seedIndex]._pos - Helix._center).phi();
    phi0        = Helix._fz0 + dfdz*(_xyzp[Helix._seedIndex]._pos.z());
    x0          = Helix._center.x();
    y0          = Helix._center.y();
    radius      = Helix._radius;

    double      dx,dy,phi,max_dist;
    Hep3Vector  shPos, hePos, strawDir, helCenter(x0, y0, 0);

    double deltaZ(0.); // , deltaX(0.), deltaY(0.);
    double distXY(0.);
    double dist  (0.), dist2(0.); // help parameter for storing strawhit position residual
    int    i_last(Helix._seedIndex), rescuedPoints(0);

    char banner[]="CalHelixFinderAlg::rescueHitsBeforeSeed";

    if (_debug > 0) {
      printf("[%s:BEGIN] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.5f dfdz = %5.6f chi2 = %5.3f \n", banner,
	     x0, y0, radius, phi0, dfdz , Helix._sxy.chi2DofCircle());
      printf("[%s] SeedIndex = %i N-points = %5.3f\n",  banner, Helix._seedIndex, Helix._sxy.qn()-1);
      if (Helix._seedIndex > 0) {
	printf("[%s] index      Z        xi      yi       xp       yp       X0        Y0         R        dfdZ  dXY(pred) dXY(seed) dZ(seed)     wt\n",banner);
	printf("[%s]-------------------------------------------------------------------------------------------------------------------------------\n",banner);
      }
    }
//-----------------------------------------------------------------------------
// given a helix candidate, move upstream and pick up points with Z < _xyzp[fSeedIndex].z
//-----------------------------------------------------------------------------
    for (int i=Helix._seedIndex-1; i>=0; --i){
      if (_xyzp[i].isOutlier())                              continue;
      shPos     = _xyzp[i]._pos;
      strawDir  = _xyzp[i]._sdir;

      // 2017-09-26 gianipez changed the weight from 1 to "inteligent".
      weight    = calculateWeight(shPos, strawDir, helCenter, radius);

      deltaZ    = shPos.z() - _xyzp[i_last]._pos.z();
      phi       = phi0 + (deltaZ)*dfdz;                     // tanLambda/radius;

      hePos     = Hep3Vector(x0 + radius*std::cos(phi),
			     y0 + radius*std::sin(phi),
			     shPos.z());

      dx    = hePos.x() - shPos.x();
      dy    = hePos.y() - shPos.y();
      dist2 = dx*dx + dy*dy;
      dist  = std::sqrt(dist2);

      if (_debug > 10) {
	printf("[%s:LOOP] %5i %9.3f %8.3f %8.3f %8.3f %8.3f %9.3f %9.3f %9.3f %8.5f %8.3f %8.3f %8.3f %8.3f",
	       banner,i,shPos.z(),shPos.x(),shPos.y(),hePos.x(),hePos.y(),
	       x0,y0,radius,dfdz,dist,distXY,deltaZ,weight);
      }

      max_dist = _distPatRec + _dfdzErr*fabs(deltaZ);
      if (dist <= max_dist) {
	i_last = i;			// index of last good point
	//	phi0   = phi;                   // CLHEP::Hep3Vector(shPos - CLHEP::Hep3Vector(x0, y0, 0.)).phi();
					// add point to the helixfithack result objet
	Helix._sxy.addPoint(shPos.x(),shPos.y(),weight);

					// mark good point
	_indicesTrkCandidate[i] = 1;
					// store distance along z-axis from the last point found relying on thr helix
	_dzTrkCandidate     [i] = deltaZ;

					// store distance from predition
	_distTrkCandidate   [i] = dist;
					// update helix parameters
	x0      = Helix._sxy.x0();
	y0      = Helix._sxy.y0();
	radius  = Helix._sxy.radius();

	helCenter.setX(x0);
	helCenter.setY(y0);
	phi0    = CLHEP::Hep3Vector(shPos - helCenter).phi();
	++rescuedPoints;

	if (_debug > 0) {
	  printf("rescued %08x %2i %12.5f %12.5f %12.5f\n",
		 *((int*) &_xyzp[i]._flag), _indicesTrkCandidate[i],
		 _xyzp[i]._pos.x(), _xyzp[i]._pos.y(), _xyzp[i]._pos.z());
	}
      }
      else {
	if (_debug > 10) {
	  printf(" missed\n");
	}
      }
    }
//-----------------------------------------------------------------------------
// update Helix info
//-----------------------------------------------------------------------------
    Helix._center.set(x0, y0, 0.0);
    Helix._radius = radius;

    _goodPointsTrkCandidate = Helix._sxy.qn() - 1; // removing the EMC cluster!

    if (_debug > 5) {
      printf("[%s:END] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.3f dfdz = %5.6f chi2 = %5.3f",
	     banner,
	     x0, y0, radius,  Helix._fz0, dfdz , Helix._sxy.chi2DofCircle());
      printf(" SeedIndex: %i N(rescued points): %i\n",Helix._seedIndex,rescuedPoints);
    }

    Helix._diag.n_rescued_points_16 = rescuedPoints;
  }

//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::printInfo(CalHelixFinderData&  Helix){
    const char banner [] = "CalHelixFinderAlg::printInfo";
    double dr(0), dx, dy;

    if (_debug > 0) {
      printf("[%s] N(points): %3.0f x0: %12.5f y0: %12.5f r: %12.5f chi2c: %12.5f phi0: %5.5f dfdz: %5.6f chi2l: %5.3f\n",
	     banner,
	     Helix._sxy.qn() - 1,
	     Helix._sxy.x0(),Helix._sxy.y0(),Helix._sxy.radius(),Helix._sxy.chi2DofCircle(),
	     Helix._fz0, Helix._dfdz , Helix._srphi.chi2DofLine());

      int np = _xyzp.size();
      for (int i=0; i<np; i++) {
	dx = _xyzp[i]._pos.x() - Helix._sxy.x0();
	dy = _xyzp[i]._pos.y() - Helix._sxy.y0();
	dr = sqrt(dx*dx+dy*dy) - Helix._sxy.radius();
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
  void CalHelixFinderAlg::filterUsingPatternRecognition(CalHelixFinderData& Helix) {

    if (Helix._seedIndex < 0) return;

    int    nActive(0), nActive_hel(0);
    int    np = _xyzp.size();

    double      straw_mean_radius(0), chi2_global_helix(0), total_weight(0);
    double      x_center(Helix._center.x()), y_center(Helix._center.y()), radius(Helix._radius), fz0(Helix._fz0), lambda(1./Helix._dfdz);
    Hep3Vector  hel_pred(0);

    for (int i=0; i<np; ++i) {
      if (_debug > 0) {
	if (i == 0) {
	  Hep3Vector pSeed = _xyzp[Helix._seedIndex]._pos;
					// if no second strawhit has been found, use the target center in the transverse plane
	  // Hep3Vector pCand(0.0, 0.0, 0.0);
	  // if (fCandIndex >= 0) pCand = _xyzp[fCandIndex]._pos;

	  printf("[CalHelixFinderAlg::filterUsingPatternRecognition]  filterUsingPatternRecognition() will set asOutlier the following hits using helix parameters\n");
	  printf("[CalHelixFinderAlg::filterUsingPatternRecognition] X0 = %5.3f Y0 = %5.3f r = %5.3f chi2N = %5.5f phi0 = %5.5f dfdz = %5.5f chi2N = %5.5f straw-hits = %i\n",
		 Helix._sxy.x0(), Helix._sxy.y0(), Helix._radius, Helix._sxy.chi2DofCircle(),
		 Helix._srphi.phi0(), Helix._srphi.dfdz(), Helix._srphi.chi2DofLine(),
		 _goodPointsTrkCandidate );// +1 for counting also the seeding hit
	  printf("[CalHelixFinderAlg::filterUsingPatternRecognition] index  shID type           X        Y         Z        dist      Dz\n");
	  printf("[CalHelixFinderAlg::filterUsingPatternRecognition] -------------------------------------------------------------------\n");
	}
      }

      double dist = _distTrkCandidate[i];
      double dz   = _dzTrkCandidate  [i];

      if (_indicesTrkCandidate[i] <= 0) _xyzp[i].setOutlier(); // *FIXME* ? ? do we want to call this an outlier ?
      else  {
	++nActive;
	double  x         = _xyzp[i]._pos.x();
	double  y         = _xyzp[i]._pos.y();
	double  z         = _xyzp[i]._pos.z();
	double  phi_pred  = fz0 + z/lambda;
	double  x_pred    = x_center + radius*cos(phi_pred);
	double  y_pred    = y_center + radius*sin(phi_pred);
	hel_pred.setX(x_pred);
	hel_pred.setY(y_pred);
	double  weight    = calculateWeight(_xyzp[i]._pos, _xyzp[i]._sdir, hel_pred, radius)*_weight3D;
	double  x_resid2  = (x - x_pred)*(x - x_pred);
	double  y_resid2  = (y - y_pred)*(y - y_pred);
	double  hitResi2  = (x_resid2 + y_resid2)*weight;
	if (hitResi2 > _chi2hel3DMax) {
	  _xyzp[i].setOutlier();  // *FIXME* ? ? do we want to call this an outlier ?
	} else {
	  ++nActive_hel;
	  chi2_global_helix = chi2_global_helix + hitResi2;
	  straw_mean_radius = straw_mean_radius + sqrt(x*x + y*y)*weight;
	  total_weight      = total_weight + weight;	
	}
      }

      if (_debug > 10) {
	Hep3Vector* shPos = &_xyzp[i]._pos;
	int is_outlier    = _xyzp[i].isOutlier();
	string type;
	if      (i == Helix._seedIndex) type = "seed";
	else if (i == Helix._candIndex) type = "cand";

	printf("[CalHelixFinderAlg::filterUsingPatternRecognition] %5i %5i %4i %4s  %8.3f %8.3f %9.3f %8.3f %8.3f\n",
	       i,_xyzp[i]._strawhit->strawIndex().asInt(),is_outlier,type.data(),shPos->x(),shPos->y(),shPos->z(),dist,dz);
      }
    }

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::filterUsingPatternRecognition:END] N( active points) = %i over N-hits = %i\n", nActive, np);
    }

    _goodPointsTrkCandidate = nActive_hel; // nActive

    Helix._diag.n_active_11       = nActive_hel;
    Helix._diag.straw_mean_radius = (nActive > 0) ? straw_mean_radius/total_weight : -1;
    Helix._diag.chi2d_helix       = (nActive > 0) ? chi2_global_helix/double(nActive_hel - 5.): -1;
  }

//-----------------------------------------------------------------------------
// leave only hits in the same hemisphere with the cluster
// 
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::filterDist() {

    static const double pi(M_PI);
    static const double twopi(2*pi);

    double clPhi(-9999.);

    int np = _xyzp.size();
    if (fCaloTime > 0) clPhi = atan2(fCaloY,fCaloX);

    if (_debug > 0) {
      printf("[CalHelixFinderAlg::filterDist]----------------------------------------------------------------------------------------------\n");
      printf("[CalHelixFinderAlg::filterDist]     i     flag        X            Y         Z            phi    nhits: %5i clusterPhi: %8.4f\n",np,clPhi);
      printf("[CalHelixFinderAlg::filterDist]----------------------------------------------------------------------------------------------\n");
    }

    int nfiltered(0);
    for (int ip=0; ip < np; ip++) {
      CalHelixPoint* hit = &_xyzp[ip];
      if (hit->use()) {

	double dphi = hit->_phi - clPhi;

	if (dphi >  pi) dphi -= twopi;
	if (dphi < -pi) dphi += twopi;

	if (fabs(dphi) > pi/2){
	  hit->setOutlier();
	  nfiltered++;
	}
      }
      if (_debug > 0) {
	printf("[CalHelixFinderAlg::filterDist] %5i %08x %12.5f %12.5f %12.5f %12.5f \n",
	       ip, *((int*) &hit->_flag), hit->x(), hit->y(),hit->z(),hit->_phi);
      }
    }

    if (_debug > 0) {
      printf("[CalHelixFinderAlg::filterDist:END]: nfiltered: %5i\n",nfiltered);
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

//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::doPatternRecognition(CalHelixFinderData& Helix) {

    int useMPVdfdz(0), useDefaultDfDz(false);
    int np = _xyzp.size();

    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition:BEGIN] fUseDefaultDfDz = %i\n",fUseDefaultDfDz);

    if (_debug2 == 0){
      _debug2 = _debug;
      _debug  = 0;
    }

    for (int i=0; i<np; i++) {
      if (_xyzp[i].isOutlier())                             continue;
//----------------------------------------------------------------------
// 2014-12-26 gianipez: don't start the search from an already used hit
// 2018-01-03 PM : for a long time, hits are never marked as used
//-----------------------------------------------------------------------------
//      if ( isHitUsed(i) == 1 )                              continue;

      if ((np-i) > _goodPointsTrkCandidate) {
	if (_debug > 5) {
	  printf("[CalHelixFinderAlg::doPatternRecognition]: calling findTrack(i=%i,Helix,useDefaltDfDz=FALSE,useMPVdfdz=%i)",i,useMPVdfdz);
	  printf(" : np=%3i _goodPointsTrkCandidate=%3i useDefaultDfDz = %i\n",np,_goodPointsTrkCandidate,useDefaultDfDz);
	}
	findTrack(i,Helix,useDefaultDfDz,useMPVdfdz);
      }
//------------------------------------------------------------------------------
// 2015-01-22 P.Murat: what happens when the very first candidate is good enough ?
//                     where is the comparison of the found candidate with the best previous one ?
// 2018-01-04 PM *FIXME*: at this point, if there is a good enough candidate, the algorithm should 
//               try to complete the phi-z search before reiterating on XY 
//-----------------------------------------------------------------------------
    }

    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition:STEP1] useDefaultDfDz = %i\n",useDefaultDfDz);
//-----------------------------------------------------------------------------
// 2014-11-09 gianipez: if no track was found requiring the recalculation of dfdz
// look for a track candidate using the default value of dfdz and the target center
//-----------------------------------------------------------------------------
    if (fUseDefaultDfDz == 0) {
      int useMPVdfdz = 0;
      for (int i=0; i<np; i++) {
	if (_xyzp[i].isOutlier())                           continue;
	if ((np -i) > _goodPointsTrkCandidate) {
	  if (_debug > 5) {
	    printf("[CalHelixFinderAlg::doPatternRecognition]: calling findTrack(i=%i,Helix,useDefaltDfDz=TRUE,useMPVdfdz=%i)",i,useMPVdfdz);
	    printf(" : np=%3i _goodPointsTrkCandidate=%3i\n",np,_goodPointsTrkCandidate);
	  }
	  findTrack(i,Helix,true,useMPVdfdz);
	}
      }
    }
//-----------------------------------------------------------------------------
// 2015-01-14 G. Pezzullo added the findDfDz procedure
//-----------------------------------------------------------------------------
    if (_debug > 5) printf("[CalHelixFinderAlg::doPatternRecognition]: ------------ calling findDfDz\n");

    int    diag_flag(1);
    if (Helix._seedIndex >= 0) {
      findDfDz(Helix,Helix._seedIndex,_indicesTrkCandidate, diag_flag);
      if (_debug > 5) printf("[CalHelixFinderAlg::doPatternRecognition]: findDfDz ----> phi0 = %5.5f dfdz = %5.5f \n",
			      _hphi0, _hdfdz);
    }
    else {
//-----------------------------------------------------------------------------
// Helix._seedIndex < 0 means that no candidate has been found
// usually it happens when the cluster is not on the trajectory or the time peak 
// has very low number of hits
// maybe we should set a threshold on the time peak size to avoid such?
//-----------------------------------------------------------------------------
      int vIndices[np];
      for (int i=0; i<np; ++i) vIndices[i] = 1;

      findDfDz(Helix, 0, vIndices, diag_flag);
      if (_debug > 5) {
	printf("[CalHelixFinderAlg::doPatternRecognition]: findDfDz called using SeedIndex = 0 and using all hits (expect outliers!) \n");
	printf("[CalHelixFinderAlg::doPatternRecognition]: findDfDz ----> phi0 = %5.5f dfdz = %5.5f \n",
	       _hphi0, _hdfdz);
      }
    }
//-----------------------------------------------------------------------------
// 3rd loop - what is its role? 
//-----------------------------------------------------------------------------
    useMPVdfdz = 1;
    for (int i=0; i<np; i++) {
      if (_xyzp[i].isOutlier())                             continue;
      if ((np -i) > _goodPointsTrkCandidate) {
	if (_debug > 5) { 
	  printf("[CalHelixFinderAlg::doPatternRecognition]: calling findTrack(i=%i,Helix,useDefaltDfDz=FALSE,useMPVdfdz=%i)",i,useMPVdfdz);
	  printf(" : np=%3i _goodPointsTrkCandidate=%3i\n",np,_goodPointsTrkCandidate);
	}
	findTrack(i,Helix,false,useMPVdfdz);
      }
    }

    if (_debug == 0){
      _debug  = _debug2;
      _debug2 = 0;
    }

    char banner[200];
    bool rc;
    int  rc1, refineHelParamRes, rs, usePhiResid, useInteligentWeight(1);

    if ((Helix._seedIndex < 0) || (Helix._sxy.qn() < 5) )   goto  PATTERN_RECOGNITION_END;

    // 2015-01-17 G. Pezzullo: rescue points with z-coordinate less than the seed hit
    if (_debug != 0) {
      printf("[CalHelixFinderAlg::doPatternRecognition]: calling rescueHitsBeforeSeed\n");
      printInfo(Helix);
    }

    rc = doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate, useInteligentWeight);

    //2017-10-05 Gianipez added the following line to make some tests
    if ((_smartTag == 1) && (Helix._srphi.qn() == 0.))      goto  PATTERN_RECOGNITION_END;

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

    refineHelParamRes = refineHelixParameters(Helix, 0, _indicesTrkCandidate, banner, _debug);

    if (refineHelParamRes >= 0) {
      Helix._center.setX(Helix._sxyw.x0());
      Helix._center.setY(Helix._sxyw.y0());
      Helix._radius    = Helix._sxyw.radius();
      Helix._sxy.init(Helix._sxyw);
      if (_debug != 0)  printInfo(Helix);
    }
//---------------------------------------------------------------------------------------
// use the results of the helix search to see if points along the track can be rescued
//---------------------------------------------------------------------------------------
    if ((Helix._srphi.qn() < 10) || (!rc)) usePhiResid = 0;
    else                                   usePhiResid = 1;

    rescueHits(Helix, 0, _indicesTrkCandidate, usePhiResid);
    
    if (Helix._sxy.qn() != Helix._srphi.qn()) rc = doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate, useInteligentWeight);

    if (_debug != 0)  printInfo(Helix);
//--------------------------------------------------------------------------------------------------------------
// 2015-03-25 G. Pezzu added the following call to findDfDz(...) in order to help the fitter on finding
// the more reliable value of dfdz which is needed for resolving the 2pi ambiguity.
// Since in the previous step we could have rescued few points, that would give us an help!
// re-evaluate the df/dz and phi0 including rescued hits and new XY parameters
//--------------------------------------------------------------------------------------------------------------
    if (Helix._srphi.qn() < 10 || (!rc)){
      rs = findDfDz(Helix, 0, _indicesTrkCandidate);
      
      if (rs == 1) {			// update Helix Z-phi part
	Helix._dfdz = _hdfdz;
	Helix._fz0  = _hphi0;
      }
    }

    rc = doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate, useInteligentWeight);

    if (rc) {
      usePhiResid = 1;
      rescueHits(Helix, 0, _indicesTrkCandidate, usePhiResid);
      if (Helix._sxy.qn() != Helix._srphi.qn()) rc = doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate, useInteligentWeight);

      if (_debug != 0)  printInfo(Helix);
      strcpy(banner,"refineHelixParameters-after-doLinearFitPhiZ");
      rc1 = refineHelixParameters(Helix, 0, _indicesTrkCandidate, banner, _debug);
      if (rc1 >=0 ){
	Helix._center.setX(Helix._sxyw.x0());
	Helix._center.setY(Helix._sxyw.y0());
	Helix._radius    = Helix._sxyw.radius();
	Helix._sxy.init(Helix._sxyw);
//2017-04-24 gianipez comment: if refineHelixParamters(...) removedone hit, we need to re-perform  the phi-z fit
//	doLinearFitPhiZ(Helix, 0, _indicesTrkCandidate, useInteligentWeight);
	if (_debug != 0)  printInfo(Helix);
      }
    }

    // 2014-11-09 gianipez changed the cleanup process. now it is faster and cleaner
    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition]: calling filterUsingPatternRecognition\n");

                                  // update hack data with last results
    Helix._diag.radius_14          = Helix._radius;
    Helix._diag.chi2_dof_circle_15 = Helix._sxy.chi2DofCircle();
    Helix._diag.z0_6               = Helix._fz0;
    Helix._diag.rdfdz_7            = Helix._dfdz*Helix._radius;
    Helix._diag.dfdz_8             = Helix._dfdz;

    //evaluate the radius without using the target center
    Helix._sxy.removePoint(0., 0., 1./900.);
    Helix._diag.dr                 = Helix._radius - Helix._sxy.radius();
    Helix._sxy.addPoint(0., 0., 1./900.);
    
    filterUsingPatternRecognition(Helix);

  PATTERN_RECOGNITION_END:;
//-----------------------------------------------------------------------------
// if running in the diagnostics mode, save state of the Xyzp (this is a deep copy)
//-----------------------------------------------------------------------------
    if (_debug != 0) printf("[CalHelixFinderAlg::doPatternRecognition]: END\n");
  }

//--------------------------------------------------------------------------------
// a straw man attempt to account for significantly different resolutions 
// along the wire and in the drift direction
//--------------------------------------------------------------------------------
  double  CalHelixFinderAlg::calculateWeight(const Hep3Vector& HitPos   ,
					     const Hep3Vector& StrawDir ,
					     const Hep3Vector& HelCenter,
					     double            Radius   ) {

    double    rs(2.5);   // straw radius, mm

    double x   = HitPos.x();
    double y   = HitPos.y();
    double dx  = x-HelCenter.x();
    double dy  = y-HelCenter.y();
    double dxn = dx*StrawDir.x()+dy*StrawDir.y();

    double costh2 = dxn*dxn/(dx*dx+dy*dy);
    double sinth2 = 1-costh2;

    double e2     = _ew*_ew*sinth2+rs*rs*costh2;
    double wt     = 1./e2;
                                                    // scale the weight for having chi2/ndof distribution peaking at 1
    wt *= _weightXY;

    return wt;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  double  CalHelixFinderAlg::calculatePhiWeight(const Hep3Vector&  HitPos   ,
						const Hep3Vector&  StrawDir ,
						const Hep3Vector&  HelCenter,
						double             Radius   ,
						int                Print    ,
						const char*        Banner   ) {
    double    rs(2.5);  // straw radius, mm

    double x  = HitPos.x();
    double y  = HitPos.y();
    double dx = x-HelCenter.x();
    double dy = y-HelCenter.y();

    double dxn = dx*StrawDir.x()+dy*StrawDir.y();

    double costh2  = dxn*dxn/(dx*dx+dy*dy);
    //    double costh  = (dx*StrawDir.x()+dy*StrawDir.y())/sqrt(dx*dx+dy*dy);
    double sinth2 = 1-costh2;

    double e2     = _ew*_ew*costh2+rs*rs*sinth2;
    double wt     = Radius*Radius/e2;
    wt           *= _weightZPhi;

    if (Print > 0) {
      double dr = calculateRadialDist(HitPos,HelCenter,Radius);
      printf("[CalHelixFinderAlg::%s] %9.3f %9.3f %10.5f %10.5f %10.5f %10.5f %12.5e %10.3f\n", 
	                       Banner, x, y, dx, dy, costh2, sinth2, e2, dr);
    }

    return wt;
  }

//--------------------------------------------------------------------------------
// calculate the radial distance of a straw hit from the helix prediction
//--------------------------------------------------------------------------------
  double  CalHelixFinderAlg::calculateRadialDist (const Hep3Vector& HitPos   ,
						  const Hep3Vector& HelCenter,
						  double            Radius   ) {
    double dx = HitPos.x()-HelCenter.x();
    double dy = HitPos.y()-HelCenter.y();
    double dr = sqrt(dx*dx+dy*dy)-Radius;

    return dr;
  }


//-----------------------------------------------------------------------------
  void   CalHelixFinderAlg::doWeightedCircleFit (::LsqSums4&     TrkSxy   ,
						 int             SeedIndex,
						 int*            IdVec    ,
						 Hep3Vector&     HelCenter,
						 double&         Radius   ,
						 double*         Weights  ,
						 int             Print    ,
						 const char*     Banner   ) {
    const Hep3Vector  *hitPos, *strawDir;
    double     wt;
    int        np = _xyzp.size();
//-----------------------------------------------------------------------------
// add calorimeter cluster with a position error of 10 mm => wt = 1/100
//-----------------------------------------------------------------------------
    TrkSxy.clear();
    TrkSxy.addPoint(fCaloX,fCaloY,1./100.);
//-------------------------------------------------------------------------------
// add stopping target center with a position error of 100 mm/sqrt(12) ~ 30mm => wt = 1/900
//-------------------------------------------------------------------------------
    TrkSxy.addPoint(0.,0.,1./900.);

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doWeightedCircleFit] BEGIN: x0 = %8.3f y0 = %8.3f radius = %8.3f chi2dof = %8.3f\n",
	     HelCenter.x(),HelCenter.y(),Radius,TrkSxy.chi2DofCircle());
      if (_debug > 10) {
	printf("[CalHelixFinderAlg::doWeightedCircleFit:LOOP] Index      X          Y         Z          wt        wireNx     wireNy\n");
      }
    }

    for (int i=SeedIndex; i<np; i++) {
      if ( IdVec[i] < 1)                                    continue;
      if ( _xyzp[i].isOutlier())                            continue;

      hitPos     = &_xyzp[i]._pos;
      strawDir   = &_xyzp[i]._sdir;

      wt         = calculateWeight(*hitPos,*strawDir,HelCenter,Radius);
      Weights[i] = wt;

      TrkSxy.addPoint(hitPos->x(),hitPos->y(),wt);
      if (_debug > 10) {
	printf("[CalHelixFinderAlg::doWeightedCircleFit:LOOP] %4i %10.3f %10.3f %10.3f %10.3e %10.4f %10.4f\n",
	       (int)_xyzp[i]._ind, hitPos->x(), hitPos->y(), hitPos->z(), wt, strawDir->x(), strawDir->y());
      }
    }
					// update helix info
    Radius  = TrkSxy.radius();
    HelCenter.setX(TrkSxy.x0());
    HelCenter.setY(TrkSxy.y0());

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doWeightedCircleFit:END] : npt = %3.0f  chi2dof = %8.3f x0 = %8.3f y0 = %8.3f radius = %8.3f\n",
	     TrkSxy.qn(),TrkSxy.chi2DofCircle(),HelCenter.x(),HelCenter.y(),Radius);
    }
  }


//-----------------------------------------------------------------------------
// this is a rather "primitive" definition of the worst hit, should do for now
//-----------------------------------------------------------------------------
  void    CalHelixFinderAlg::searchWorstHitWeightedCircleFit(int               SeedIndex,
							     int*              IdVec,
							     const Hep3Vector& HelCenter,
							     double&           Radius,
							     double*           Weights,
							     int&              Iworst,
							     double&           HitChi2Worst)
  {
    HitChi2Worst  = _hitChi2Max;
    Iworst        = -1;

    int        np = _xyzp.size();
    double     dr, hitChi2;

    for (int i=SeedIndex; i<np; i++) {
      if (IdVec[i] < 1)                                     continue;
      if (_xyzp[i].isOutlier())                             continue;

      dr      = calculateRadialDist(_xyzp[i]._pos,HelCenter,Radius);
      hitChi2 = dr*dr*Weights[i];
					// store info aout the radial residual
      if (SeedIndex == 0) {
	_distTrkCandidate[i] = hitChi2;
      }

      if (hitChi2 > HitChi2Worst) {
	HitChi2Worst = hitChi2;
	Iworst       = i;
      }
    }
  }

//--------------------------------------------------------------------------------
// IWorst is always defined
// returns the index of the hit which provides the highest contribute to the chi2
//--------------------------------------------------------------------------------
  void    CalHelixFinderAlg::cleanUpWeightedCircleFit(LsqSums4&     TrkSxy,
						      int           SeedIndex,
						      int*          IdVec,
						      double*       Weights,
						      int&          IWorst)
  {
    LsqSums4   sxy;
    double     chi2, chi2_min (-1.), x, y;
    int        np = _xyzp.size();

    IWorst     = -1;

    for (int i=SeedIndex; i<np; i++) {
					// ignore hits rejected by the helix search
      if ( IdVec[i] < 1)                                    continue;
      if (_xyzp[i].isOutlier())                             continue;

      sxy.init(TrkSxy);

      x  = _xyzp[i].x();
      y  = _xyzp[i].y();

      sxy.removePoint(x, y, Weights[i]);

      chi2  = sxy.chi2DofCircle();

      if ((chi2 < chi2_min) || (i == SeedIndex)) {
	chi2_min    = chi2;
	IWorst      = i;
      }
    }
  }

//-----------------------------------------------------------------------------
// use hits only, at this point the cluster is no longer needed 
//-----------------------------------------------------------------------------
  int CalHelixFinderAlg::refineHelixParameters(CalHelixFinderData& Trk,
					       int                 SeedIndex,
					       int*                IndexVec,
					       const char*         Banner,
					       int                 Print  ) {
    double     x, y, r;
    double     hitChi2Worst;

    ::LsqSums4 sxyw;
    int        pointsRemoved(0);
    int        iworst(-1);
    double     wtWorst;
    double     chi2, chi2_min;
    int        np = _xyzp.size();
					// use a local copy of IndexVec
    int       idVec[np];
    for (int i=0; i<np; ++i){
      idVec[i] = IndexVec[i];
    }

    Hep3Vector hitPos, strawDir, helCenter;

    double     weights[np];
    int        rc(0);                   // success-oriented initialization :)
					// initialize helix 
    r  = Trk._sxy.radius();
    helCenter.setX( Trk._sxy.x0());
    helCenter.setY( Trk._sxy.y0());

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::refineHelixParameters] BEGIN               x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f \n",
	     Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle());
      printf("[CalHelixFinderAlg::refineHelixParameters] i       X        Y        dx        dy         costh        sinth2         e2     radial-dist\n");
    }

    doWeightedCircleFit (sxyw,SeedIndex,idVec,helCenter,r,weights,Print,Banner);
//-----------------------------------------------------------------------------
// recalcute weights using most recent helix parameters
//-----------------------------------------------------------------------------
    doWeightedCircleFit (sxyw,SeedIndex,idVec,helCenter,r,weights,Print,Banner);

    searchWorstHitWeightedCircleFit(SeedIndex,idVec,helCenter,r,weights,iworst,hitChi2Worst);

    chi2     = sxyw.chi2DofCircle();
    chi2_min = chi2;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::refineHelixParameters] npt = %3.0f x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f iworst=%3i chi2Worst = %8.3f\n",
	     sxyw.qn(),sxyw.x0(),sxyw.y0(),sxyw.radius(),sxyw.chi2DofCircle(),iworst,hitChi2Worst);
    }

    if ((chi2 <= _chi2xyMax) && (hitChi2Worst <= _hitChi2Max)) goto F_END;
//-----------------------------------------------------------------------------
// one of teh chi2's is above the threshold, cleanup is needed
//-----------------------------------------------------------------------------
    if (_debug > 5) printf("[CalHelixFinderAlg::refineHelixParameters] : START CLEANUP\n");  
  NEXT_ITERATION:;

    cleanUpWeightedCircleFit(sxyw,SeedIndex,idVec,weights,iworst);

    if (iworst >= 0) {
      x       = _xyzp[iworst]._pos.x();
      y       = _xyzp[iworst]._pos.y();
      wtWorst = weights[iworst];
					// remove point from the track, this is why need to return weights
      sxyw.removePoint(x, y, wtWorst);

      if (_debug > 5) {
	printf("[CalHelixFinderAlg::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
	       sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle(), _chi2xyMax, iworst);
      }
					// mark point as outlier
      idVec[iworst] = 0;
					// update helix info
      r  = sxyw.radius();
      helCenter.setX( sxyw.x0());
      helCenter.setY( sxyw.y0());
					// now update helix
      doWeightedCircleFit (sxyw,SeedIndex,idVec,helCenter,r,weights,0,Banner);

					// update the chi2 value
      chi2_min = sxyw.chi2DofCircle();

      ++pointsRemoved;
    }
//-----------------------------------------------------------------------------
// recalculate the worst radial residual
//-----------------------------------------------------------------------------
  CHECK_RESIDUALS: ;
    searchWorstHitWeightedCircleFit(SeedIndex,idVec,helCenter,r,weights,iworst,hitChi2Worst);
//-----------------------------------------------------------------------------
// if a hit contributes chi2 > _hitCHi2Max, remove it and go back looking for the next such hit
//-----------------------------------------------------------------------------
    if (iworst >= 0) {
      x       = _xyzp[iworst]._pos.x();
      y       = _xyzp[iworst]._pos.y();
      wtWorst = weights[iworst];

      //remove point from the track
      sxyw.removePoint(x, y, wtWorst);
      if (_debug > 5) {
	printf("[CalHelixFinderAlg::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
	       sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle(), _chi2xyMax, iworst);

      }
      //mark point as outlier
      idVec[iworst] = 0;

      //update helix info
      r  = sxyw.radius();
      helCenter.setX(sxyw.x0());
      helCenter.setY(sxyw.y0());

      // now update helix
      doWeightedCircleFit (sxyw, SeedIndex, idVec,  helCenter,  r,  weights, 0, Banner);

      //update the chi2 value
      chi2_min = sxyw.chi2DofCircle();

      ++pointsRemoved;
                                                            goto CHECK_RESIDUALS;
    }

    if ((chi2_min >= _chi2xyMax) && (iworst >= 0)) {
//-----------------------------------------------------------------------------
// still bad chi2, repeat the cleanup cycle
//-----------------------------------------------------------------------------
      if (sxyw.qn() > 10.)                                  goto NEXT_ITERATION;
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
	     sxyw.qn(), sxyw.x0(), sxyw.y0(), sxyw.radius(), sxyw.chi2DofCircle(), pointsRemoved);
    }

    if (rc >= 0) {
//-----------------------------------------------------------------------------
// update circle parameters
//-----------------------------------------------------------------------------
      Trk._sxyw.init(sxyw);
      Trk._cw.setX(Trk._sxyw.x0());
      Trk._cw.setY(Trk._sxyw.y0());
      Trk._rw    = Trk._sxyw.radius();
      Trk._chi2w = Trk._sxyw.chi2DofCircle();

      for (int i=0; i<np; ++i) {
	IndexVec[i] = idVec[i];
      }
    }

    return rc;
  }


//-----------------------------------------------------------------------------
// The previous helix search may have  thrown away point which instead have
// small radial residual, so this function is devoted for rescueing these
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::rescueHits(CalHelixFinderData& Helix          ,
				     int                 SeedIndex      ,
				     int*                IndexVec       ,
				     int                 UsePhiResiduals) {

    const char  banner[] = "rescueHits";
    double      wt, e2, x, y, r;
    double      phiwt(-9999.);

    Hep3Vector  hitPos, strawDir, helCenter, hel_pred(0);

    int         np = _xyzp.size();
    double      weights[np];

    double      dfdz, phi0, dphi, dphiChi2(0.0), phi_pred;

    ::LsqSums4  sxy;
    int         n_added_points(0);
    int         ibest(-1);
    double      wtBest, phiwtBest;
    double      chi2, chi2_min, dr, hitChi2, drChi2;

    double      x_pred(0), y_pred(0), weight_hel(0), x_resid2(0), y_resid2(0);
    //set  dfdz and phi0
    dfdz = Helix._dfdz;
    phi0 = Helix._fz0;

    //update helix info
    r  = Helix._sxy.radius();
    helCenter.setX( Helix._sxy.x0());
    helCenter.setY( Helix._sxy.y0());

    doWeightedCircleFit (Helix._sxy, SeedIndex, IndexVec,  helCenter,  r,  weights);

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::%s:BEGIN] x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f phi0 = %9.6f dfdz = %9.6f chi2 = %8.3f\n",
	     banner,
	     Helix._sxy.x0(), Helix._sxy.y0(), Helix._sxy.radius(), Helix._sxy.chi2DofCircle(),
	     phi0, dfdz, Helix._srphi.chi2DofLine());
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
      chi2_min    = _hitChi2Max;//2.*_hitChi2Max;
    } else{
      chi2_min    = _hitChi2Max;
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
      wt = calculateWeight    (hitPos,strawDir,helCenter,r);

      drChi2  = (dr*dr)*wt;

      // 2015-03-25 G.Pezzu added the request of "distance" also for the phi

      if ((UsePhiResiduals == 1) && (_phiCorrectedDefined)) {
	phi_pred = hitPos.z()*dfdz + phi0;
	dphi     = phi_pred - _phiCorrected[i];
	phiwt    = calculatePhiWeight(hitPos, strawDir, helCenter, r, 0, banner);
	dphiChi2 = dphi*dphi*phiwt;

	// calculate distance from predicted point
	x         = hitPos.x();
	y         = hitPos.y();
	x_pred    = helCenter.x() + r*cos(phi_pred);
	y_pred    = helCenter.y() + r*sin(phi_pred);
	hel_pred.setX(x_pred);
	hel_pred.setY(y_pred);
	weight_hel = calculateWeight(hitPos, strawDir, hel_pred, r)*_weight3D;
	x_resid2  = (x - x_pred)*(x - x_pred);
	y_resid2  = (y - y_pred)*(y - y_pred);
	hitChi2   = (x_resid2 + y_resid2)*weight_hel;//drChi2 + dphiChi2;
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
	Helix._srphi.addPoint(hitPos.z(), _phiCorrected[i], phiwt);

	if (Helix._sxy.chi2DofCircle() < _chi2xyMax){
	  if (UsePhiResiduals == 1){
	    if (Helix._srphi.chi2DofLine() < _chi2zphiMax){
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

	Helix._sxy.removePoint  (x, y, wt);
	Helix._srphi.removePoint(hitPos.z(), _phiCorrected[i], phiwt);

      }
    }

    if (ibest >= 0){
      x  = _xyzp[ibest]._pos.x();
      y  = _xyzp[ibest]._pos.y();
                                       //add point from the track
      Helix._sxy.addPoint(x, y, wtBest);

      if (UsePhiResiduals == 1){
      	Helix._srphi.addPoint(_xyzp[ibest]._pos.z(), _phiCorrected[ibest], phiwtBest);
	dfdz  = Helix._srphi.dfdz(); 
	phi0  = Helix._srphi.phi0();
      }

      if (_debug > 5) {
	printf("[CalHelixFinderAlg::%s:PT2] x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %6.3f chi2Maxxy = %6.3f index point added = %i straw-id = %6i hitChi2 = %6.3f x = %8.3f y = %8.3f z = %9.3f\n",
	       banner,
	       Helix._sxy.x0(), Helix._sxy.y0(), Helix._sxy.radius(), Helix._sxy.chi2DofCircle(), _chi2xyMax, ibest,
	       _xyzp[ibest]._strawhit->strawIndex().asInt(), chi2_min,
	       x, y, _xyzp[ibest]._pos.z());
      }
					// mark point as active
      IndexVec[ibest] = 1;

      r  = Helix._sxy.radius();
      helCenter.setX( Helix._sxy.x0());
      helCenter.setY( Helix._sxy.y0());
                                        // update helix
      doWeightedCircleFit (Helix._sxy, SeedIndex, IndexVec,  helCenter,  r,  weights);
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
      dr      = calculateRadialDist(hitPos,helCenter,r);
      hitChi2 = dr*dr/e2;
					// store residual
      if (SeedIndex == 0){
	_distTrkCandidate[i] = hitChi2;
      }
    }
//-----------------------------------------------------------------------------
// update circle parameters
//-----------------------------------------------------------------------------
    Helix._center.setX(Helix._sxy.x0());
    Helix._center.setY(Helix._sxy.y0());
    Helix._radius  = Helix._sxy.radius();
    Helix._chi2    = Helix._sxy.chi2DofCircle();

  F_END:;
    if (_debug > 5 ) {
      printf("[CalHelixFinderAlg::%s:END] N(added) = %i chi2 = %5.5f\n",banner,n_added_points,Helix._chi2);
    }

    if (SeedIndex ==0){
      Helix._diag.n_rescued_points_16 = n_added_points;
    }
  }



//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::findTrack(int                 SeedIndex     ,
				    CalHelixFinderData& Helix         ,
				    bool                UseDefaultDfDz,
				    int                 UseMPVDfDz    ) {
//-----------------------------------------------------------------------------
// points used in the pattern recongition:
// ----------------------------------------------------------------
// center : center of the helix, center.z = z(hit with SeedIndex)
// p1     : center of the stopping target
// p2     : point in the position SeedIndex on the vector Xyzp
// p3     : postion of the EMC cluster
//-----------------------------------------------------------------------------
    double     radius, phi0, dx, dy, phi, Chi2;
    Hep3Vector center, shPos, hePos;
//-----------------------------------------------------------------------------
// mode0GoodPoints: number of points belonging to a trajectory when dfdz is not
//                  re-calculated using the function calculateDfDz()
// mode1GoodPoints: number of points belonging to a trajectory when dfdz is
//                  re-computed using calculateDfDz()
//-----------------------------------------------------------------------------
    int  Mode(0), mode0GoodPoints(0), mode1GoodPoints(0), rescuedPoints(0);
//-----------------------------------------------------------------------------
// initialize a vector to store the straw hit candidate for dfdz recalculation
// create a temporary array for storing indices of the good point which belong
// to a track candidate
//-----------------------------------------------------------------------------
    int np = _xyzp.size();

    int hitIndex     [np];  // index in _xyzp of the i-th hit of the candidate
    int markIndexList[np];

    double     dzList[np], distList[np];

    for(int i=0; i<np; ++i){
      markIndexList[i] =     0;
      dzList[i]        = -9999;
      distList[i]      = -9999;
    }

    markIndexList[SeedIndex] = 1;
    hitIndex[0]              = SeedIndex;
//---------------------------------------------------------------------
// define constrains on the z coordinate of the strawhit candidate for re-calculating dfdz
// If the candidate and the seeding strawhit are too close, the dfdz calculated could be 
// affected by several effects which lead to a wrong estimate
// We are asking that the candidate straw must be at a distance along
// the z axes greater than tollMin and less than tollMax.
// These parameters still need to be optimized
//-----------------------------------------------------------------------------
    double tollMin(100.) ; // , tollMax(500.);

    int    goodPoint(-1);               // index of the strawhit candidate for dfdz and helix parameters recalculation
					// parameters used to calculate the strawhit position residuals
    double weight(1.), wtarget(0.1);
    double deltaZ(0.), dist(0.), dist2(0.);
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
    double dfdz_end, phi0_end, radius_end;

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
    CalHelixPoint* seedHit = &_xyzp[SeedIndex];
    CalHelixPoint* lastHit = seedHit;

    Hep3Vector p1(0.,0.,0.);	           // target, z(ST) = 5971. - 10200. is not used
    Hep3Vector p2(seedHit->_pos);          // seed hit
    Hep3Vector p3(fCaloX,fCaloY,fCaloZ);   // cluster
    
    calculateTrackParameters(p1,p2,p3,center,radius,phi0,dfdz);
    
    double     tollMax = 2.*M_PI/dfdz;
//------------------------------------------------------------------------------
// helix parameters, in particular, phi0, are defined at Z=p2.z()
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE 
//------------------------------------------------------------------------------
    if (UseMPVDfDz ==1 ) dfdz = _hdfdz;			// _mpDfDz;

    int lastIndex = -9999;

    ::LsqSums4 sxy;
    ::LsqSums4 srphi;

    sxy.addPoint(p2.x(), p2.y(), 1.     );  // seed hit
    sxy.addPoint(p3.x(), p3.y(), 1.     );  // EMC cluster position
    sxy.addPoint(    0.,     0., wtarget);  // Target center in the transverse plane, with small weight

    std::string name("CalHelixFinderAlg::findTrack");

    hitIndex[0] = SeedIndex;
    int NPoints = 1;    // nhits, associated with the track, sxy has NPoints+2 or NPoints+1 points

    double  z_phi0 = p2.z();
    
    for (int i=SeedIndex+1; i<np; i++) {
      CalHelixPoint* hit = &_xyzp[i];
      if (hit->isOutlier())                             continue;
      shPos  = hit->_pos;
      double hitz   = hit->z();
//-----------------------------------------------------------------------------
// dfdz = tanLambda/radius; phi0 is the last found hit phi
//-----------------------------------------------------------------------------
      deltaZ = hitz-lastHit->z(); // distance from the last found hit
      phi    = phi0 + deltaZ*dfdz;                 

      hePos.set(center.x()+radius*cos(phi),center.y()+radius*sin(phi),hitz);

					// residuals in XY
      dx     = hePos.x() - hit->x();
      dy     = hePos.y() - hit->y();
      dist2  = dx*dx + dy*dy;
      dist   = std::sqrt(dist2);

      if( _debug > 10){
	if( i==SeedIndex+1) {
	  printf("[%s:LOOP]  findTrack() starts with helix parameters derived from these points \n",name.data());
	  printf("[%s:LOOP]   point  type      X         Y         Z       xyzp-index \n", name.data());
	  printf("[%s:LOOP] ----------------------------------------------------------\n", name.data());
	  printf("[%s:LOOP] seeding        %9.3f %9.3f %9.3f %5i\n",name.data(),p2.x(),p2.y(),p2.z(),SeedIndex );
	  printf("[%s:LOOP] candidate      %9.3f %9.3f %9.3f %5i\n",name.data(),p1.x(),p1.y(),p1.z(),lastIndex);
	  printf("[%s:LOOP] emc cluster    %9.3f %9.3f %9.3f %5i\n",name.data(),p3.x(),p3.y(),p3.z(),	 -1);
	  printf("[%s:LOOP]----------------------------------------------------------------------------------------------------------------------------------------\n",name.data());
	  printf("[%s:LOOP]   i      Z        xi       yi       xp       yp    dXYpred  dXYseed   dZseed    X0       Y0        R        phi      dfdz    chi2 added\n",name.data());
	  printf("[%s:LOOP]----------------------------------------------------------------------------------------------------------------------------------------\n",name.data());
	}
					      // dist betw the straw hit and the seed hit in the transverse plane

	double dx  = std::fabs(seedHit->x() - hit->x());
	double dy  = std::fabs(seedHit->y() - hit->y());
	double dxy = std::sqrt(dx*dx+dy*dy);
	double chi2   = sxy.chi2DofCircle();

	printf("[%s:LOOP] %3i %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f %8.5f %8.3f",
	       name.data(),i,hit->z(),hit->x(),hit->y(),hePos.x(),hePos.y(),dist,dxy,deltaZ,center.x(),center.y(),radius,phi,dfdz,chi2) ;
      }
//-----------------------------------------------------------------------------
// dxy_max: running search window accounts for the finite extrapolation accuracy
//-----------------------------------------------------------------------------
      double dxy_max = _distPatRec + _dfdzErr*deltaZ;
      if (dist <= dxy_max) {
	hitIndex[NPoints] = i;
	NPoints++;
//-----------------------------------------------------------------------------
// Mode = 0: helix parameters evaluated using stopping_target+fitst_hit(seed)+cluster
//           dphi/dz is fixed and so far set to the most proable value for CE
//-----------------------------------------------------------------------------
	sxy.addPoint(hit->x(),hit->y(),weight);
	if (Mode == 1) {
					// dfdz has already been evaluated, update XY part of the helix
	  center.setX(sxy.x0());
	  center.setY(sxy.y0());
	  radius = sxy.radius();
	}

	phi0    = atan2(hit->y()-center.y(),hit->x()-center.x());  // *DOUBLE_CHECK*
	z_phi0  = hitz;			                           // *DOUBLE_CHECK*
	lastHit = hit; 

	markIndexList[i] = 1;		// mark point as good
	dzList       [i] = deltaZ;      // Z-distance from the last point
	distList     [i] = dist;        // distance from prediction in XY

	if      ( Mode == 0                    ) ++mode0GoodPoints;
	else if ((Mode == 1) && (i<= lastIndex)) ++mode1GoodPoints;

	double dzFromSeed = lastHit->z() - seedHit->z();              // expected to be positive (non-negative)
	if ((dzFromSeed > tollMin) && (dzFromSeed < tollMax)) {
//-----------------------------------------------------------------------------
// goodPoint - index of the first hit separated in Z from the seed by > 10 cm
// after such found, the target center is removed and the circle parameters 
// recalculated using the cluster, the seed hit and the 'goodPoint' hit
// an additional requirement is that at the recalculation time there are 3 or more 
// hits found in total 
//-----------------------------------------------------------------------------
	  if (removeTarget) goodPoint = i;
	}
      }

      if (_debug > 10) printf (" %3i\n",markIndexList[i]);

      if ((goodPoint >= 0) && (goodPoint != lastIndex) && (NPoints >= 2)) {
//-----------------------------------------------------------------------------
// the first point separated from the seed one by more than 10 cm has been found
// recalculate helix parameters: for XY part use accumulated sxy sums
// replace stopping target with the hit
//-----------------------------------------------------------------------------
//	sxy.removePoint(0.,0.,wtarget);
	p1 = _xyzp[goodPoint]._pos;

	center.setX(sxy.x0());
	center.setY(sxy.y0());
	radius = sxy.radius();
//-----------------------------------------------------------------------------
// now calculate more accuratelly the value of dfdz using just the two strawhit positions
// change in the circle parameterization changes the phi0 value
//-----------------------------------------------------------------------------
	phi0 = atan2(hit->y()-center.y(),hit->x()-center.x());

	if (UseMPVDfDz == 0) {
	  //	  calculateDfDz(phi_0,phi_1,p2.z(),p1.z(),dfdz);
	  calculateDphiDz_2(hitIndex,NPoints,center.x(),center.y(),dfdz);
	}
	else if (UseMPVDfDz ==1) {
	  dfdz = _hdfdz;
	}

	if (_debug > 10) {
	  printf("[%s:DEF2]   i     Z       X        Y   type\n", name.data());
	  printf("[%s:DEF2] ----------------------------------------------------\n", name.data());
	  printf("[%s:DEF2] %3i %9.3f %8.3f %8.3f seed \n", name.data(),SeedIndex,seedHit->z(),seedHit->x(),seedHit->y());
	  printf("[%s:DEF2] %3i %9.3f %8.3f %8.3f seed \n", name.data(),goodPoint,hit->z()    ,hit->x()    ,hit->y()    );
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
	  p1.set(0.,0.,0.);
	  dfdz = _mpDfDz;
	}
	else {
//-----------------------------------------------------------------------------
// dPhi/Dz makes sense, exclude the stopping target and change the search mode
//-----------------------------------------------------------------------------
	  removeTarget = false;
	  Mode         = 1;
	  lastIndex    = goodPoint;
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
    CalHelixFinderData tmp1HelFitRes(Helix);

    tmp1HelFitRes._sxy.init(sxy);
    tmp1HelFitRes._radius = sxy.radius();
    tmp1HelFitRes._center.set(sxy.x0(), sxy.y0(), 0.0);

    radius_end = sxy.radius();

    CalHelixFinderData tmp2HelFitRes(tmp1HelFitRes);

    tmp2HelFitRes._dfdz   = dfdz;
    tmp2HelFitRes._fz0    = phi0 - dfdz*z_phi0; // *DOUBLE_CHECK*
    int rc = refineHelixParameters(tmp1HelFitRes, SeedIndex, markIndexList);
//-----------------------------------------------------------------------------
// if weighted XY fit didn't converge, there is nothing else one can do, return
//-----------------------------------------------------------------------------
    if (rc < 0) return;

    tmp2HelFitRes._center.set(tmp1HelFitRes._cw.x(), tmp1HelFitRes._cw.y(), 0.0);
    tmp2HelFitRes._radius = tmp1HelFitRes._rw;
    radius_end            = tmp1HelFitRes._rw;
    sxy.init(tmp1HelFitRes._sxyw); // *FIXME*  dont really need this, just rename tmp2HelFitRes and use it
					// doWeightedCircleFit still adds the ST and the cluster
    Chi2    = sxy.chi2DofCircle();
    NPoints = sxy.qn()-2;               //  *FIXME*  in principle, the fit can remove ST as well as the cluster
					// diagnostics
    radiusRes[1] = tmp2HelFitRes._radius;
//-----------------------------------------------------------------------------
// 2015-01-22 G. Pezzullo and P. Murat; update the dfdz value using all hits
//-----------------------------------------------------------------------------
    int rs = findDfDz(tmp2HelFitRes, SeedIndex, markIndexList);
    if (rs ==1 ) {
      tmp2HelFitRes._dfdz = _hdfdz;
      tmp2HelFitRes._fz0  = _hphi0;
					// fill diag vector
      dfdzRes[1]  = _hdfdz;
      dphi0Res[1] = _hphi0;
    }
//-----------------------------------------------------------------------------
// 2015-01-23 G. Pezzu and P. Murat: when it fails, doLinearFitPhiZ returns negative value
//                                   in this case, use the previous value for dfdz and phi0
//-----------------------------------------------------------------------------
    bool rcPhiZ = doLinearFitPhiZ(tmp2HelFitRes, SeedIndex, markIndexList);

    if (rcPhiZ) {
      dfdz_end   = tmp2HelFitRes._dfdz;
      phi0_end   = tmp2HelFitRes._fz0;
      srphi.init(tmp2HelFitRes._srphi); // *FIXME* dont really need this, just rename tmp2HelFitRes and use it
					// fill diagnostic vector
      dfdzRes [2] = tmp2HelFitRes._dfdz;
      dphi0Res[2] = tmp2HelFitRes._fz0;

      NPoints = 0;
      for (int i=SeedIndex; i<np; ++i) {
	if ( markIndexList[i] > 0)
	  ++NPoints;
      }
    }
    else {
      dfdz_end = _hdfdz;
      phi0_end = _hphi0;
    }

    if (_debug > 10) {
      printf("[%s] strawhit type     X        Y        Z     index     \n", name.data());
      printf("[%s] ----------------------------------------------------\n", name.data());
      printf("[%s]    seeding   %9.3f %9.3f %9.3f   %i  \n", name.data(),p2.x(), p2.y(), p2.z(), SeedIndex);
      printf("[%s]   candidate  %9.3f %9.3f %9.3f   %i  \n", name.data(),p1.x(), p1.y(), p1.z(), goodPoint);
      printf("[%s]  emc cluster %9.3f %9.3f %9.3f       \n", name.data(),p3.x(), p3.y(), p3.z());
      printf("[%s] NPoints = %i x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.6fdfdz = %5.6f chi2 = %5.3f \n", 
	     name.data(),NPoints,sxy.x0(), sxy.y0(), radius_end, phi0_end, dfdz_end , sxy.chi2DofCircle());
    }

    if (mode1GoodPoints > 0) rescuedPoints = mode1GoodPoints - mode0GoodPoints ;
    else                     rescuedPoints = -1;
//-----------------------------------------------------------------------------
// all points in the time window checked
//----------------------------------------------------------------------
    if (NPoints > 2) Chi2 = sxy.chi2DofCircle();
    else             Chi2 = -1;

    if (_debug > 5) {
      printf("[%s:END] SeedIndex:%3i Mode= %i UseDefaultDfDz:%1i UseMPVDfDz:%1i NPoints= %3i Chi2=%5.3f dfdz= %5.5f \n",
	     name.data(),SeedIndex,Mode,UseDefaultDfDz,UseMPVDfDz,NPoints,Chi2,dfdz_end);
    }

    // can execution really come here with Mode == 0 ? - YES! not sure, why
    //    if ((Mode == 1) || UseDefaultDfDz) {
      if (( NPoints >  _goodPointsTrkCandidate) ||
	  ((NPoints == _goodPointsTrkCandidate) && (Chi2< _chi2TrkCandidate))) {
//-----------------------------------------------------------------------------
// found candidate is better, than the best previous one
//-----------------------------------------------------------------------------
	//	fSeedIndex              = SeedIndex;
	_goodPointsTrkCandidate = NPoints;
	_chi2TrkCandidate       = Chi2;
//-----------------------------------------------------------------------------
// reset the vector holding the informations about:
// -> hit belonging to the track candidate
// -> distance in the X-Y plane from the prediction
// -> distance from the seeding hit along the z-axes
//-----------------------------------------------------------------------------
	for (int i=0; i<kMaxNHits; ++i) {
	  _indicesTrkCandidate[i] = -9999;
	  _distTrkCandidate   [i] = -9999;
	  _dzTrkCandidate     [i] = -9999;
	}

	for (int i=SeedIndex; i<np; ++i){
	  _indicesTrkCandidate[i] = markIndexList[i];
	  _distTrkCandidate   [i] = distList     [i];
	  _dzTrkCandidate     [i] = dzList       [i];
	}

	Helix._nPoints   = NPoints;
	Helix._seedIndex = SeedIndex;
	Helix._candIndex = lastIndex;
//-----------------------------------------------------------------------------
// helix parameters and LSQ sums, phi is defined at z=0
//-----------------------------------------------------------------------------
	Helix._center.set(sxy.x0(),sxy.y0(), 0.0);
	Helix._radius = radius_end;	// _radius;
	Helix._fz0    = phi0_end;	// _phi0;
	Helix._dfdz   = dfdz_end;	// _dfdz;
	Helix._sxy.init(sxy);
	Helix._srphi.init(srphi);
//-----------------------------------------------------------------------------
// if needed, fill diagnostics information for histogramming
//-----------------------------------------------------------------------------
	if (_diag > 0) {
	  int loopId(1);
	  if (UseDefaultDfDz == 0) {
	    if (UseMPVDfDz) loopId = 2;
	    else            loopId = 0;
	  }

	  Helix._diag.loopId_4           = loopId;
	  Helix._diag.radius_5           = Helix._radius;
	  Helix._diag.n_rescued_points_9 = rescuedPoints;

	  double dz = p1.z() - p2.z();

	  Helix._diag.dz_10              = (Mode == 1) ? dz : -1.;
	  Helix._diag.n_active_11        = NPoints;
	  Helix._diag.chi2_dof_circle_12 = sxy.chi2DofCircle();
	  Helix._diag.chi2_dof_line_13   = srphi.chi2DofLine();

	  Helix._diag.dfdzres_17 = dfdzRes[0];
	  Helix._diag.dfdzres_18 = dfdzRes[1];
	  Helix._diag.dfdzres_19 = dfdzRes[2];

	  Helix._diag.dr_20 = radiusRes[0];
	  Helix._diag.dr_21 = radiusRes[1];

	  Helix._diag.dphi0res_22 = dphi0Res[0];
	  Helix._diag.dphi0res_23 = dphi0Res[1];
	  Helix._diag.dphi0res_24 = dphi0Res[2];

	  int j=0;
	  for (int i=SeedIndex; i<np; ++i){
	    if (_indicesTrkCandidate[i] != 1) continue;
	    if (j < Helix.maxIndex()) {
	      Helix._diag.dist[j] = _distTrkCandidate[i];
	      Helix._diag.dz  [j] = _dzTrkCandidate  [i];
	      ++j;
	    }
	    else {
	      printf("ERROR in CalHelixFinderAlg::findTrack : index out limits. IGNORE; \n");
	    }
	  }
	}
      }
      //    }
  }

//-----------------------------------------------------------------------------
// helix parameters are defined at Z=p2.z, Phi0 corresponds to p2
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::calculateTrackParameters(const Hep3Vector&   p1       ,
						   const Hep3Vector&   p2       ,
						   const Hep3Vector&   p3       ,
						   Hep3Vector&         Center   ,
						   double&             Radius   ,
						   double&             Phi0     ,
						   double&             DfDz23) {
    Center.setZ(p2.z());

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
    Center.setX(x0);
    double y0 = m*x0 + c;   //(c - t) * m / (m - k) + t;
    Center.setY(y0);
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

    Phi0        = std::atan2(dy2,dx2);
//-----------------------------------------------------------------------------
// this assumes that the helix is right-handed, *FIXME*
// make sure that we are lookign for a particle which makes the number of turns
// close to the expected 
//-----------------------------------------------------------------------------
    double dphi32 = std::atan2(dy3,dx3) - Phi0;
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
    if (_debug > 5) {
//-----------------------------------------------------------------------------
// in debug mode also want to print the helix parameters, calculate them
//-----------------------------------------------------------------------------
      double d0     = sqrt(x0*x0+y0*y0)-Radius;
      double phi00  = atan2(y0,x0)+M_PI/2;   // for negatively charged particle
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
  void  CalHelixFinderAlg::calculateDphiDz_2(const int* HitIndex, int NHits, double X0, double Y0, double& DphiDz) {
    LsqSums2 sphiz;

    double phi, phi0, phiCl;
    
    phiCl = atan2(fCaloY-Y0,fCaloX-X0);
    if (phiCl < 0) phiCl += 2*M_PI;

    for (int i=0; i<NHits; i++) {
      int ind         = HitIndex[i];
      Hep3Vector* pos = &_xyzp[ind]._pos;
      phi = atan2(pos->y()-Y0,pos->x()-X0);
      if (i == 0) phi0 = phi;

      if (phi-phi0 >  M_PI) phi -= 2*M_PI;
      if (phi-phi0 < -M_PI) phi += 2*M_PI;

      sphiz.addPoint(pos->z(),phi);

      if (_debug > 10) printf("[CalHelixFinderAlg::calculateDphiDz_2:LOOP] i,ind,phi,z=%3i %3i %8.5f %9.3f\n",i,ind,pos->z(),phi);
    }
//-----------------------------------------------------------------------------
// define straight line phi = phi0+dPhi/Dz*z , where phi0 = phi(z=0)
//-----------------------------------------------------------------------------
    phi0   = sphiz.yMean();
    DphiDz = sphiz.dydx();

    if (_debug > 5) printf("[CalHelixFinderAlg::calculateDphiDz_2:END] fCaloZ,phiCl: %9.3f %8.5f phi0,DphiDz = %9.5f %9.5f \n",fCaloZ,phiCl,phi0,DphiDz);
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

    helx->print("from CalHelixFinderAlg::plotZPhi");
  }

//-----------------------------------------------------------------------------
// 'TrkIndex' = 0 : mark all hits as active
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::saveResults(std::vector<CalHelixPoint>& Xyzp ,
				      CalHelixFinderData&         Helix,
				      int                         Index) {
    _results[Index]._xyzp  = Xyzp;
    _results[Index]._helix = Helix;
  }


//-----------------------------------------------------------------------------
// 2017-01-25 P.Murat: print
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::printXYZP(CalHelixFinderData& Helix) {
    int n = _xyzp.size();

    printf("[CalHelixFinderAlg::printXYZP]-----------------------------------------------------------------------------------------\n");
    printf("[CalHelixFinderAlg::printXYZP]     i index strawID   flag    Used     X         Y         Z      _debug: %5i nhits: %5i\n",_debug,n);
    printf("[CalHelixFinderAlg::printXYZP]-----------------------------------------------------------------------------------------\n");
    
    const vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();

    for(int i=0; i<n; ++i) {
      CalHelixPoint* pt = &_xyzp[i];

      int loc    = shIndices[i];	 // index in shcol of i-th timecluster hit
    
      const StrawHit& sh          = Helix.shcol()->at(loc);

      printf("[CalHelixFinderAlg::printXYZP] %5i %5i %5i   %08x   %2i %9.3f %9.3f %9.3f \n",
	     i, loc,  sh.strawIndex().asInt(), *((int*) &pt->_flag), pt->use(), pt->_pos.x(), pt->_pos.y(), pt->_pos.z());
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::resolve2PiAmbiguity(const CLHEP::Hep3Vector& Center, double DfDz, double Phi0){
  
    const CLHEP::Hep3Vector* pos;
    double                   z, phi, phi_ref, dphi;

    int np = _xyzp.size();

    for(int i=0; i<np; ++i){
      pos = &_xyzp[i]._pos;
      z   = pos->z();

      phi = CLHEP::Hep3Vector(*pos - Center).phi();
      phi = TVector2::Phi_0_2pi(phi);
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
      _phiCorrected[i] = phi;
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
    //    _dfdz           = -9999.;
//-----------------------------------------------------------------------------
// quality paramters used for doing comparison between several track candidates
//-----------------------------------------------------------------------------
    _goodPointsTrkCandidate = 0;
    _chi2TrkCandidate       = 1e10;
//-----------------------------------------------------------------------------
// vector which holds indices of the strawhit far from the predicted position
//-----------------------------------------------------------------------------
    for(int i=0; i<kMaxNHits; ++i){
      _indicesTrkCandidate[i] = -9999;
      _distTrkCandidate   [i] = -9999;
      _dzTrkCandidate     [i] = -9999;
      _phiCorrected       [i] = -9999;
    }
  }

}
