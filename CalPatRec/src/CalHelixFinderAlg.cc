///////////////////////////////////////////////////////////////////////////////
// Helix fit to straw hits
//
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
#include "Offline/CalPatRec/inc/CalHelixFinderAlg.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "art_root_io/TFileService.h"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/BFieldGeom/inc/BFieldConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
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

#include "Offline/CalPatRec/inc/CalHelixFinderAlg.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"

using CLHEP::HepVector;
using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;

namespace mu2e {

//   // comparison functor for ordering points
//   struct radcomp : public std::binary_function<VALERR, VALERR, bool> {
//     bool operator()(VALERR const& r1, VALERR const& r2) { return r1._val < r2._val; }
//   };

  // comparison functor for sorting by z
  struct zcomp {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1._pos.z() < p2._pos.z(); }
  };

  // comparison functor for sorting byuniquePanel ID
  struct panelcomp {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1.strawId().uniquePanel() < p2.strawId().uniquePanel(); }
  };

//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::defineHelixParams(CalHelixFinderData& Helix) const {

    static const float pi(M_PI), halfpi(pi/2.0);

    HepVector pvec(5,0), perr(5,0);

    // the helix fit introduces a radial bias due to an asymmetry in the detector (more phase space for
    // noise hits outside the circle than inside.  correct for it.
    float radius = Helix._radius; //  + _rbias;
    // omega is the inverse transverse radius of the particle's circular motion.
    // It is signed by the particle angular momentum about the cirle center.
    // This CANNOT be deduced geometrically, so must be supplied as an ad-hoc assumption
    float amsign = Helix._helicity == Helicity::poshel ? 1 : -1;//copysign(1.0,-Helix._tpart.charge()*bz());
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
    float dphi = deltaPhi(Helix._fz0+amsign*halfpi,pvec[HelixTraj::phi0Index]);
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

  CalHelixFinderAlg::CalHelixFinderAlg(const Config& config) :
    _diag(config.diag()),
    _debug(config.debug()),
    _debug2(config.debug2()),
    _hsel(config.hsel()),
    _bkgsel(config.bkgsel()),
    _maxHitEnergy(config.maxHitEnergy()),
    _minNHits(config.minNHits()),
    _absMpDfDz(config.absMpDfDz()),
    _initDfDz(config.initDfDz()),
    _dzOverHelPitchCut(config.dzOverHelPitchCut()),
    _maxDfDz(config.maxDfDz()),
    _minDfDz(config.minDfDz()),
    _sigmaPhi(config.sigmaPhi()),
    _weightXY(config.weightXY()),
    _targetcon(config.targetcon()),
    _weightZPhi(config.weightZPhi()),
    _weight3D(config.weight3D()),
    _maxXDPhi(config.maxXDPhi()),
    _maxPanelToHelixDPhi(config.maxPanelToHelixDPhi()),
    _distPatRec(config.distPatRec()),
    _minDeltaNShPatRec(config.minDeltaNShPatRec()),
    _mindist(config.mindist()),
    _pmin(config.pmin()),
    _pmax(config.pmax()),
    _tdmin(config.tdmin()),
    _tdmax(config.tdmax()),
    _xyweights(config.xyweights()),
    _zweights(config.zweights()),
    _filter(config.filter()),
    _plotall(config.plotall()),
    _usetarget(config.usetarget()),
    _maxZTripletSearch(config.maxZTripletSearch()),
    _bz(0.),
    _nHitsMaxPerPanel(config.nHitsMaxPerPanel()),
    _hitChi2Max(config.hitChi2Max()),
    _chi2xyMax(config.chi2xyMax()),
    _chi2zphiMax(config.chi2zphiMax()),
    _chi2hel3DMax(config.chi2hel3DMax()),
    _dfdzErr(config.dfdzErr()),
    _maxNHitsRatio(config.maxNHitsRatio()),
    _procAllTCs(config.procAllTCs()),
    _slopeRatioLimit(config.slopeRatioLimit()){

    float minarea(config.minArea());
    _minarea2    = minarea*minarea;

    std::vector<std::string> bitnames;
    bitnames.push_back("Outlier");
    bitnames.push_back("OtherBackground");
    //mu2e::ComboHit::_useflag = StrawHitFlag(bitnames);
  }

//-----------------------------------------------------------------------------
  CalHelixFinderAlg::~CalHelixFinderAlg() {
  }

//-----------------------------------------------------------------------------
  float CalHelixFinderAlg::bz() const {
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

  void CalHelixFinderAlg::setCaloCluster(CalHelixFinderData& Helix) {
    //check presence of a cluster
    const CaloCluster* cl = Helix._timeCluster->caloCluster().get();
    if (cl == NULL){
      const std::vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();

      if((shIndices.size() !=0) && (_procAllTCs)){
        float  zMax(-1e10);
        size_t zMaxIndex(0);
        for (size_t i=0; i<shIndices.size(); ++i){
          StrawHitIndex   ll   = shIndices[i];
          float           hitZ = Helix.chcol()->at(ll).pos().z();
          if (hitZ > zMax){
            zMax      = hitZ;
            zMaxIndex = i;
          }
        }
        StrawHitIndex   loc = shIndices[zMaxIndex];
        const ComboHit& ch  = Helix.chcol()->at(loc);

        fCaloTime = ch.correctedTime();
        fCaloX    = ch.pos().x();
        fCaloY    = ch.pos().y();
        fCaloZ    = ch.pos().z();
      }
      else{
        fCaloTime = -9999.;
        fCaloX    = -9999.;
        fCaloY    = -9999.;
        fCaloZ    = -9999.;
      }

      return;
    }
    //fill the calorimeter cluster info
    Hep3Vector  gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskID(),cl->cog3Vector());
    Hep3Vector  tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
    fCaloTime        = cl->time();
    fCaloX           = tpos.x();
    fCaloY           = tpos.y();
    float     offset = _calorimeter->caloInfo().getDouble("diskCaseZLength")/2. + (_calorimeter->caloInfo().getDouble("BPPipeZOffset") + _calorimeter->caloInfo().getDouble("BPHoleZLength")+ _calorimeter->caloInfo().getDouble("FEEZLength"))/2. - _calorimeter->caloInfo().getDouble("FPCarbonZLength") - _calorimeter->caloInfo().getDouble("FPFoamZLength");
    fCaloZ           = tpos.z()-offset;
  }


//-----------------------------------------------------------------------------
//2016-12-26 gianipez added the following function to make CalHelixFinderAlg compatible with TimeCluster obj
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::findHelix(CalHelixFinderData& Helix) {

    // fTimeCluster = TimePeak;
    //check presence of a cluster
    // const CaloCluster* cl = TimePeak->caloCluster().get();
    // if (cl == NULL)   return false;

    // //fill the calorimeter cluster info
    // Hep3Vector         gpos = _calorimeter->geomUtil().diskToMu2e(cl->diskID(),cl->cog3Vector());
    // Hep3Vector         tpos = _calorimeter->geomUtil().mu2eToTracker(gpos);
    // fCaloTime = cl->time();
    // fCaloX    = tpos.x();
    // fCaloY    = tpos.y();
    // fCaloZ    = tpos.z();
//-----------------------------------------------------------------------------
//  compute the allowed radial range for this fit
//-----------------------------------------------------------------------------
    float pb = fabs((CLHEP::c_light*1e-3)/(bz()*Helix._tpart.charge()));
    _rmin = _pmin/(pb*sqrt(1.0+_tdmax*_tdmax));
    _rmax = _pmax/(pb*sqrt(1.0+_tdmin*_tdmin));
//-----------------------------------------------------------------------------
//  particle charge, field, and direction affect the pitch range
//-----------------------------------------------------------------------------
    _dfdzsign = Helix._helicity == Helicity::poshel ? 1 : -1;// copysign(1.0,-Helix._tpart.charge()*Helix._fdir.dzdt()*bz());

    _smax     = _dfdzsign/(_rmax*_tdmax);
    _smin     = _dfdzsign/(_rmin*_tdmin);

    _mpDfDz   = _dfdzsign*_absMpDfDz;
//-----------------------------------------------------------------------------
// call down
//-----------------------------------------------------------------------------
    bool retval = fitHelix(Helix);

    return retval;
  }

//-----------------------------------------------------------------------------
// called internally; in the diagnostics mode save several states of _xyzp
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::fitHelix(CalHelixFinderData& Helix) {
    bool retval(false);
                                        // initialize internal array of hits, print if requested
    //    fillXYZP(Helix); //2019-01-18: gianipez moved this into the CalHelixFinder_module to exploit the helicity loop-search
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
    else if ((Helix._nZPhiSh < _minNHits) || (Helix._szphi.chi2DofLine() > _chi2zphiMax) ||
             (fabs(Helix._szphi.dfdz()) < _minDfDz) || (fabs(Helix._szphi.dfdz()) > _maxDfDz)) {
      Helix._fit = TrkErrCode(TrkErrCode::fail,4); // phi-z reconstruction failure
    }
    else {
      //-----------------------------------------------------------------------------
      // success, form output
      //-----------------------------------------------------------------------------
      Helix._goodhits.clear();

      PanelZ_t*      panelz(0);
      FaceZ_t*       facez(0);

      for (int f=StrawId::_ntotalfaces-1; f>=0; --f){
        facez     = &Helix._oTracker[f];
        for (int p=0; p<FaceZ_t::kNPanels; ++p){
          panelz = &facez->panelZs[p];
          int  nhits = panelz->nChHits();
          for (int i=0; i<nhits; ++i){
            int index = panelz->idChBegin + i;// facez->evalUniqueHitIndex(f,p,i);
            if (Helix._hitsUsed[index] != 1)                continue;
            Helix._goodhits.push_back(index); // HitInfo_t(f,p,index));
          }
        }
      }

      defineHelixParams(Helix);
      retval = true;
    }

    return retval;
  }

//-----------------------------------------------------------------------------
  int CalHelixFinderAlg::findDfDz(CalHelixFinderData& Helix,
                                  HitInfo_t          SeedIndex,
                                   // int *IndexVec,
                                  int                 Diag_flag) {
    //    return findDfDz_1(Helix, SeedIndex, IndexVec, Diag_flag);
    return findDfDz_2(Helix, SeedIndex, Diag_flag);
  }

//----------------------------------------------------------------------------------------
// 2015-01-13  calculate track DphiDz using histogrammed distribution of the dfdz residuals
//----------------------------------------------------------------------------------------
  int CalHelixFinderAlg::findDfDz_1(CalHelixFinderData& Helix, HitInfo_t SeedIndex, int  Diag_flag) {

//     float phi, phi_ref(-1e10), z, z_ref, dphi, dz, dzOverHelPitch;

//     CLHEP::Hep3Vector* center = &Helix._center;
//     CLHEP::Hep3Vector pos_ref;

//     _hDfDzRes->Reset();
//     _hPhi0Res->Reset();
//                                         // 2015 - 03 -30 G. Pezzu changed the value of tollMax.
//                                         // using the initial value of dfdz we can set it more accuratelly:
//                                         // tollMax = half-helix-step = Pi / dfdz
//     float tollMin(100.);
// //-----------------------------------------------------------------------------
// // 2017-09-26 gianipez fixed a bug: in case the Helix phi-z fit didn't converge yet,
// // Helix._dfdz is set to -1e6, so we need to make a check here!
// // this is a tempOrary fix that doesn't take into account the particle helicity. FIX ME!
// //-----------------------------------------------------------------------------
//     float helix_dfdz(_mpDfDz);
//     // 2017-11-14 gianipez: findDfDz shoudl use the dfdz value obtained only from the linearFit
//     if (Helix._szphi.qn() >= 10) helix_dfdz = Helix._szphi.dfdz();
//     //    if (Helix._dfdz > 0) helix_dfdz =  Helix._dfdz;
//     float tollMax = 2.*M_PI / helix_dfdz;

//     if (_debug > 5) {
//       printf("[CalHelixFinderAlg::findDfDz:BEGIN] x0 = %9.3f y0 = %9.3f Helix._radius = %9.3f",
//              center->x(), center->y(), Helix._radius);
//       printf("helix_dfdz = %9.6f Helix._nStrawHits = %3i tollMax = %8.6f\n",helix_dfdz, Helix._nStrawHits, tollMax);
//     }

//     int       nstations, nhits[30];
//     float    phiVec[30], zVec[30], weight(0), weight_cl(0);
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
//     for (int p=SeedIndex.Panel; p<FaceZ_t::kNTotalPanels; ++p){
//       panelz = &Helix._oTracker[p];
//       int  nhitsPerPanel  = panelz->fNHits;
//       int  seedPanelIndex(0);
//       if (nhitsPerPanel == 0)                                                            continue;
//       if (p==SeedIndex.Panel) seedPanelIndex = SeedIndex.PanelHitIndex;

//       for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
//         mu2e::ComboHit* hit = &panelz->_chHitsToProcess.at(i);
//         int index = p*FaceZ_t::kNMaxHitsPerPanel + i;
//         if (Helix._hitsUsed[index] != 1)                     continue;

//         int ist = hit->_straw->id().getStation();                   // station number
//         phi     = polyAtan2(hit->y()-center->y(),hit->x()-center->x()); // atan2 returns its result in [-pi,pi], convert to [0,2pi]
//         if (phi < 0) phi += 2*M_PI;
//         zVec  [ist] += hit->z();
//         //-----------------------------------------------------------------------------
//         // make sure there all hits within the station get close values of phi, although a
//         // common 2pi ambiguity is still unresolved
//         //-----------------------------------------------------------------------------
//         if (nhits[ist] == 0) phiVec[ist] = phi;
//         else {
//           while (phi-phiVec[ist] >  M_PI) phi -= 2*M_PI;
//           while (phi-phiVec[ist] < -M_PI) phi += 2*M_PI;

//           phiVec[ist] = (phiVec[ist]*nhits[ist]+phi)/(nhits[ist]+1);
//         }
//         nhits [ist] += 1;
//       }
//     }

//     for (int i=0; i<nstations; i++) {
//       if (nhits[i] > 0) {
//         zVec  [i] = zVec  [i]/nhits[i];
//       }
//     }

//     if (_debug >5) {
//       printf("[CalHelixFinderAlg::findDfDz] StationID  nhits       z        phi\n");
//       for (int i=0; i<nstations; i++) {
//         if (nhits[i] > 0) printf("[CalHelixFinderAlg::findDfDz] %5i %6i    %9.3f %8.5f\n", i,nhits[i],zVec[i],phiVec[i]);
//       }
//     }

//     int i0(-1), first_point(1);
//                                         // add the cluster phi
//     float zCl   = fCaloZ;
//     float phiCl = polyAtan2(fCaloY-center->y(),fCaloX-center->x());
//     if (phiCl < 0) phiCl += 2*M_PI;

//     for (int i=0; i<nstations; i++) {
//       if (nhits[i] == 0)                                    continue;
//                                         // i0: fist station with hits
//       if (first_point) {
//         i0          = i;
//         first_point = 0;
//       }

//       phi_ref = phiVec[i];
//       z_ref   = zVec  [i];

//       for(int j=i+1; j<nstations; ++j) {
//         if (nhits[j] == 0)                                  continue;
//         phi = phiVec[j];
//         z   = zVec  [j];
//         dz  = z - z_ref;

//         dzOverHelPitch = dz/tollMax - int(dz/tollMax);
//         weight         = nhits[i] + nhits[j];

//         //        if ((phi_ref > -9999) && (dzOverHelPitch < _dzOverHelPitchCut) && (dz > tollMin)) {
//         if (dz > tollMin) {
//           dphi = phi-phi_ref;
//            while (dphi >  M_PI) dphi -= 2*M_PI;
//            while (dphi < -M_PI) dphi += 2*M_PI;
// //-----------------------------------------------------------------------------
// // add 2*PI to take into account the fact we are in the second loop
// // FIX ME: what to do if we are in the third loop?
// //-----------------------------------------------------------------------------
//           if (dz > tollMax) dphi += 2*M_PI*int(dz/tollMax);

//           float dphidz = dphi/dz;
//           while (dphidz < 0.) {
//             dphi    = dphi+2.*M_PI;
//             dphidz  = dphi/dz;
//           }
//           _hDfDzRes->Fill(dphidz, weight);

//           float tmpphi0 = phi_ref - dphidz*z_ref;
//           tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

//           if (_debug > 5) {
//             printf("[CalHelixFinderAlg::findDfDz:1] z_ref: %9.3f z: %9.3f dz: %9.3f",z_ref,z,dz);
//             printf(" phi_ref: %9.5f phi: %9.5f dphi: %9.5f dz/HelPitch: %10.3f dphi/dz: %9.5f phi0 = %9.6f\n",
//                    phi_ref,phi,dphi, dzOverHelPitch, dphidz, tmpphi0);
//           }
// //-----------------------------------------------------------------------------
// // in case dfdz is out of limits, set tmpphi0 as negative
// //-----------------------------------------------------------------------------
//           if ((dphidz < _minDfDz) || (dphidz >  _maxDfDz)) tmpphi0 = -1;
//           _hPhi0Res->Fill(tmpphi0, weight);
//         }
//       }
// //-----------------------------------------------------------------------------
// // include the calorimeter cluster phi
// //-----------------------------------------------------------------------------
//       dz             = zCl - z_ref;
//       dzOverHelPitch = dz/tollMax - int(dz/tollMax);
//       weight_cl      =  nhits[i];

//       //      if ((dzOverHelPitch < _dzOverHelPitchCut) && (dz > tollMin)) {
//       if (dz > tollMin) {
//         dphi  = phiCl - phi_ref;
//         dphi  = TVector2::Phi_0_2pi(dphi);
// //-----------------------------------------------------------------------------
// // add 2 pi for taking into account the fact we are in the second loop
// // *FIX ME*: what if we are in the third loop?
// //-----------------------------------------------------------------------------
//         if (dz > tollMax) dphi += 2*M_PI*int(dz/tollMax);

//         float dphidz = dphi/dz;
//         while (dphidz < 0.) {
//           dphi   += 2.*M_PI;
//           dphidz = dphi/dz;
//         }

//         float tmpphi0 = phi_ref - dphidz*z_ref;
//         tmpphi0        = TVector2::Phi_0_2pi(tmpphi0);

//         if (_debug > 5){
//           printf("[CalHelixFinderAlg::findDfDz:2] z_ref: %9.3f z: %9.3f dz: %9.3f",z_ref,zCl,dz);
//           printf(" phi_ref: %9.5f phi: %9.5f dphi: %9.5f dz/HelPitch: %10.3f dphi/dz: %9.5f phi0 = %9.6f\n",
//                  phi_ref,phiCl,dphi, dzOverHelPitch, dphidz, tmpphi0);
//         }

//         if (dzOverHelPitch < _dzOverHelPitchCut ) {
//           _hDfDzRes->Fill(dphidz, weight_cl);
//           if ((dphidz < _minDfDz) || (dphidz >  _maxDfDz)) tmpphi0 = -1;
//           _hPhi0Res->Fill(tmpphi0, weight_cl);
//         }
//       }
//     }
// //-----------------------------------------------------------------------------
// // 2015 - 04- 02 G. Pezzu changed the way the maximum is searched
// // since sometimes a 2pi ambiguity creates two peaks in the histogram
// // we want to use the second, because it is the correct one
// //-----------------------------------------------------------------------------
//     float  maxContent = _hDfDzRes->GetMaximum() - 0.001;
//     int      maxBin    = _hDfDzRes->FindLastBinAbove(maxContent);//GetMaximumBin();
//     _hdfdz             = _hDfDzRes->GetBinCenter(maxBin);//_hDfDzRes->GetMean();
//     float dfdzmean    = _hDfDzRes->GetMean();
//     int    nentries    = _hDfDzRes->GetEntries();
//     int    overflows   = _hDfDzRes->GetBinContent(0)  + _hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1);

//     maxContent         = _hPhi0Res->GetMaximum() - 0.001;
//     maxBin             = _hPhi0Res->FindLastBinAbove(maxContent);//GetMaximumBin();

//     float mpvphi0     = _hPhi0Res->GetBinCenter(maxBin); //_hPhi0Res->GetMean();
//     float menaphi0    = _hPhi0Res->GetMean();
//     int    nentriesphi = _hPhi0Res->GetEntries();

//     _hphi0 = mpvphi0;  // 2018-01-05: *FLOAT_CHECK*

//     if (_debug > 5) {
//       printf("[CalHelixFinderAlg::findDfDz:DFDZ] nent: %3i mpvDfDz: %9.6f meanDphiDz: %9.6f under: %3.0f over: %3.0f ENTRIES:",
//              nentries, _hdfdz, dfdzmean,
//              _hDfDzRes->GetBinContent(0),_hDfDzRes->GetBinContent(_hDfDzRes->GetNbinsX()+1)
//              );
//       for (int i=0; i<_hDfDzRes->GetNbinsX(); i++) {
//         printf(" %3.0f",_hDfDzRes->GetBinContent(i+1));
//       }
//       printf("\n");

//       printf("[CalHelixFinderAlg::findDfDz:PHI0] nent: %3i mpvPhi0: %9.6f meanPhi0  : %9.6f under: %3.0f over: %3.0f ENTRIES:",
//              nentriesphi, mpvphi0,  menaphi0,
//              _hPhi0Res->GetBinContent(0),_hPhi0Res->GetBinContent(_hPhi0Res->GetNbinsX()+1)
//              );
//       for (int i=0; i<_hPhi0Res->GetNbinsX(); i++) {
//         printf(" %3.0f",_hPhi0Res->GetBinContent(i+1));
//       }
//       printf("\n");
//     }
// //-----------------------------------------------------------------------------
// // Part 2: try to perform a more accurate estimate - straight line fit
// //-----------------------------------------------------------------------------
//     float z0, phi0, dphidz, pred;

//     z0     = 0.    ;
//     phi0   = _hphi0;
//     dphidz = _hdfdz;
//     //    _sdfdz = -1;

//     if (_debug > 5) {
//       float tmpphi0=phi0+dphidz*z0;
//       printf("[CalHelixFinderAlg::findDfDz:PART2] phi0 = %9.6f dfdz = %9.6f\n", tmpphi0, dphidz);
//     }
// //--------------------------------------------------------------------------------
// // 2015-03-25 G. Pezzu changed the way the 2PI ambiguity is resolved
// //--------------------------------------------------------------------------------
//     LsqSums4 szphi;

//     weight = 1./(_sigmaPhi*_sigmaPhi);

//     float xdphi, zLast(z0), zdist;

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

//     float errCl    = 2./30.;
//     float weightCl = 1./(errCl*errCl);

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
//              0, zCl, dz, phiCl, phiCl, szphi.dfdz(), dphi, xdphi, dphidz);
//     }

//     // 2015-07-06 Gianipez added the following line for avoiding infinite loops
//     if ( i0 < 0 ) goto NEXT_STEP;

//     for (int i=i0; i<nstations; i++) {
//       if (nhits[i] > 0) {
//         z    = zVec[i];
//         dz   = z-z0;
//         pred = phi0 + dz*dphidz;
//         phi  = phiVec[i];
//         dphi = phi - pred;//pred - phi;

//         while (dphi > M_PI){
//           phi -= 2*M_PI;//+= 2*M_PI;
//           dphi = phi - pred;//pred - phi;
//         }
//         while (dphi < -M_PI){
//           phi += 2*M_PI;//-= 2*M_PI;
//           dphi = phi - pred;//pred - phi;
//         }

//         xdphi = fabs(dphi)/_sigmaPhi;

//         if (xdphi < 2.*_maxXDPhi){
//           szphi.addPoint(z, phi, weight);

//           zdist = z - zLast;

//           if ( (szphi.qn() >= 3.) && (zdist > 500.)){
//             z0     = 0.;
//             phi0   = szphi.phi0();
//             dphidz = szphi.dfdz();
//           }
//         }

//         if (_debug > 10) {
//           float tmpDfDz = szphi.dfdz();//, Helix._szphi.chi2DofLine());
//           printf("[CalHelixFinderAlg::findDfDz:LOOP] %3i %9.3f %9.3f %9.5f %9.5f %9.6f %9.6f %9.6f %9.6f\n",
//                  i, z, dz, phiVec[i], phi, tmpDfDz, dphi, xdphi, dphidz);
//         }
//       }
//     }

//   NEXT_STEP:;

//     if (szphi.qn() >= 3.) {
//       _hdfdz = szphi.dfdz();                // sigxy/sigxx;
//       _hphi0 = szphi.phi0();                // ymean - xmean*sigxy/sigxx;
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
//              _hphi0, -1.);
//       printf(" FIT: szphi.dfdz() = %9.5f szphi.phi0() = %9.6f chi2 = %9.3f qn = %6.0f\n", szphi.dfdz(),
//              szphi.phi0(), szphi.chi2DofLine(), szphi.qn());
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
  int CalHelixFinderAlg::findDfDz_2(CalHelixFinderData& Helix, HitInfo_t SeedIndex, int  Diag_flag) {

    float phi, phi_ref(-1e10), z_ref, dphi, dz;

    //    float hist[20], minX(0), maxX(0.01), stepX(0.0005), nbinsX(20); // make it 20 bins
    float hist[50], minX(0), maxX(0.025), stepX(0.0005), nbinsX(50); // make it 20 bins: gianipez test 2019-09-23

    XYZVectorF* center = &Helix._center;
    XYZVectorF  pos_ref;
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
    float    phiVec[30], zVec[30], weight(0);

    // np        = _xyzp.size();
    int       nPlanesPerStation(2);
    nstations = StrawId::_nplanes/nPlanesPerStation;//_tracker->nStations();

    for (int i=0; i<nstations+1; i++) {
      phiVec[i] = 0;
      zVec  [i] = 0;
      nhits [i] = 0;
    }

    for (int i=0; i<nbinsX; i++) hist[i] = 0;
//-----------------------------------------------------------------------------
// calorimeter cluster - point number nstations+1
//-----------------------------------------------------------------------------

    float zCl     = fCaloZ;
    float phiCl   = polyAtan2(fCaloY-center->y(),fCaloX-center->x());
    if (phiCl < 0) phiCl += 2*M_PI;

    phiVec[nstations] = phiCl;
    zVec  [nstations] = zCl;
    nhits [nstations] = 1;
    //-----------------------------------------------------------------------------
    // Step 1: for each station with track candidate hits, calculate average phi per station
    //-----------------------------------------------------------------------------
    PanelZ_t* panelz(0);
    FaceZ_t*  facez(0);

    for (int f=SeedIndex.face; f>=0; --f){
      facez     = &Helix._oTracker[f];

      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;

      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];
        int  nhitsPerPanel  = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhitsPerPanel == 0)                                      continue;
        //        if ( (f == SeedIndex.face) && (p==SeedIndex.panel) ) seedPanelIndex = SeedIndex.panelHitIndex;
        if ( (f == SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
          int           index = panelz->idChBegin + i;
          mu2e::ComboHit* hit = &Helix._chHitsToProcess[index];

          // int index =  facez->evalUniqueHitIndex(f,p,i);
          if (Helix._hitsUsed[index] != 1 )                         continue;

          int ist = hit->strawId().station();//_straw->id().getStation();                   // station number
          phi     = polyAtan2(hit->pos().y()-center->y(),hit->pos().x()-center->x()); // atan2 returns its result in [-pi,pi], convert to [0,2pi]
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
      }//end panels loop
    }//end face loop

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

        if(std::fabs(dz) < 10e-10)         continue;
        if((dz < 0.) && (!Helix._timeCluster->hasCaloCluster())) {
          dz = -dz;
        }

        float dphidz =dphi/dz*_dfdzsign; //HERE

        weight = nhits[i] + nhits[j];
//-----------------------------------------------------------------------------
// calculate N potential choices for 2*PI ambiguity resolution
//-----------------------------------------------------------------------------
        int n(0), nmax(0), nmin(0), nchoices = 0;

        float x = dphidz + n*2*M_PI/dz;
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
          float x = dphidz + n*2*M_PI/dz;
          int bin = (x-minX)/stepX;
          hist[bin] += weight;
        }
      }
    }
//-----------------------------------------------------------------------------
// the 'histogram' is filled, find a peak
//-----------------------------------------------------------------------------
    int ixmax = int(maxX/stepX);

    float swmax(0), sw, xmp(0);
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
    else                         _hdfdz = xmp*_dfdzsign;
//-----------------------------------------------------------------------------
// last step - determine phi0 = phi(z=0)
//-----------------------------------------------------------------------------
    float phi0(0), sdphi(0);
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
        float dphi_min = dphi;

        int n= 0;
        while (1) {
          n += 1;
          float dphi = phiVec[i]+2*M_PI*n-(phi0+zVec[i]*_hdfdz);
          if (fabs(dphi) < fabs(dphi_min)) dphi_min = dphi;
          else break;
        }

        n=0;
        while (1) {
          n -= 1;
          float dphi = phiVec[i]+2*M_PI*n-(phi0+zVec[i]*_hdfdz);
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
//--------------------------------------------------------------------------------
// function to count the nukber of used hits in units of ComboHits and Strawhits
//--------------------------------------------------------------------------------
  void   CalHelixFinderAlg::countUsedHits           (CalHelixFinderData& Helix,
                                                     HitInfo_t           SeedIndex,
                                                     int&                NComboHits,
                                                     int&                NPoints){
    FaceZ_t*        facez(0);
    PanelZ_t*       panelz(0);
    mu2e::ComboHit* hit(0);
    int             index(0);

    for (int f=SeedIndex.face; f>=0;  --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];
        int  nhitsPerPanel  = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhitsPerPanel == 0)                                 continue;
        if ( (f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
          index = panelz->idChBegin + i;
          hit   = &Helix._chHitsToProcess[index];
          if (Helix._hitsUsed[index] > 0 )  {
            ++NComboHits;
            NPoints += hit->nStrawHits();
          }
        }
      }//endl panels loop
    }//end face loop
  }
//-----------------------------------------------------------------------------
// subroutine used in the ::doLinearFitPhiZ function
//
//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::findGoodFaceHitInFitPhiZ (CalHelixFinderData& Helix,
                                                    PhiZFitInfo_t&      PhiZInfo,
                                                    HitInfo_t&          GoodFaceHit,
                                                    float&              FaceHitChi2){
      PanelZ_t*      panelz(0);
      int            firstPanel(0);
      float         dphi, err, xdphi, phi_ref;

      FaceZ_t*       facez     = &Helix._oTracker[PhiZInfo.faceId];
      float         z         = Helix._zFace[PhiZInfo.faceId];
      XYZVectorF         helCenter = Helix._center;
      float         radius    = Helix._radius;
      float         weight    = PhiZInfo.weight;

      if (Helix._sxy.qn() > 0) {
        helCenter = XYZVectorF( Helix._sxy.x0(), Helix._sxy.y0(), 0);
        radius    = Helix._sxy.radius();
      }

      if (PhiZInfo.faceId == PhiZInfo.seedIndex.face) firstPanel = PhiZInfo.seedIndex.panel;

      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];
        int       nhits  = panelz->nChHits();
        int       seedPanelIndex(0);
        if (nhits == 0)                                                                    continue;
        if ( (PhiZInfo.faceId == PhiZInfo.seedIndex.face )  &&
             (p               == PhiZInfo.seedIndex.panel)  &&
             (PhiZInfo.seedIndex.panelHitIndex >=0) ) seedPanelIndex = PhiZInfo.seedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhits; ++i){
          int             index = panelz->idChBegin + i;
          mu2e::ComboHit* hit   = &Helix._chHitsToProcess[index];
          PhiZInfo.dz           = z - PhiZInfo.zlast;
          // predicted value of phi
          phi_ref  = z*PhiZInfo.dfdz + PhiZInfo.phi0;
          // resolve 2PI ambiguity
          resolve2PiAmbiguity(hit, helCenter, phi_ref, dphi);

          dphi     = fabs(dphi);
          err      = _sigmaPhi;

          if (PhiZInfo.useInteligentWeight == 1){
            weight  = calculatePhiWeight(*hit, helCenter, radius, 0, PhiZInfo.banner);
            err     = 1./sqrt(weight);
          }
          hit->_zphiWeight = weight;

          xdphi = dphi/err;

          if (_debug > 10) {
            printf("[CalHelixFinderAlg::doLinearFitPhiZ:LOOP] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
                   *((int*) &hit->_flag), Helix._hitsUsed[index],
                   hit->strawId().straw(), i,
                   z, hit->_hphi, dphi, xdphi, PhiZInfo.zlast, PhiZInfo.dz,
                   PhiZInfo.dfdz, Helix._szphi.dfdz(), Helix._szphi.chi2DofLine());
          }

          if (Helix._hitsUsed[index] != 1)                     continue;

          if ( (PhiZInfo.doCleanUp == 1) && (xdphi > _maxXDPhi) ) {
            Helix._hitsUsed[index] = 0;
            ++PhiZInfo.nPointsRemoved;
            continue;
          }

          if (xdphi < FaceHitChi2){
            FaceHitChi2               = xdphi;
            GoodFaceHit.face          = PhiZInfo.faceId;
            GoodFaceHit.panel         = p;
            GoodFaceHit.panelHitIndex = index;
          }
        }

      }//end panels loop


  }

//-----------------------------------------------------------------------------
// subroutine used in the ::doLinearFitPhiZ function
//
//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::findWorstChi2HitInFitPhiZ (CalHelixFinderData& Helix,
                                                     PhiZFitInfo_t&      PhiZInfo,
                                                     HitInfo_t&          WorstFaceHit,
                                                     float&             HitChi2){
    PanelZ_t*      panelz(0);
    FaceZ_t*       facez(0);
    float         z, phi, weight(PhiZInfo.weight), chi2;
    ::LsqSums4     szphi;

    for (int f=PhiZInfo.seedIndex.face; f>=0; --f){
      facez     = &Helix._oTracker[f];
      z         = Helix._zFace[f];
      int            firstPanel(0);
      if (f == PhiZInfo.seedIndex.face) firstPanel = PhiZInfo.seedIndex.panel;
      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];
        int       nhits  = panelz->nChHits();
        int       seedPanelIndex(0);
        if (nhits == 0)                                         continue;
        if ( (f      == PhiZInfo.seedIndex.face )  &&
             (p      == PhiZInfo.seedIndex.panel)  &&
             (PhiZInfo.seedIndex.panelHitIndex >=0) ) seedPanelIndex = PhiZInfo.seedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhits; ++i){
          int             index = panelz->idChBegin + i;
          mu2e::ComboHit* hit   = &Helix._chHitsToProcess[index];
          if (Helix._hitsUsed[index] != 1)            continue;

          szphi = Helix._szphi;
          phi      = hit->_hphi;

          if (PhiZInfo.useInteligentWeight == 1){
            weight = hit->_zphiWeight;
          }

          szphi.removePoint(z, phi, weight);
          chi2 = szphi.chi2DofLine();
          if (chi2 < HitChi2) {
            WorstFaceHit.face          = f;
            WorstFaceHit.panel         = p;
            WorstFaceHit.panelHitIndex = index;
            HitChi2                    = chi2;
          }

        }//end panel-hits loop

      }//end panels loop
    }//end faces loop

  }

//-----------------------------------------------------------------------------
// subroutine used in the ::doLinearFitPhiZ function
//
//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::findWorstResidHitInFitPhiZ (CalHelixFinderData& Helix,
                                                      PhiZFitInfo_t&      PhiZInfo,
                                                      HitInfo_t&          WorstFaceHit,
                                                      float&             Resid){
    PanelZ_t*      panelz(0);
    FaceZ_t*       facez(0);
    float         z, phi, weight(PhiZInfo.weight), xdphi, dphi, err;

    for (int f=PhiZInfo.seedIndex.face; f>=0; --f){
      facez     = &Helix._oTracker[f];
      z         = Helix._zFace[f];
      int            firstPanel(0);
      if (f == PhiZInfo.seedIndex.face) firstPanel = PhiZInfo.seedIndex.panel;
      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];
        int       nhits  = panelz->nChHits();
        int       seedPanelIndex(0);
        if (nhits == 0)                                         continue;
        if ( (f      == PhiZInfo.seedIndex.face )  &&
             (p      == PhiZInfo.seedIndex.panel)  &&
             (PhiZInfo.seedIndex.panelHitIndex >=0) ) seedPanelIndex = PhiZInfo.seedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhits; ++i){
          int             index = panelz->idChBegin + i;
          mu2e::ComboHit* hit   = &Helix._chHitsToProcess[index];
          if (Helix._hitsUsed[index] != 1)            continue;
          phi      = z*Helix._szphi.dfdz() + Helix._szphi.phi0();
          dphi     = fabs(hit->_hphi - phi);
          err      = _sigmaPhi;

          if (PhiZInfo.useInteligentWeight == 1){
            weight = hit->_zphiWeight;//calculatePhiWeight(*hit,  helCenter, radius, 0, PhiZInfo.banner);
            err    = 1./sqrt(weight);
          }

          xdphi = dphi/err;

          if ( xdphi > Resid) {
            WorstFaceHit.face          = f;
            WorstFaceHit.panel         = p;
            WorstFaceHit.panelHitIndex = index;

            Resid                = xdphi;
          }

        }//end panel-hits loop

      }//end panels loop
    }//end faces loop

  }
//-----------------------------------------------------------------------------
// subroutine used in the ::doLinearFitPhiZ function
//
//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::addCaloClusterToFitPhiZ(CalHelixFinderData& Helix){
    float dfdz      = Helix._dfdz;
    float phi0      = Helix._fz0;
    XYZVectorF helCenter = Helix._center;
    if (Helix._sxy.qn() > 0) {
      helCenter = XYZVectorF( Helix._sxy.x0(), Helix._sxy.y0(), 0);
    }

    float zCl   = fCaloZ;
    float dx    = (fCaloX - helCenter.x());
    float dy    = (fCaloY - helCenter.y());
    float phiCl = polyAtan2(dy, dx);
    if (phiCl < 0) phiCl = phiCl + 2*M_PI;

    float deltaPhi = zCl*dfdz + phi0 - phiCl;
    while (deltaPhi > M_PI){
      phiCl   += 2*M_PI;
      deltaPhi = zCl*dfdz + phi0 - phiCl;
    }
    while (deltaPhi < -M_PI){
      phiCl   -= 2*M_PI;
      deltaPhi = zCl*dfdz + phi0 - phiCl;
    }

//check residual before adding to the LSqsum
    float xdphi  = fabs(deltaPhi)/_sigmaPhi;

// weight_cl of 10 corresponds to an uncertainty of 0.1 in phi(cluster),
// which means sigma(R-phi) of about 2-3 cm, about right
// hit weight is determined by the assumed phi error of _sigmaPhi=0.1
    if ( xdphi < 2.*_maxXDPhi ) {
      float weight_cl = 784.;//10.0;
      Helix._szphi.addPoint(zCl,phiCl,weight_cl);
      Helix._nZPhiSh += 1;
    }

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ] %08x %2i %6i %3i %12.5f %12.5f %10.5f %10.3f %10.3f %10.3f %10.5f %10.5f %5.3f\n",
             0, 1, 0, 0,  zCl, phiCl, deltaPhi, xdphi, 0., 0., dfdz, 0., 0.);
    }

  }

//-----------------------------------------------------------------------------
// function that performs the fit to a line in the phi-z space.
// the function performs also the resolution of the 2pi ambiguity
//-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::doLinearFitPhiZ(CalHelixFinderData& Helix    ,
                                          HitInfo_t           SeedIndex,
                                          int                 UseInteligentWeight,
                                          int                 DoCleanUp           ) {

    auto              hitsUsed  = Helix._hitsUsed;

    bool              success(false);
    int               nPointsRemoved(0);

    float            phi(0.0), z(0.0), weight(0.0);
    PanelZ_t*         panelz(0);
    FaceZ_t*          facez(0);
    mu2e::ComboHit*   hit(0);
//-----------------------------------------------------------------------------
// gianipez: procedure for aligning the phi vector
//-----------------------------------------------------------------------------
    ::LsqSums4        szphi;
    int               count(0), minNFitHits(5);
    float             chi2min, deltaPhi, dphi_max(0);

    HitInfo_t         iworst;//(-1,-1);

    const char        banner[200] = "doLinearFitPhiZ";
    if (UseInteligentWeight == 0){
      weight = 1./(_sigmaPhi*_sigmaPhi);
    }

    PhiZFitInfo_t     phiZInfo;
    phiZInfo.dfdz                 = Helix._dfdz;
    phiZInfo.phi0                 = Helix._fz0;
    phiZInfo.seedIndex            = SeedIndex;
    phiZInfo.weight                  = weight;
    phiZInfo.useInteligentWeight  = UseInteligentWeight;
    phiZInfo.dz                   = 0;
    phiZInfo.zlast                = 0;
    phiZInfo.nPointsRemoved       = 0;
    phiZInfo.doCleanUp            = DoCleanUp;
    phiZInfo.banner               = banner;

//reset the Lsqsum
    Helix._szphi.clear();
    Helix._nZPhiSh = 0;
//-----------------------------------------------------------------------------
// add the cluster phi to the LSq sum
//-----------------------------------------------------------------------------
    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ:BEGIN] phi0 = %10.6f dfdz = %10.6f chi2N = %10.3f DoCleanup = %i\n",
             Helix._fz0,  Helix._dfdz, 0.,DoCleanUp);
      printf("[CalHelixFinderAlg::doLinearFitPhiZ]    flag   A   shID   i       z         ");
      printf("    phi         dphi      xdphi      zlast        dz      dphidz  szphidfdz  chi2\n");
    }

    if(Helix._timeCluster->hasCaloCluster()){
      addCaloClusterToFitPhiZ(Helix);
    }

    count = 0;

    for (int f=SeedIndex.face; f>=0; --f){
      // reset the chi2 values calulated in the previous panel
      HitInfo_t      goodFaceHit;
      float          faceHitChi2(1e10);

      phiZInfo.faceId = f;
      findGoodFaceHitInFitPhiZ(Helix, phiZInfo, goodFaceHit, faceHitChi2);

      if (goodFaceHit.face < 0)                                  continue;
      hit    = &Helix._chHitsToProcess[goodFaceHit.panelHitIndex];
      z      = Helix._zFace[goodFaceHit.face];

      Helix._szphi.addPoint(z,hit->_hphi,hit->_zphiWeight);
      Helix._nZPhiSh    += hit->nStrawHits();

      ++count;

      if (count == 1) {//FIXME! investigate if it is needed or not
              phiZInfo.zlast = z;
        phiZInfo.dz    = 0.;
      }

      //update the dfdz and phi0 if...
      if ( (count>=minNFitHits) &&      //FIXME!
           (faceHitChi2 < 2.) &&
           ( (fabs(phiZInfo.dfdz - Helix._szphi.dfdz()) < 8.e-4) ) &&//  || //require that the new value of dfdz is
                                 //close to the starting one. update dfdz only if:
           ((Helix._szphi.dfdz()*_dfdzsign) > 0.) && //{                    // 1. the points browsed are more the half
           (-1.*phiZInfo.dz >=_mindist ) ){
        phiZInfo.dfdz  = Helix._szphi.dfdz();                     //    delta hits could have moved dfdz to negative value!
        phiZInfo.phi0  = Helix._szphi.phi0();                     // 2. and require dfdz to be positivie! scattered hits or
        phiZInfo.zlast = z;


      }
    }//end face loop
    _phiCorrectedDefined = 1;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doLinearFitPhiZ:BEFORE_CLEANUP] Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f points removed = %4i\n",
             Helix._szphi.phi0(),Helix._szphi.dfdz(), Helix._szphi.chi2DofLine(), nPointsRemoved);
    }
    //-----------------------------------------------------------------------------
    // perform a cleanup in RZ
    //-----------------------------------------------------------------------------
    if ( (DoCleanUp == 1) && (Helix._szphi.qn()>=minNFitHits)){
      if ( Helix._szphi.chi2DofLine() > _chi2zphiMax) {
      NEXT_ITERATION:;
        //reset the coordinates of the worst hit
        iworst.face          = -1;
        iworst.panel         = -1;
        iworst.panelHitIndex = -1;
        chi2min              = 1e10;

        findWorstChi2HitInFitPhiZ(Helix, phiZInfo, iworst, chi2min);

        if ((iworst.panel >= 0) && (Helix._nZPhiSh > _minNHits)) {
          hit     = &Helix._chHitsToProcess[iworst.panelHitIndex];
          Helix._hitsUsed[iworst.panelHitIndex] = 0;

          z   = Helix._zFace[iworst.face];
          phi = hit->_hphi;

          Helix._szphi.removePoint(z, phi, hit->_zphiWeight);
          Helix._nZPhiSh -= hit->nStrawHits();

          chi2min = Helix._szphi.chi2DofLine();
          if (_debug > 5) {
            printf("[CalHelixFinderAlg::doLinearFitPhiZ_removed:LOOP2] %6i %5.3f     %5.3f chi2 = %5.3f  \n", iworst.panelHitIndex, z, phi, chi2min);//FIXME! remove indexworst
          }
        }

      CHECK_RESIDUALS:;
        dphi_max    = _maxXDPhi;
        //reset the coordinates of the worst hit
        iworst.face          = -1;
        iworst.panel         = -1;
        iworst.panelHitIndex = -1;

        findWorstResidHitInFitPhiZ(Helix, phiZInfo, iworst, dphi_max);

        //remove the point
        if(iworst.panel>=0 && Helix._nZPhiSh > _minNHits){
          int   index = iworst.panelHitIndex;
          hit         = &Helix._chHitsToProcess[index];

          Helix._hitsUsed[index] = 0;

          z           = Helix._zFace[iworst.face];
          phi         = hit->_hphi;

          Helix._szphi.removePoint(z, phi, hit->_zphiWeight);
          Helix._nZPhiSh -= hit->nStrawHits();

          chi2min     = Helix._szphi.chi2DofLine();

          if (_debug > 5) {
            printf("[CalHelixFinderAlg::doLinearFitPhiZ:REMOVED] %6i %5.3f     %5.3f chi2 = %5.3f  \n", iworst.panelHitIndex, z, phi, chi2min);
          }
          goto CHECK_RESIDUALS;
        }

        if(Helix._nZPhiSh <= _minNHits) chi2min = Helix._szphi.chi2DofLine();

        if ( (chi2min >= _chi2zphiMax) ||
             (iworst.panel>=0 )) {

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
    if ( (Helix._szphi.qn() < minNFitHits) ||
         ((Helix._szphi.qn() >= minNFitHits) && (Helix._szphi.dfdz()*_dfdzsign < 0.)) ) {
      success = false;
    }
    else if (success) {                               // update helix results
      Helix._fz0  = Helix._szphi.phi0();
      Helix._dfdz = Helix._szphi.dfdz();
    }

    if ((SeedIndex.face == 0 ) && (SeedIndex.panel == 0) && (SeedIndex.panelHitIndex == 0) && (_diag > 0)) {
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
      Helix._diag.phi0_6           = Helix._szphi.phi0();
      Helix._diag.rdfdz_7          = Helix._szphi.dfdz()* Helix._diag.n_rescued_points_9;
      Helix._diag.dfdz_8           = Helix._szphi.dfdz();
      Helix._diag.chi2_dof_line_13 = Helix._szphi.chi2DofLine();

      if (success) {
        int h=0;

        for (int f=SeedIndex.face; f>=0; --f){
          facez     = &Helix._oTracker[f];
          int  firstPanel(0);
          if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
          z        = Helix._zFace[f];
          for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
            panelz = &facez->panelZs[p];
            int  nhitsPerPanel  = panelz->nChHits();
            int  seedPanelIndex(0);
            if (nhitsPerPanel == 0)                                                            continue; //Could the error be here?
            if ( (f == SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0) ) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

            for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
              int   index = panelz->idChBegin + i;
              hit = &Helix._chHitsToProcess[index];

              if (Helix._hitsUsed[index] != 1)                      continue;
              phi      = z* Helix._dfdz + Helix._fz0;
              deltaPhi = hit->_hphi - phi;

               if (h < Helix.maxIndex()) {
                Helix._diag.resid[h] = deltaPhi;
                ++h;
              }

              else {
                printf (" ERROR: too many hits. Ignore \n");
              }
            } //There seems to be an error here in where the brackets are placed.
          }
        }//end loop over the panels
      }//end loop pver the faces
    }

    if (_debug > 5) {

      printf("[CalHelixFinderAlg::doLinearFitPhiZ:END] retval = %d Helix: phi_0 = %5.3f dfdz = %5.5f chi2N = %5.3f szphi: phi_0 = %5.3f dfdz = %5.5f\n",
             success ? 1:0, Helix._szphi.phi0(),Helix._szphi.dfdz(), Helix._szphi.chi2DofLine(), szphi.phi0(), szphi.dfdz() );

      if (_debug > 10) {
        printf("[CalHelixFinderAlg::doLinearFitPhiZ:END2]    flag   A   shID       z             phi      deltaPhi\n");

        int lastFace  = SeedIndex.face;

        for (int f=StrawId::_ntotalfaces-1; f>= lastFace; --f){
          facez     = &Helix._oTracker[f];
          z         = Helix._zFace[f];
          int    lastPanel(0);
          if (f == SeedIndex.face) lastPanel = SeedIndex.panel;

          for (int p=FaceZ_t::kNPanels-1; p>=lastPanel;  --p){
            panelz = &facez->panelZs[p];
            int  nhitsPerPanel  = panelz->nChHits();
            if (nhitsPerPanel == 0)                                                            continue;
            if ( (f == SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) nhitsPerPanel = SeedIndex.panelHitIndex - panelz->idChBegin;

            for (int  i=nhitsPerPanel-1;i>=0; --i){
              int    index = panelz->idChBegin + i;
              hit      = &Helix._chHitsToProcess[index];
              phi      = z* Helix._dfdz + Helix._fz0;
              deltaPhi = hit->_hphi - phi;

              //              int index =  facez->evalUniqueHitIndex(f,p,i);//p*FaceZ_t::kNMaxHitsPerPanel + i;
              printf("[CalHelixFinderAlg::doLinearFitPhiZ:END2] %08x %2i %6i %12.5f %12.5f %12.5f\n",
                     *((int*) &hit->_flag), Helix._hitsUsed[index],
                     hit->strawId().straw()/*_strawhit->strawIndex().asInt()*/,  z, hit->_hphi, deltaPhi);
            }
          }//end panel loop
        }//end face loop
      }
    }

    if (!success) {
      Helix._hitsUsed = hitsUsed;
    }

    return success;
  }

//-----------------------------------------------------------------------------
// 12-09-2013 gianipez modified this procedure to avoid the doubling of the
// same stereohitposition
// points in filled array are ordered in Z coordinate
//-------------------------------------------------------------------------
  void CalHelixFinderAlg::fillFaceOrderedHits(CalHelixFinderData& Helix) {

//-----------------------------------------------------------------------------
// set the CaloCluster of CalHelixFinderAlg: this info is stored in the TimeCluster
//-----------------------------------------------------------------------------
    setCaloCluster(Helix);

    static const float pi(M_PI);
    static const float twopi(2*pi);

    float clPhi(-9999.);

    if (fCaloTime > 0) clPhi = polyAtan2(fCaloY,fCaloX);

    const vector<StrawHitIndex>& shIndices = Helix._timeCluster->hits();
    ChannelID cx, co;

    int     size           = Helix._timeCluster->nhits();
    int     nFiltPoints(0), nFiltStrawHits(0);
    //--------------------------------------------------------------------------------
    // if (Helix.shpos() != 0) {
    int loc;
    StrawHitFlag flag;

    //sort the hits by z coordinate
    ComboHitCollection ordChCol;
    ordChCol.reserve(size);

    if (_debug >0 ){
      printf("-----------------------------------------------------------------------------------/n");
      printf("[CalHelixFinderAlg::fillFaceOrderedHits]  Plane   Panel  Layer Straw     x          y           z  \n");
    }

    for (int i=0; i<size; ++i) {
      loc = shIndices[i];
      const ComboHit& ch  = Helix.chcol()->at(loc);
      flag = ch.flag();
      //-----------------------------------------------------------------------------
      // select hits: don't reuse straw hits
      //-----------------------------------------------------------------------------
      int good_hit = flag.hasAllProperties(_hsel  );
      int bkg_hit  = flag.hasAnyProperty  (_bkgsel);
      // int used_hit = flag.hasAnyProperty  (StrawHitFlag::calosel);
      // if (good_hit && (! bkg_hit) && (! used_hit)) {
      if (good_hit && (! bkg_hit) ) {

        if (ch.energyDep() > _maxHitEnergy)                 continue;

        //skip the hit if it doesn't rely on the semi-plane where the calo-lcuster is
        if (_filter) {
          float chPhi = polyAtan2(ch.pos().y(), ch.pos().x());
          float dphi  = chPhi - clPhi;

          if (dphi >  pi) dphi -= twopi;
          if (dphi < -pi) dphi += twopi;

          if (fabs(dphi) > pi/2)                            continue;
        }

        if (_debug >0 ){
          const mu2e::Straw* straw;

          straw = &_tracker->getStraw(ch.strawId());//ndex());

          printf("[CalHelixFinderAlg::fillFaceOrderedHits]  %5i  %5i   %5i   %5i   %8.3f   %8.3f    %10.3f\n",
                 straw->id().getPlane(),
                 straw->id().getPanel(),
                 straw->id().getLayer(),
                 straw->id().getStraw(),
                 ch.pos().x(), ch.pos().y(), ch.pos().z());
        }

        ordChCol.push_back(ComboHit(ch));
      }
    }
    std::sort(ordChCol.begin(), ordChCol.end(),panelcomp());

    for (unsigned i=0; i<ordChCol.size(); ++i) {
      ComboHit& ch = ordChCol[i];

      cx.Station                 = ch.strawId().station();//straw.id().getStation();
      cx.Plane                   = ch.strawId().plane() % 2;//straw.id().getPlane() % 2;
      cx.Face                    = ch.strawId().face();
      cx.Panel                   = ch.strawId().panel();//straw.id().getPanel();

      // get Z-ordered location
      Helix.orderID(&cx, &co);

      int os       = co.Station;
      int of       = co.Face;
      //FIXME!!! COMPARE OF WITH ch.sid().uniqueFace()
      int op       = co.Panel;

      int       stationId = os;
      int       faceId    = of + stationId*StrawId::_nfaces*FaceZ_t::kNPlanesPerStation;//FaceZ_t::kNFaces;
      //        int       panelId   = op + faceId*FaceZ_t::kNPanels;//PerFace;
      FaceZ_t*  fz        = &Helix._oTracker[faceId];
      PanelZ_t* pz        = &fz->panelZs[op];

      if ((of < 0) || (of >  StrawId::_nfaces*FaceZ_t::kNPlanesPerStation  )) printf(" >>> ERROR: wrong face    number: %i\n",of);
      if ((op < 0) || (op >= FaceZ_t::kNPanels )) printf(" >>> ERROR: wrong panel   number: %i\n",op);

      Helix._chHitsToProcess.push_back(mu2e::ComboHit(ch));

      if (pz->idChBegin < 0 ){
        pz->idChBegin = Helix._chHitsToProcess.size() - 1;
        pz->idChEnd   = Helix._chHitsToProcess.size();
      } else {
        pz->idChEnd   = Helix._chHitsToProcess.size();
      }

      if (fz->idChBegin < 0 ){
        fz->idChBegin = Helix._chHitsToProcess.size() - 1;
        fz->idChEnd   = Helix._chHitsToProcess.size();
      } else {
        fz->idChEnd   = Helix._chHitsToProcess.size();
      }

      ++nFiltPoints;
      nFiltStrawHits += ch.nStrawHits();
    }

    Helix._nFiltPoints    = nFiltPoints;
    Helix._nFiltStrawHits = nFiltStrawHits;

    if (_debug > 0) printXYZP(Helix);
  }




//----------------------------------------------------------------------------
//2015-01-17 G. Pezzullo: the following procedure looks the hit with
// z-coordinate smaller then the seeding one and calculates distance from
// prediction in order to check if they are good or outliers
//----------------------------------------------------------------------------
  void CalHelixFinderAlg::rescueHitsBeforeSeed(CalHelixFinderData& Helix){
    FaceZ_t*    facez  = &Helix._oTracker[Helix._seedIndex.face];
    PanelZ_t*   panelz = &facez->panelZs[Helix._seedIndex.panel];
    mu2e::ComboHit* hit = &Helix._chHitsToProcess[Helix._seedIndex.panelHitIndex];

    float      weight(-1), radius, phi0, dfdz, x0, y0;
    dfdz        = Helix._dfdz;
    phi0        = Helix._fz0 + dfdz*Helix._zFace[Helix._seedIndex.face];
    x0          = Helix._center.x();
    y0          = Helix._center.y();
    radius      = Helix._radius;

    float      dx,dy,phi,max_dist;
    XYZVectorF      shPos, hePos, strawDir, helCenter(x0, y0, 0);

    float      deltaZ(0.);
    float      distXY(0.);
    float      dist(0.), dist2(0.); // help parameter for storing strawhit position residual
    int         rescuedStrawHits(0), rescuedPoints(0);

    char banner[]="CalHelixFinderAlg::rescueHitsBeforeSeed";

    if (_debug > 0) {
      printf("[%s:BEGIN] x0 = %5.3f y0 = %5.3f radius = %5.3f phi0 = %5.5f dfdz = %5.6f chi2 = %5.3f \n", banner,
             x0, y0, radius, phi0, dfdz , Helix._sxy.chi2DofCircle());
      printf("[%s] SeedIndex = %i N-points = %5.3f\n",  banner, Helix._seedIndex.panel, Helix._sxy.qn()-1);//FIXME!
      if (Helix._seedIndex.panel >= 0) {//FIXME!
        printf("[%s] index      Z        xi      yi       xp       yp       X0        Y0         R        dfdZ  dXY(pred) dXY(seed) dZ(seed)     wt\n",banner);
        printf("[%s]-------------------------------------------------------------------------------------------------------------------------------\n",banner);
      }
    }
//-----------------------------------------------------------------------------
// given a helix candidate, move upstream and pick up points with Z < _xyzp[fSeedIndex].z
//-----------------------------------------------------------------------------
    float           lastFacez = Helix._zFace[Helix._seedIndex.face];

    int   firstFace = Helix._seedIndex.face;
    if (firstFace<0) firstFace = 0;
    for (int f=firstFace; f<StrawId::_ntotalfaces;  ++f){
      facez     = &Helix._oTracker[f];
      float          faceHitChi2(1e10);
      HitInfo_t      goodFaceHit;

      int       panelIndex(FaceZ_t::kNPanels-1);
      //-----------------------------------------------------------------------------
      // dfdz = tanLambda/radius; phi0 is the last found hit phi
      //-----------------------------------------------------------------------------
      deltaZ    = Helix._zFace[f] - lastFacez;
      phi       = phi0 + (deltaZ)*dfdz;
      //evaluate the helix prediction using the z coordinate of the panel
      hePos     = XYZVectorF(x0 + radius*std::cos(phi),
                         y0 + radius*std::sin(phi),
                         Helix._zFace[f]);
      //check the Panel-phi wrt to the DS center
      float  hePosPhi = polyAtan2(hePos.y(), hePos.x());
      if (hePosPhi < 0) hePosPhi = hePosPhi + 2*M_PI;

      for (int p=panelIndex; p>=0;  --p){
        panelz = &facez->panelZs[p];
        int  nhitsPerPanel  = panelz->nChHits();

        if (nhitsPerPanel == 0)                                              continue;
        if ( (f==Helix._seedIndex.face) && (p==Helix._seedIndex.panel) && (Helix._seedIndex.panelHitIndex >=0) ) nhitsPerPanel = Helix._seedIndex.panelHitIndex - panelz->idChBegin;//the seedHit is already clusterized!


        float  deltaPhi = hePosPhi - Helix._phiPanel[f*FaceZ_t::kNPanels + p];
        if ( deltaPhi > M_PI ) deltaPhi -= 2*M_PI;
        if ( deltaPhi < -M_PI) deltaPhi += 2*M_PI;
        if ( fabs(deltaPhi) > _maxPanelToHelixDPhi)                             continue;

        for (int  i=nhitsPerPanel-1;i>=0; --i){
          int   index = panelz->idChBegin + i;
          if (Helix._hitsUsed[index] >= 1)                    continue;
          hit       = &Helix._chHitsToProcess[index];
          shPos     = hit->_pos;
          strawDir  = hit->vDir();

          dx        = hePos.x() - shPos.x();
          dy        = hePos.y() - shPos.y();
          dist2     = dx*dx + dy*dy;
          dist      = std::sqrt(dist2);

          if (_debug > 10) {
            printf("[%s:LOOP] %5i %9.3f %8.3f %8.3f %8.3f %8.3f %9.3f %9.3f %9.3f %8.5f %8.3f %8.3f %8.3f %8.3f",
                   banner,i,shPos.z(),shPos.x(),shPos.y(),hePos.x(),hePos.y(),
                   x0,y0,radius,dfdz,dist,distXY,deltaZ,weight);
          }

          max_dist = _distPatRec + _dfdzErr*std::fabs(deltaZ);
          if (dist <= max_dist) {
            if (dist < faceHitChi2){
              goodFaceHit.face          = f;
              goodFaceHit.panel         = p;
              goodFaceHit.panelHitIndex = index;
              faceHitChi2               = dist;
            }
          }else {
            if (_debug > 10) {
              printf(" missed\n");
            }
          }
        }//end loop pver the hits within the panel
      }//end loop over the panels within a face

      if (goodFaceHit.panel < 0)                    continue;
      //get the best hit found
      int  index = goodFaceHit.panelHitIndex;
      hit        = &Helix._chHitsToProcess[index];

      lastFacez  = Helix._zFace[f];

      // add point to the helixfithack result objet
      weight     = calculateWeight(*hit, helCenter, radius);
      Helix._sxy.addPoint(hit->pos().x(),hit->pos().y(),weight);
      Helix._nXYSh += hit->nStrawHits();

      // update helix parameters
      x0      = Helix._sxy.x0();
      y0      = Helix._sxy.y0();
      radius  = Helix._sxy.radius();

      helCenter.SetX(x0);
      helCenter.SetY(y0);

      // float dx  = (hit->pos().x() - helCenter.x());
      // float dy  = (hit->pos().y() - helCenter.y());
      phi0       =  hit->_hphi;//polyAtan2(dy,dx);
      //      hit->_hphi = phi0;

      // 2019-02-05: gianipez comment; for future development, resolve the 2pi ambig once you find a new hit.
      // this is not necessary if the 2pi ambig is resolved for all the hits in the function
      // ::doLinearFitPhiZ(...)
      //      resolve2PiAmbiguity(hit, dfdz, Helix._dfdz, Helix._phi0);


      //update hit info
      hit->_xyWeight   = weight;
      Helix._hitsUsed[index] = 1;

      //      float dzFromSeed = facez->z - seedFacez->z;         // expected to be negative
      // hit->_dzFromSeed  = dzFromSeed;
      // hit->_drFromPred  = faceHitChi2;//[t];

      rescuedStrawHits += hit->nStrawHits();
      ++rescuedPoints;

      if (_debug > 0) {
        printf("rescued %08x %2i %12.5f %12.5f %12.5f\n",
               *((int*) &hit->_flag), Helix._hitsUsed[index],
               hit->pos().x(), hit->pos().y(), hit->pos().z());
      }

    }//end loop over the faces

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
      printf(" SeedIndex: %i N(rescued points): %i\n",Helix._seedIndex.panel,rescuedPoints);//FIXME!
    }

    if (_diag) Helix._diag.n_rescued_points_16 = rescuedPoints;
  }

//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::printInfo(CalHelixFinderData&  Helix){
    const char banner [] = "CalHelixFinderAlg::printInfo";
    float dr(0), dx, dy;

    if (_debug > 0) {
      printf("[%s] N(points): %3i x0: %12.5f y0: %12.5f r: %12.5f chi2c: %12.5f phi0: %5.5f dfdz: %5.6f chi2l: %5.3f\n",
             banner,
             Helix._nStrawHits,
             Helix._sxy.x0(),Helix._sxy.y0(),Helix._sxy.radius(),Helix._sxy.chi2DofCircle(),
             Helix._fz0, Helix._dfdz , Helix._szphi.chi2DofLine());

      FaceZ_t*        facez(0);
      PanelZ_t*       panelz(0);
      mu2e::ComboHit* hit(0);

      for (int f=StrawId::_ntotalfaces-1; f>=0; --f){
        facez     = &Helix._oTracker[f];
        for (int p=0; p<FaceZ_t::kNPanels; ++p){
          panelz =  &facez->panelZs[p];//&Helix._oTracker[p];
          int  nhits          = panelz->nChHits();
          for (int i=0; i<nhits; ++i){
            hit = &Helix._chHitsToProcess[panelz->idChBegin + i];
            dx = hit->pos().x() - Helix._sxy.x0();
            dy = hit->pos().y() - Helix._sxy.y0();
            dr = sqrt(dx*dx+dy*dy) - Helix._sxy.radius();
            printf("[%s] %08x %6i %3i %6i %12.5f %12.5f %12.5f %10.3f\n",banner,
                   *((int*) &hit->_flag),  int(hit->index()), i, hit->strawId().straw(),
                   hit->pos().x(), hit->pos().y(), hit->pos().z(), dr
                   );//FIXME!
          }
        }//end loop over the panels within the face
      }//end loop opver the faces
    }
  }

//-----------------------------------------------------------------------------
// this routine simply checks '_indicesTrkCandidate' array and for negative
// indices sets the 'outlier' flag to the corresponding 'xyzp'
// no actual check of residuals is performed
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::filterUsingPatternRecognition(CalHelixFinderData& Helix) {

    if (Helix._seedIndex.panel < 0) return;

    int            nActive(0), nActive_hel(0);
    int            nSh = Helix._nFiltStrawHits;

    float         straw_mean_radius(0), chi2_global_helix(0), total_weight(0);
    float         x_center(Helix._center.x()), y_center(Helix._center.y()), radius(Helix._radius);
    float         fz0(Helix._fz0), lambda(1./Helix._dfdz);
    XYZVectorF         hel_pred(0., 0., 0.);

    PanelZ_t*      panelz(0);
    FaceZ_t*       facez(0);

    mu2e::ComboHit* hit(0); bool isFirst(true);

    for (int f=StrawId::_ntotalfaces-1; f>=0; --f){
      facez     = &Helix._oTracker[f];
      for (int p=0; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];// &Helix._oTracker[p];
        int  nhits     = panelz->nChHits();
        if (nhits == 0)             continue;

        for (int i=0; i<nhits; ++i){
          int   index = panelz->idChBegin + i;
          hit = &Helix._chHitsToProcess[index];
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

          //          int index = facez->evalUniqueHitIndex(f,p,i);
          if (Helix._hitsUsed[index] != 1)        hit->_flag.merge(StrawHitFlag::outlier); // *FIXME* ? ? do we want to call this an outlier ?
          else  {
            // ++nActive;
            nActive          += hit->nStrawHits();
            float  x         = hit->pos().x();
            float  y         = hit->pos().y();
            float  z         = Helix._zFace[f];
            float  phi_pred  = fz0 + z/lambda;
            float  x_pred    = x_center + radius*cos(phi_pred);
            float  y_pred    = y_center + radius*sin(phi_pred);
            hel_pred.SetX(x_pred);
            hel_pred.SetY(y_pred);
            float  weight    = calculateWeight(*hit, hel_pred, radius)*_weight3D;
            float  x_resid2  = (x - x_pred)*(x - x_pred);
            float  y_resid2  = (y - y_pred)*(y - y_pred);
            float  hitResi2  = (x_resid2 + y_resid2)*weight;
            if (hitResi2 > _chi2hel3DMax) {
              hit->_flag.merge(StrawHitFlag::outlier);  // *FIXME* ? ? do we want to call this an outlier ?
            } else {
              // ++nActive_hel;
              nActive_hel      += hit->nStrawHits();
              chi2_global_helix = chi2_global_helix + hitResi2;
              straw_mean_radius = straw_mean_radius + sqrt(x*x + y*y)*weight;
              total_weight      = total_weight + weight;
            }
          }

          if (_debug > 10) {
            XYZVectorF*     shPos = &hit->_pos;
            int         is_outlier    = hit->_flag.hasAllProperties(StrawHitFlag::outlier);
            string      type;
            if      ((f == Helix._seedIndex.face) && (p == Helix._seedIndex.panel) && (i == Helix._seedIndex.panelHitIndex)) type = "seed";
            else if ((f == Helix._seedIndex.face) && (p == Helix._candIndex.panel) && (i == Helix._candIndex.panelHitIndex)) type = "cand";
            float dist(0);// = hit->_drFromPred;
            float dz  (0);// = hit->_dzFromSeed;//FIXME!
            printf("[CalHelixFinderAlg::filterUsingPatternRecognition] %5i %5i %4i %4s  %8.3f %8.3f %9.3f %8.3f %8.3f\n",
                   i,hit->strawId().straw(),is_outlier,type.data(),shPos->x(),shPos->y(),shPos->z(),dist,dz);
          }
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
      Helix._diag.chi2d_helix       = (nActive > 0) ? chi2_global_helix/float(nActive_hel - 5.): -1;
    }
  }

//-----------------------------------------------------------------------------
  float CalHelixFinderAlg::deltaPhi(float phi1, float phi2){
    static const float pi(M_PI);
    static const float twopi(2*pi);
    float dphi = fmod(phi2-phi1,twopi);
    if(dphi>pi)dphi -= twopi;
    if(dphi<-pi)dphi += twopi;
    return dphi;
  }


//--------------------------------------------------------------------------------
  void CalHelixFinderAlg::searchBestTriplet   (CalHelixFinderData& Helix, CalHelixFinderData& TmpHelix, int UseMPVdfdz){
    int       nSh = Helix._nFiltStrawHits;
    int       nHitsTested(0);

    FaceZ_t*  facez(0);
    PanelZ_t* panelz(0);

    if (_debug > 5) { printf("[CalHelixFinderAlg::searchBestTriplet] BEGIN, nSh = %d\n",nSh); }
    for (int f=StrawId::_ntotalfaces-1; f>=0; --f){
      ///if (Helix._zFace[f] < _maxZTripletSearch)     break;
      facez     = &Helix._oTracker[f];
      for (int p=0; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];
        int       nhits  = panelz->nChHits();
        for (int i=0; i<nhits; ++i){
          if (Helix._nStrawHits > (nSh - nHitsTested))   continue;
          if ((nSh - nHitsTested) < _minNHits        )   continue;
          //clear the info of the tmp object used to test the triplet
          TmpHelix.clearResults();

          HitInfo_t          seed(f,p,panelz->idChBegin + i);
          findTrack(seed,TmpHelix,UseMPVdfdz);

          nHitsTested += Helix._chHitsToProcess[panelz->idChBegin + i].nStrawHits();

          //compare tripletHelix with bestTripletHelix
          //2019-02-08: gianipez chanceg the logic;
          //2019-02-15: gianipez put the old logic back. FIXME!
          if (( TmpHelix._nStrawHits >  Helix._nStrawHits) ||
              ((TmpHelix._nStrawHits == Helix._nStrawHits) && (TmpHelix._helixChi2 < Helix._helixChi2))) {
            // int   deltaNSh = TmpHelix._nStrawHits -  Helix._nStrawHits;
            // if ( ( deltaNSh >=  _minDeltaNShPatRec)  ||
            //      ( deltaNSh>=0 && (deltaNSh-_minDeltaNShPatRec < 0) && (TmpHelix._helixChi2 < Helix._helixChi2)) ||
            //      ((TmpHelix._nStrawHits == Helix._nStrawHits) && (TmpHelix._helixChi2 < Helix._helixChi2)) ) {
            Helix = TmpHelix;
          }
          if (_debug > 5) {
            printf("[CalHelixFinderAlg::doPatternRecognition]: calling findTrack(i=%i,Helix,useDefaltDfDz=FALSE,useMPVdfdz=%i)",panelz->idChBegin +i,UseMPVdfdz);
            printf(" : np=%3i _goodPointsTrkCandidate=%3i\n",nSh,Helix._nStrawHits);
          }
        }//end loop over the hits on the panel
      }//end panels loop
    }//end faces loop
    //}//end loop over the fake caloclusters
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

    _findTrackLoopIndex = 1;                 // debugging
    if (_debug != 0) {
      printf("[CalHelixFinderAlg::doPatternRecognition]: before first SearchBestTriplet \n");
      printInfo(Helix);
    }
    searchBestTriplet(Helix, tripletHelix);
    if (_debug != 0) {
      printf("[CalHelixFinderAlg::doPatternRecognition]: after first SearchBestTriplet \n");
      printInfo(Helix);
    }
    //-----------------------------------------------------------------------------
    // 2014-11-09 gianipez: if no track was found requiring the recalculation of dfdz
    // look for a track candidate using the default value of dfdz and the target center
    //-----------------------------------------------------------------------------
    _findTrackLoopIndex = 2;                 // *DEBUGGING*
    if (fUseDefaultDfDz == 0) {
      searchBestTriplet(Helix, tripletHelix, useMPVdfdz);
      if (_debug != 0) {
        printf("[CalHelixFinderAlg::doPatternRecognition]: after second SearchBestTriplet \n");
        printInfo(Helix);
      }
   }

    if (_debug == 0){
      _debug  = _debug2;
      _debug2 = 0;
    }

    char banner[200];
    bool rc;
    int  rs, usePhiResid;

    float circleHits(0.0), phiHits(0.0), nHitsRatio(0.0);

    if ((Helix._seedIndex.panel < 0) || (Helix._nXYSh < _minNHits) ) goto  PATTERN_RECOGNITION_END;

    // 2015-01-17 G. Pezzullo: rescue points with z-coordinate less than the seed hit
    if (_debug != 0) {
      printf("[CalHelixFinderAlg::doPatternRecognition]: calling rescueHitsBeforeSeed\n");
      printInfo(Helix);
    }

    circleHits = Helix._sxy.qn();

    rc = doLinearFitPhiZ(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), useIntelligentWeight);

    phiHits = Helix._szphi.qn();

    nHitsRatio = circleHits/(phiHits + 1e-6);

    if (_diag > 0) {
      Helix._diag.nHitsRatio = nHitsRatio;
    }

    if (nHitsRatio > _maxNHitsRatio) goto PATTERN_RECOGNITION_END;
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

    refineHelixParameters(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), banner, _debug);

    if (_debug != 0)  printInfo(Helix);
//---------------------------------------------------------------------------------------
// use the results of the helix search to see if points along the track can be rescued
//---------------------------------------------------------------------------------------
    if ((Helix._nZPhiSh < _minNHits) || (!rc)) usePhiResid = 0;
    else                                       usePhiResid = 1;

    rescueHits(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), usePhiResid);

    if ((Helix._nXYSh - 1) != Helix._nZPhiSh) rc = doLinearFitPhiZ(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), useIntelligentWeight);//the factor "-1" takes into account that the XY fit includes the target center

    if (_debug != 0)  printInfo(Helix);
//--------------------------------------------------------------------------------------------------------------
// 2015-03-25 G. Pezzu added the following call to findDfDz(...) in order to help the fitter on finding
// the more reliable value of dfdz which is needed for resolving the 2pi ambiguity.
// Since in the previous step we could have rescued few points, that would give us an help!
// re-evaluate the df/dz and phi0 including rescued hits and new XY parameters
//--------------------------------------------------------------------------------------------------------------
    if (Helix._nZPhiSh < _minNHits || (!rc)){
      rs = findDfDz(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1));

      if (rs == 1) {                        // update Helix Z-phi part
        Helix._dfdz = _hdfdz;
        Helix._fz0  = _hphi0;
      }
    }

    rc = doLinearFitPhiZ(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), useIntelligentWeight);

    if (rc) {
      usePhiResid = 1;
      rescueHits(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), usePhiResid);
      if ((Helix._nXYSh - 1) != Helix._nZPhiSh) rc = doLinearFitPhiZ(Helix, HitInfo_t(StrawId::_ntotalfaces-1,0,-1), useIntelligentWeight);  //the factor "-1" takes into account that the XY fit includes the target center

      if (_debug != 0)  printInfo(Helix);
      strcpy(banner,"refineHelixParameters-after-doLinearFitPhiZ");
      refineHelixParameters(Helix,HitInfo_t(StrawId::_ntotalfaces-1,0,-1),banner,_debug);
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
  float  CalHelixFinderAlg::calculateWeight(const mu2e::ComboHit& Hit      ,
                                             const XYZVectorF&         HelCenter,
                                             float                Radius   ) {

    float    transErr = 5./sqrt(12.);
    //scale the error based on the number of the strawHits that are within the mu2e::ComboHit
    if (Hit.nStrawHits() > 1) transErr *= 1.5;
    float    transVar = transErr*transErr;

    float x   = Hit.pos().x();
    float y   = Hit.pos().y();
    float dx  = x-HelCenter.x();
    float dy  = y-HelCenter.y();
    float dxn = dx*Hit.vDir().x()+dy*Hit.vDir().y();

    float costh2 = dxn*dxn/(dx*dx+dy*dy);
    float sinth2 = 1-costh2;

    float e2     = Hit.wireVar()*sinth2+transVar*costh2;
    float wt     = 1./e2;
                                                    // scale the weight for having chi2/ndof distribution peaking at 1
    wt *= _weightXY;

    return wt;
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  float  CalHelixFinderAlg::calculatePhiWeight(const mu2e::ComboHit& Hit      ,
                                                const XYZVectorF&         HelCenter,
                                                float                Radius   ,
                                                int                   Print    ,
                                                const char*           Banner   ) {
    //    float    transErr = 5./sqrt(12.);
    //scale the error based on the number of the strawHits that are within the mu2e::ComboHit
    //    if (Hit.nStrawHits() > 1) transErr *= 1.5;
    //    float    transVar = transErr*transErr;

    float x  = Hit.pos().x();
    float y  = Hit.pos().y();
    float dx = x-HelCenter.x();
    float dy = y-HelCenter.y();

//-----------------------------------------------------------------------------
// if dr(dx,dy) is orthogonal to the wire, costh = 1
//-----------------------------------------------------------------------------
    float dxn    = dx*Hit.vDir().x()+dy*Hit.vDir().y();
    float costh2 = dxn*dxn/(dx*dx+dy*dy);
    float sinth2 = 1-costh2;

    //    float e2     = Hit.wireVar()*costh2+transVar*sinth2;
    float e2     = Hit.wireVar()*costh2+Hit.transVar()*sinth2;
    float wt     = Radius*Radius/e2;
    wt           *= _weightZPhi;

    if (Print > 0) {
      float dr = calculateRadialDist(Hit.pos(),HelCenter,Radius);
      printf("[CalHelixFinderAlg::%s] %9.3f %9.3f %10.5f %10.5f %10.5f %10.5f %12.5e %10.3f\n",
                               Banner, x, y, dx, dy, costh2, sinth2, e2, dr);
    }

    return wt;
  }

//--------------------------------------------------------------------------------
// calculate the radial distance of a straw hit from the helix prediction
//--------------------------------------------------------------------------------
  float  CalHelixFinderAlg::calculateRadialDist (const XYZVectorF& HitPos   ,
                                                  const XYZVectorF& HelCenter,
                                                  float        Radius   ) {
    float dx = HitPos.x()-HelCenter.x();
    float dy = HitPos.y()-HelCenter.y();
    float dr = sqrt(dx*dx+dy*dy)-Radius;

    return dr;
  }


//-----------------------------------------------------------------------------
  void   CalHelixFinderAlg::doWeightedCircleFit (CalHelixFinderData& Helix,
                                                 HitInfo_t           SeedIndex,
                                                 XYZVectorF&         HelCenter,
                                                 float&              Radius   ,
                                                 int                 Print    ,
                                                 const char*         Banner   ) {
    float     wt;
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
    if (_targetcon == 1){
      Helix._sxy.addPoint(0.,0.,1./900.);
      Helix._nXYSh += 1;
      Helix._nComboHits += 1;
    }

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doWeightedCircleFit] BEGIN: x0 = %8.3f y0 = %8.3f radius = %8.3f chi2dof = %8.3f\n",
             HelCenter.x(),HelCenter.y(),Radius,Helix._sxy.chi2DofCircle());
      if (_debug > 10) {
        printf("[CalHelixFinderAlg::doWeightedCircleFit:LOOP] Index      X          Y         Z          wt        wireNx     wireNy\n");
      }
    }

    PanelZ_t*      panelz(0);
    FaceZ_t*       facez(0);

    mu2e::ComboHit* hit   (0);

    for (int f=SeedIndex.face; f>=0; --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;

      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];//&Helix._oTracker[p];
        int  nhits          = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhits == 0)                                                                    continue;
        if ( (f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhits; ++i){
          int     index = panelz->idChBegin + i;
          hit = &Helix._chHitsToProcess[index];
          //          int index = facez->evalUniqueHitIndex(f,p,i);
          if (Helix._hitsUsed[index] != 1)                     continue;

          wt             = calculateWeight(*hit,HelCenter,Radius);
          hit->_xyWeight = wt;

          Helix._sxy.addPoint(hit->_pos.x(),hit->_pos.y(),wt);
          Helix._nXYSh      += hit->nStrawHits();
          Helix._nComboHits += 1;

          if (_debug > 10) {
            printf("[CalHelixFinderAlg::doWeightedCircleFit:LOOP] %4i %10.3f %10.3f %10.3f %10.3e %10.4f %10.4f\n",
                   (int)hit->index(), hit->_pos.x(), hit->_pos.y(), hit->_pos.z(), wt, hit->vDir().x(), hit->vDir().y());
          }
        }
      }//end panels loop
    }
                                        // update helix info
    Radius  = Helix._sxy.radius();
    HelCenter.SetX(Helix._sxy.x0());
    HelCenter.SetY(Helix._sxy.y0());

    //2018-10-04: gianipez added the follwoing line to store the new valuse of circle center and radius in the Helix
    Helix._center.SetX(Helix._sxy.x0());
    Helix._center.SetY(Helix._sxy.y0());
    Helix._radius  = Helix._sxy.radius();


    if (_debug > 5) {
      printf("[CalHelixFinderAlg::doWeightedCircleFit:END] : npt = %3.0f  chi2dof = %8.3f x0 = %8.3f y0 = %8.3f radius = %8.3f\n",
             Helix._sxy.qn(),Helix._sxy.chi2DofCircle(),HelCenter.x(),HelCenter.y(),Radius);
    }
  }


//-----------------------------------------------------------------------------
// this is a rather "primitive" definition of the worst hit, should do for now
//-----------------------------------------------------------------------------
  void    CalHelixFinderAlg::searchWorstHitWeightedCircleFit(CalHelixFinderData& Helix,
                                                             HitInfo_t           SeedIndex,
                                                             const XYZVectorF&       HelCenter,
                                                             float&             Radius,
                                                             HitInfo_t&          Iworst,
                                                             float&             HitChi2Worst)
  {
    HitChi2Worst         = _hitChi2Max;
    Iworst.face          = -1;
    Iworst.panel         = -1;
    Iworst.panelHitIndex = -1;

    float     dr, hitChi2;

    mu2e::ComboHit* hit(0);
    FaceZ_t*        facez(0);
    PanelZ_t*       panelz(0);
    //    FaceZ_t*       seedFacez = &Helix._oTracker[Helix._seedIndex.face];


    for (int f=SeedIndex.face; f>=0;  --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];//Helix._oTracker[p];
        int  nhitsPerPanel  = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhitsPerPanel == 0)                                                            continue;
        if ((f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
          int     index = panelz->idChBegin + i;
          hit = &Helix._chHitsToProcess[index];
          if (Helix._hitsUsed[index] != 1)                    continue;

          dr      = calculateRadialDist(hit->_pos,HelCenter,Radius);
          hitChi2 = dr*dr*hit->_xyWeight;

          // store info out the radial residual
          // if ((SeedIndex.face == 0) && (SeedIndex.panel == 0) && (SeedIndex.panelHitIndex == 0)) {
          // hit->_drFromPred  = fabs(dr);//hitChi2;
          // float dzFromSeed = facez->z - seedFacez->z;              // expected to be positive (non-negative)
          //            hit->_dzFromSeed  = dzFromSeed;//FIXME!
          // }

          if (hitChi2 > HitChi2Worst) {
            HitChi2Worst         = hitChi2;
            Iworst.face          = f;
            Iworst.panel         = p;
            Iworst.panelHitIndex = index;
          }
        }
      }//end panels loop
    }//end faces loop
  }

//--------------------------------------------------------------------------------
// IWorst is always defined
// returns the index of the hit which provides the highest contribute to the chi2
//--------------------------------------------------------------------------------
  void    CalHelixFinderAlg::cleanUpWeightedCircleFit(CalHelixFinderData& Helix,
                                                      HitInfo_t          SeedIndex,
                                                      HitInfo_t&         IWorst)
  {
    LsqSums4   sxy;
    float     chi2, chi2_min (-1.), x, y;

    //reset the coordinates of the worst hit found previousl
    IWorst.face          = -1;
    IWorst.panel         = -1;
    IWorst.panelHitIndex = -1;

    mu2e::ComboHit* hit(0);
    PanelZ_t*       panelz(0);
    FaceZ_t*        facez(0);

    for (int f=SeedIndex.face; f>=0; --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;

      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];//&Helix._oTracker[p];
        int  nhitsPerPanel  = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhitsPerPanel == 0)                                                             continue;
        if ( (f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
          int     index = panelz->idChBegin + i;
          hit = &Helix._chHitsToProcess[index];
          //          int index = facez->evalUniqueHitIndex(f,p,i);
          if (Helix._hitsUsed[index] != 1)                    continue;

          sxy = Helix._sxy;

          x  = hit->_pos.x();
          y  = hit->_pos.y();

          sxy.removePoint(x, y, hit->_xyWeight);

          chi2  = sxy.chi2DofCircle();

          if ((chi2 < chi2_min) || ( (i == SeedIndex.panelHitIndex) && (p == SeedIndex.panel) && (f == SeedIndex.face)) ) {
            chi2_min             = chi2;
            IWorst.face          = f;
            IWorst.panel         = p;
            IWorst.panelHitIndex = index;
          }
        }
      }//end panels loop
    }//end faces loop
  }

//-----------------------------------------------------------------------------
// use hits only, at this point the cluster is no longer needed
//-----------------------------------------------------------------------------
  int CalHelixFinderAlg::refineHelixParameters(CalHelixFinderData& Trk,
                                               HitInfo_t           SeedIndex,
                                               const char*         Banner,
                                               int                 Print  ) {
    auto           hitsUsed = Trk._hitsUsed;
    float         x, y, r, r_start;
    float         hitChi2Worst;

    int            pointsRemoved(0);

    HitInfo_t     iworst;//(-1, -1);
    float         wtWorst;
    float         chi2, chi2_min;

    XYZVectorF          hitPos, strawDir, helCenter, helCenter_start;
    mu2e::ComboHit* hit(0);
    //    FaceZ_t*        facez(0);
    //    PanelZ_t*       panelz(0);

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
    // doWeightedCircleFit (Trk,SeedIndex,helCenter,r,Print,Banner);

    //now initialize the LsqSum4 variable
    // sxyw = Trk._sxy;

    searchWorstHitWeightedCircleFit(Trk,SeedIndex,helCenter,r,iworst,hitChi2Worst);

    chi2     = Trk._sxy.chi2DofCircle();
    chi2_min = chi2;

    if (_debug > 5) {
      printf("[CalHelixFinderAlg::refineHelixParameters] npt = %3.0f x0 = %8.3f y0 = %8.3f radius = %8.3f chi2 = %8.3f iworst=%3i chi2Worst = %8.3f\n",
             Trk._sxy.qn(),Trk._sxy.x0(),Trk._sxy.y0(),Trk._sxy.radius(),Trk._sxy.chi2DofCircle(),iworst.panel,hitChi2Worst);
    }

    if ((chi2 <= _chi2xyMax) && (hitChi2Worst <= _hitChi2Max)) goto F_END;
//-----------------------------------------------------------------------------
// one of the chi2's is above the threshold, cleanup is needed
//-----------------------------------------------------------------------------
    if (_debug > 5) printf("[CalHelixFinderAlg::refineHelixParameters] : START CLEANUP\n");
  NEXT_ITERATION:;

    cleanUpWeightedCircleFit(Trk,SeedIndex,iworst);

    if (iworst.panel >= 0) {
      // facez   = &Trk._oTracker[iworst.face];
      // panelz  = &facez->panelZs[iworst.panel];
      hitUsedIndex = iworst.panelHitIndex;
      hit     = &Trk._chHitsToProcess[hitUsedIndex];
      x       = hit->_pos.x();
      y       = hit->_pos.y();
      wtWorst = hit->_xyWeight;

                                        // remove point from the track, this is why need to return weights
      Trk._sxy.removePoint(x, y, wtWorst);
      Trk._nXYSh -=hit->nStrawHits();

      //      hitUsedIndex = facez->evalUniqueHitIndex(iworst);
      Trk._hitsUsed[hitUsedIndex] = 0;

      Trk._nComboHits -= 1;
      Trk._nStrawHits -= hit->nStrawHits();
      ++pointsRemoved;

      if (_debug > 5) {
        printf("[CalHelixFinderAlg::refineHelixParameters]  x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %5.5f chi2Maxxy = %5.5f index point removed = %i\n",
               Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(), _chi2xyMax, iworst.panel);//FIXME!
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
    if (iworst.panel >= 0) {
      // facez   = &Trk._oTracker[iworst.face];
      // panelz  = &facez->panelZs[iworst.panel];//panelz  = &Trk._oTracker[iworst.Panel];
      hitUsedIndex = iworst.panelHitIndex;
      hit     = &Trk._chHitsToProcess[hitUsedIndex];
      x       = hit->_pos.x();
      y       = hit->_pos.y();
      wtWorst = hit->_xyWeight;
                                        // remove point from the track and mark it
      Trk._sxy.removePoint(x, y, wtWorst);
      Trk._nXYSh -=hit->nStrawHits();

      //      hitUsedIndex = facez->evalUniqueHitIndex(iworst);
      Trk._hitsUsed[hitUsedIndex] = 0;

      Trk._nComboHits -= 1;
      Trk._nStrawHits -= hit->nStrawHits();
      ++pointsRemoved;

      if (_debug > 5) {
        printf("[CalHelixFinderAlg::refineHelixParameters:REMOVE] iworst=%3i (x0,y0,R) = (%8.3f, %8.3f, %8.3f) chi2 = %8.3f chi2Maxxy = %8.3f\n",
               iworst.panel, Trk._sxy.x0(), Trk._sxy.y0(), Trk._sxy.radius(), Trk._sxy.chi2DofCircle(), _chi2xyMax);//FIXME!
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

    if ((chi2_min >= _chi2xyMax) && (iworst.panel >= 0)) {
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
      // Trk._sxyw = sxyw;
      Trk._center.SetX(Trk._sxy.x0());
      Trk._center.SetY(Trk._sxy.y0());
      Trk._radius = Trk._sxy.radius();
      //      Trk._chi2   = Trk._sxy.chi2DofCircle();

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
                                     HitInfo_t           SeedIndex      ,
                                     int                 UsePhiResiduals) {

    const char  banner[] = "rescueHits";
    float      wt, /*e2,*/ x, y, r;
    float      phiwt(-9999.);

    XYZVectorF      hitPos, strawDir, helCenter, hel_pred(0.,0.,0.);

    float      dfdz, phi0, dphi, dphiChi2(0.0), phi_pred;

    ::LsqSums4  sxy;
    int         n_added_points(0);
    HitInfo_t   ibest;//(-1,-1);
    float      wtBest, phiwtBest;
    float      chi2, chi2_min, dr, hitChi2, drChi2;

    float      x_pred(0), y_pred(0), weight_hel(0), x_resid2(0), y_resid2(0);

    FaceZ_t*    facez(0);
    PanelZ_t*   panelz(0);

    mu2e::ComboHit* hit(0);

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

    ibest.face          = -1;
    ibest.panel         = -1;
    ibest.panelHitIndex = -1;

    wtBest      = -1;
    phiwtBest   = -1;

    for (int f=SeedIndex.face; f>=0; --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
      if (isFaceUsed(Helix, facez))                continue;
      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];//Helix._oTracker[p];
        int  nhits          = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhits == 0)                                       continue;
        if ((f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex  - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhits; ++i){
          int    index = panelz->idChBegin + i;
          hit = &Helix._chHitsToProcess[index];
          //          int index = facez->evalUniqueHitIndex(f,p,i);
          if (Helix._hitsUsed[index] >= 1)                    continue;

          hitPos    = hit->_pos;
          strawDir  = hit->vDir();

          dr = calculateRadialDist(hitPos,helCenter,r);
          wt = calculateWeight    (*hit,helCenter,r);

          drChi2  = (dr*dr)*wt;

          if ((UsePhiResiduals == 1) && (_phiCorrectedDefined)) {
            phi_pred = Helix._zFace[f]*dfdz + phi0;
            dphi     = phi_pred - hit->_hphi;
            phiwt    = calculatePhiWeight(*hit, helCenter, r, 0, banner);
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
            Helix._szphi.addPoint(Helix._zFace[f], hit->_hphi, phiwt);

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

                  ibest.face          = f;
                  ibest.panel         = p;
                  ibest.panelHitIndex = index;

                  wtBest      = wt;
                  phiwtBest   = phiwt;
                }
              }else {
                chi2_min    = hitChi2;

                ibest.face          = f;
                ibest.panel         = p;
                ibest.panelHitIndex = index;

                wtBest      = wt;
              }
            }

            Helix._sxy.removePoint  (x, y, wt);
            Helix._szphi.removePoint(Helix._zFace[f], hit->_hphi, phiwt);

          }
        }
      }//end panels loop
    }//ends faces loop

    if (ibest.panel >= 0){

      //evaluate the number of hits already in use in the same panel
      int    nHitsUsed(0), idSamePanel(-1);
      panelz = &facez->panelZs[ibest.panel];
      for (int l=0; l<panelz->nChHits(); ++l){
        int   id = panelz->idChBegin + l;
        if (Helix._hitsUsed[id] == 1) {
          ++nHitsUsed;
          idSamePanel = l;
        }
      }

      int    index = ibest.panelHitIndex;
      hit    = &Helix._chHitsToProcess[index];

      // now check two condition we want to skip:
      //  1) the hit we want to add belongs to a panel where we already reached the maximum number of hits allowed
      //  2) the hit belongs to the same layer where  we already reached the maximum number of hits allowed
      if (nHitsUsed >= 1){      //we don't want more than 1 ComboHit per Face
        Helix._hitsUsed[index] = 10;
                                     goto NEXT_ITERATION;
      }

      if ( idSamePanel>=0 ){
        Helix._hitsUsed[index] = 10;
                                     goto NEXT_ITERATION;
      }

      x      = hit->pos().x();
      y      = hit->pos().y();
                                       //add point from the track
      Helix._sxy.addPoint(x, y, wtBest);
      int    nHitSh = hit->nStrawHits();
      Helix._nXYSh += nHitSh;

      if (UsePhiResiduals == 1){
              Helix._szphi.addPoint(Helix._zFace[ibest.face], hit->_hphi, phiwtBest);
        Helix._nZPhiSh += nHitSh;
        dfdz  = Helix._szphi.dfdz();
        phi0  = Helix._szphi.phi0();
      }

      if (_debug > 5) {
        printf("[CalHelixFinderAlg::%s:PT2] x0 = %8.3f y0 = %8.3f radius = %8.3f  chi2 = %6.3f chi2Maxxy = %6.3f index point added = %i straw-id = %6i hitChi2 = %6.3f x = %8.3f y = %8.3f z = %9.3f\n",
               banner,
               Helix._sxy.x0(), Helix._sxy.y0(), Helix._sxy.radius(), Helix._sxy.chi2DofCircle(), _chi2xyMax, ibest.panel,
               hit->strawId().straw()/*_strawhit->strawIndex().asInt()*/, chi2_min,
               x, y, Helix._zFace[ibest.face]);//FIXME!
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
    // for (int f=SeedIndex.face; f<StrawId::_ntotalfaces; ++f){
    //   facez     = &Helix._oTracker[f];
    //   int firstPanel = 0;
    //   if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
    //   for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
    //         panelz = &facez->panelZs[p];
    //         int  nhits          = panelz->nChHits();
    //         int  seedPanelIndex(0);
    //         if (nhits == 0)                                                                    continue;
    //         if ( (f == SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

    //         for (int i=seedPanelIndex; i<nhits; ++i){
    //           int index = panelz->idChBegin + i;
    //           if (Helix._hitsUsed[index] != 1)               continue;

    //           hit     = &Helix._chHitsToProcess[index];
    //           wt      = hit->_xyWeight;
    //           e2      = 1./wt;
    //           hitPos  = hit->_pos;
    //           dr      = calculateRadialDist(hitPos,helCenter,r);
    //           hitChi2 = dr*dr/e2;
    //           // store residual
    //           //          if ( (SeedIndex.face == 0) && (SeedIndex.panel == 0) && (SeedIndex.panelHitIndex == 0)){
    //             // hit->_drFromPred = fabs(dr);
    //             // float dzFromSeed = facez->z - seedFacez->z;//seedPanelz->z;              // expected to be positive (non-negative)
    //             // hit->_dzFromSeed  = dzFromSeed;//FIXME!
    //             // }
    //         }
    //   }//end panels loop
    // }//end faces loop
//-----------------------------------------------------------------------------
// update circle parameters
//-----------------------------------------------------------------------------
    Helix._center.SetX(Helix._sxy.x0());
    Helix._center.SetY(Helix._sxy.y0());
    Helix._radius  = Helix._sxy.radius();
    //    Helix._chi2    = Helix._sxy.chi2DofCircle();

  F_END:;
    if (_debug > 5 ) {
      printf("[CalHelixFinderAlg::%s:END] N(added) = %i chi2 = %5.5f\n",banner,n_added_points,Helix._sxy.chi2DofCircle());
    }

    if ((SeedIndex.face == 0) && (SeedIndex.panel == 0) && (SeedIndex.panelHitIndex == 0)){
      Helix._diag.n_rescued_points_16 = n_added_points;
    }
  }



//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::findTrack(HitInfo_t&          SeedIndex     ,
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
    float         radius, phi0, dx, dy, phi, Chi2;
    XYZVectorF         center, shPos, hePos;
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
    mu2e::ComboHit* seedHit    = &Helix._chHitsToProcess[SeedIndex.panelHitIndex];
    float           seedFacez  = Helix._zFace[SeedIndex.face];

    //mark the seed-hit as used
    int            index      = SeedIndex.panelHitIndex;
    Helix._hitsUsed[index] = 1;

    FaceZ_t*       facez(0);
    PanelZ_t*      panelz(0);

//---------------------------------------------------------------------
// define constrains on the z coordinate of the strawhit candidate for re-calculating dfdz
// If the candidate and the seeding strawhit are too close, the dfdz calculated could be
// affected by several effects which lead to a wrong estimate
// We are asking that the candidate straw must be at a distance along
// the z axes greater than tollMin and less than tollMax.
// These parameters still need to be optimized
//-----------------------------------------------------------------------------
    float tollMin(100.) ; // , tollMax(500.);

                                        // parameters used to calculate the strawhit position residuals
    float weight(1.), wtarget(0.1);
    float deltaZ(0.), dist(0.), dist2(0.);
//----------------------------------------------------------------------//
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE //
//----------------------------------------------------------------------//
    float dfdz_end(-1e10), phi0_end(-1e10), radius_end(-1e10);

                                        // two flags are needed:
    bool removeTarget(true);            // avoid the recalculation of dfdz
                                        // and helix parameters in case when
                                        // others strawhit candidates are found
    float dfdz = _mpDfDz;               // tanLambda/radius (set to most probable);
//----------------------------------------------------------------------
// calculate helix paramters using the center of the stopping target,
// the EMC cluster which seeded the CalTimePeak and the seeding strawhit.
// The z coordinate of the target center is set to 0 because in the formula
// inside calculateTrackParameters(...) the z coordinate is not used
//-----------------------------------------------------------------------------
    float          lastFacez    = seedFacez;//Helix._zFace[SeedIndex.face];
    float          faceHitChi2  = 1e10;
    std::string name("CalHelixFinderAlg::findTrack");


    XYZVectorF p1(0.,0.,0.);               // target, z(ST) = 5971. - 10200. is not used
    XYZVectorF p2(seedHit->_pos);          // seed hit
    XYZVectorF p3(fCaloX,fCaloY,fCaloZ);   // cluster
    if( _debug > 5){
      printf("[%s]              X       Y      Z \n",name.data());
      printf("[%s]  SeedHit %6.2f %6.2f %6.2f \n",name.data(), p2.x(), p2.y(), p2.z() );
      printf("[%s]  CaloCl  %6.2f %6.2f %6.2f \n",name.data(), p3.x(), p3.y(), p3.z() );
    }
    XYZVectorF pDiff = p2-p3;
    if (!Helix._timeCluster->hasCaloCluster() && pDiff.R() < 1e-10)     return;
    if (!calculateTrackParameters(p1,p2,p3,center,radius,phi0,dfdz))    return;

//--------------------------------------------------------------------------------
// gianipez test 2019-09-28
// let's try to evaluate the dfdz NOW!
//--------------------------------------------------------------------------------
    if (_initDfDz == 1){
      int res = findDfDz(Helix, SeedIndex);
      if (res ==1 ) {
        dfdz = _hdfdz;
      }
    }
    float     tollMax = fabs(2.*M_PI/dfdz);
//------------------------------------------------------------------------------
// helix parameters, in particular, phi0, are defined at Z=p2.z()
// 2014-11-05 gianipez set dfdz equal to the most probable value for CE
//------------------------------------------------------------------------------
    if (UseMPVDfDz ==1 ) {
      dfdz    = _hdfdz;                        // _mpDfDz;
      tollMax = fabs(2.*M_PI/dfdz);
    }

    HitInfo_t lastIndex;//(-1,-1);

    ::LsqSums4 sxy;
    ::LsqSums4 szphi;

    sxy.addPoint(p2.x(), p2.y(), 1.     );  // seed hit
    sxy.addPoint(p3.x(), p3.y(), 1.     );  // EMC cluster position
    sxy.addPoint(    0.,     0., wtarget);  // Target center in the transverse plane, with small weight

    int  NPoints = seedHit->nStrawHits();     // nhits, associated with the track, sxy has NPoints+2 or NPoints+1 points
    int  NComboHits(1);

    float    z_phi0 = p2.z();

    mu2e::ComboHit* hit(0);

    for (int f=SeedIndex.face; f>=0;  --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;

      //-----------------------------------------------------------------------------
      // dfdz = tanLambda/radius; phi0 is the last found hit phi
      //-----------------------------------------------------------------------------
      deltaZ = Helix._zFace[f] - lastFacez;
      phi    = phi0 + deltaZ*dfdz;
      // phi    = phi0 - deltaZ*dfdz;
      //evaluate the helix prediction using the z coordinate of the panel
      hePos.SetXYZ(center.x()+radius*cos(phi),center.y()+radius*sin(phi),Helix._zFace[f]);//facez->z);

      //check the Panel-phi wrt to the DS center
      float  hePosPhi = polyAtan2(hePos.y(), hePos.x());
      if (hePosPhi < 0) hePosPhi = hePosPhi + 2*M_PI;

      HitInfo_t      goodFaceHit, tripletHit;

      faceHitChi2  = 1e10;

      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];//Helix._oTracker[p];
        int  nhits          = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhits == 0)                                                                        continue;

        if( _debug > 10){
          if( (f == SeedIndex.face) && (p==SeedIndex.panel) ) {
            printf("[%s:LOOP]  findTrack() starts with helix parameters derived from these points \n",name.data());
            printf("[%s:LOOP]   point  type      X         Y         Z       xyzp-index \n", name.data());
            printf("[%s:LOOP] ----------------------------------------------------------\n", name.data());
            printf("[%s:LOOP] seeding        %9.3f %9.3f %9.3f %5i\n",name.data(),p2.x(),p2.y(),p2.z(),SeedIndex.panel );
            printf("[%s:LOOP] candidate      %9.3f %9.3f %9.3f %5i\n",name.data(),p1.x(),p1.y(),p1.z(),lastIndex.panel );
            printf("[%s:LOOP] emc cluster    %9.3f %9.3f %9.3f %5i\n",name.data(),p3.x(),p3.y(),p3.z(),             -1);
            printf("[%s:LOOP]----------------------------------------------------------------------------------------------------------------------------------------\n",name.data());
            printf("[%s:LOOP]  P     Z        xi       yi       xp       yp    dXYpred  dXYseed   dZseed    X0       Y0        R        phi      dfdz    chi2    \n",name.data());
            printf("[%s:LOOP]----------------------------------------------------------------------------------------------------------------------------------------\n",name.data());
          }
        }

        if ((f == SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex + 1 - panelz->idChBegin;

        float  deltaPhi = hePosPhi - Helix._phiPanel[f*FaceZ_t::kNPanels + p];//panelz->phi;
        if ( deltaPhi > M_PI ) deltaPhi -= 2*M_PI;
        if ( deltaPhi < -M_PI) deltaPhi += 2*M_PI;
        if ( fabs(deltaPhi) > _maxPanelToHelixDPhi)                    continue;

        for (int i=seedPanelIndex; i<nhits; ++i){
          hit = &Helix._chHitsToProcess[panelz->idChBegin + i];

          // hit->_dzFromSeed = 0;
          // hit->_drFromPred = 0;

          shPos  = hit->_pos;

          // residuals in XY
          dx              = hePos.x() - hit->pos().x();
          dy              = hePos.y() - hit->pos().y();
          dist2           = dx*dx + dy*dy;
          dist            = std::sqrt(dist2);

          if( _debug > 10){
            // dist betw the straw hit and the seed hit in the transverse plane
            float dx  = std::fabs(seedHit->pos().x() - hit->pos().x());
            float dy  = std::fabs(seedHit->pos().y() - hit->pos().y());
            float dxy = std::sqrt(dx*dx+dy*dy);
            float chi2   = sxy.chi2DofCircle();

            printf("[%s:LOOP] %3i %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f %8.5f %8.3f\n",
                   name.data(),p,hit->pos().z(),hit->pos().x(),hit->pos().y(),hePos.x(),hePos.y(),dist,dxy,deltaZ,center.x(),center.y(),radius,phi,dfdz,chi2) ;
          }
          //-----------------------------------------------------------------------------
          // dxy_max: running search window accounts for the finite extrapolation accuracy
          //-----------------------------------------------------------------------------
          //          float dxy_max = _distPatRec + _dfdzErr*deltaZ;
          float dxy_max = _distPatRec + _dfdzErr*std::fabs(deltaZ);
          if (dist <= dxy_max) {
            if (dist < faceHitChi2){
              faceHitChi2               = dist;
              //set the coordiante of the hit
              goodFaceHit.face          = f;
              goodFaceHit.panel         = p;
              goodFaceHit.panelHitIndex = panelz->idChBegin + i;
            }
          }
        }

      }//end panels loop

      //check if a hit close to the predixction was found
      if (goodFaceHit.face < 0)                 continue;

      hit       = &Helix._chHitsToProcess[goodFaceHit.panelHitIndex];

      //mark the hit as used
      index = goodFaceHit.panelHitIndex;  //facez->evalUniqueHitIndex(goodFaceHit);
      Helix._hitsUsed[index] = 1;

      ++NComboHits;
      NPoints += hit->nStrawHits();

      //-----------------------------------------------------------------------------
      // Mode = 0: helix parameters evaluated using stopping_target+fitst_hit(seed)+cluster
      //           dphi/dz is fixed and so far set to the most proable value for CE
      //-----------------------------------------------------------------------------
      sxy.addPoint(hit->pos().x(),hit->pos().y(),weight);
      if (Mode == 1) {
        // dfdz has already been evaluated, update XY part of the helix
        center.SetX(sxy.x0());
        center.SetY(sxy.y0());
        radius = sxy.radius();
      }

      phi0      = polyAtan2(hit->pos().y()-center.y(),hit->pos().x()-center.x());  // *FLOAT_CHECK*
      z_phi0    = Helix._zFace[goodFaceHit.face];
      lastFacez = z_phi0;//Helix._zFace[goodFaceHit.face];

      if      ( Mode == 0 ) ++mode0GoodPoints;
      //      else if ((Mode == 1) && (panelHitIndex[t]<= lastIndex.PanelHitIndex)) ++mode1GoodPoints;//FIXME!

      float dzFromSeed = -1.*(z_phi0 - seedFacez);
      // hit->_dzFromSeed = dzFromSeed;
      // hit->_drFromPred = faceHitChi2;//panelHitChi2[t];

      if (/*(panelHitChi2[t] < dist_min) &&*/ (dzFromSeed > tollMin) && (dzFromSeed < tollMax)) {
        //-----------------------------------------------------------------------------
        // goodPoint - index of the first hit separated in Z from the seed by > 10 cm
        // after such found, the target center is removed and the circle parameters
        // recalculated using the cluster, the seed hit and the 'goodPoint' hit
        // an additional requirement is that at the recalculation time there are 3 or more
        // hits found in total
        //-----------------------------------------------------------------------------
        if (removeTarget) {
          //            goodPoint = panelHitIndex[t];
          tripletHit = goodFaceHit;
          //          dist_min   = faceHitChi2;//panelHitChi2 [t];
        }
      }


      if ((tripletHit.face >= 0) && (NComboHits >= 2) && (UseMPVDfDz == 0)) {
        //-----------------------------------------------------------------------------
        // the first point separated from the seed one by more than 10 cm has been found
        // recalculate helix parameters: for XY part use accumulated sxy sums
        // replace stopping target with the hit
        //-----------------------------------------------------------------------------
        //        sxy.removePoint(0.,0.,wtarget);
        //        hit = &panelz->_chHitsToProcess.at(goodPoint);
        p1  = hit->_pos;

        center.SetX(sxy.x0());
        center.SetY(sxy.y0());
        radius = sxy.radius();
        //-----------------------------------------------------------------------------
        // now calculate more accuratelly the value of dfdz using just the two strawhit positions
        // change in the circle parameterization changes the phi0 value
        //-----------------------------------------------------------------------------
        phi0 = polyAtan2(hit->pos().y()-center.y(),hit->pos().x()-center.x());

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
          printf("[%s:DEF2] %3i %9.3f %8.3f %8.3f seed \n", name.data(),SeedIndex.panel,seedHit->pos().z(),seedHit->pos().x(),seedHit->pos().y());//FIXME!
          printf("[%s:DEF2] %3i %9.3f %8.3f %8.3f seed \n", name.data(), tripletHit.panel,hit->pos().z()    ,hit->pos().x()    ,hit->pos().y()    );
        }
        //-----------------------------------------------------------------------------
        // what to do if dfdz is negative? - the case of negative helicity is not covered yet
        //-----------------------------------------------------------------------------
        if ((fabs(dfdz) > _maxDfDz) || (fabs(dfdz) < _minDfDz)) {
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

          lastIndex    = tripletHit;
        }
      }

    }//end face loop

    if (NPoints < _minNHits) return;
    //-----------------------------------------------------------------------------
    // 3 or more points have been found on a helix candidate, update Chi2
    //-----------------------------------------------------------------------------
    Chi2 = sxy.chi2DofCircle();
    //-----------------------------------------------------------------------------
    // temporary variables to store dfdz values out of the method 'calculateDfDz(...)'
    //-----------------------------------------------------------------------------
    float dfdzRes  [3] = {   -1.,    -1.,    -1.};
    float dphi0Res [3] = {-9999., -9999., -9999.};
    float radiusRes[2] = {   -1.,    -1.};

    if (_diag > 0) {
      if (UseMPVDfDz == 0) dfdzRes[0] = dfdz;
      dphi0Res [0] = phi0 - dfdz*z_phi0;
      radiusRes[0] = sxy.radius();
    }
//-----------------------------------------------------------------------------
// initialize only the xy part, z-phi part is not needed here
//-----------------------------------------------------------------------------
    Helix._sxy = sxy;
    Helix._nXYSh      = NPoints;
    Helix._radius     = sxy.radius();
    Helix._center.SetXYZ(sxy.x0(), sxy.y0(), 0.0);
    Helix._nStrawHits = NPoints;
    Helix._nComboHits = NComboHits;
    Helix._dfdz       = dfdz;
    Helix._fz0        = phi0 - dfdz*z_phi0; // *FLOAT_CHECK*

    radius_end = Helix._radius;
    //breakpoint -- radius before refine $
    int rc = refineHelixParameters(Helix, SeedIndex);
//-----------------------------------------------------------------------------
// if weighted XY fit didn't converge, there is nothing else one can do, return
//-----------------------------------------------------------------------------
    if (rc < 0) return;
    //breakpoint -- radius after refine $
    // Helix._center.SetXYZ(Helix._cw.x(), Helix._cw.y(), 0.0);
    // Helix._radius  = Helix._rw;
    radius_end     = Helix._radius;

    // Helix._sxy = Helix._sxyw;
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
      dfdz_end    = Helix._dfdz;
      phi0_end    = Helix._fz0;
                                        // fill diagnostic vector
      dfdzRes [2] = Helix._dfdz;
      dphi0Res[2] = Helix._fz0;

      NPoints     = 0;
      NComboHits  = 0;

      //FIXME! implement a function
      countUsedHits(Helix, SeedIndex, NComboHits, NPoints);

      // for (int f=SeedIndex.face; f<StrawId::_ntotalfaces;  ++f){
      //         facez     = &Helix._oTracker[f];
      //         int  firstPanel(0);
      //         if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
      //         for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
      //           panelz = &facez->panelZs[p];
      //           int  nhitsPerPanel  = panelz->nChHits();
      //           int  seedPanelIndex(0);
      //           if (nhitsPerPanel == 0)                                 continue;
      //           if ( (f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

      //           for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
      //             index = panelz->idChBegin + i;
      //             hit   = &Helix._chHitsToProcess[index];
      //             if (Helix._hitsUsed[index] > 0 )  {
      //               ++NComboHits;
      //               NPoints += hit->nStrawHits();
      //             }
      //           }
      //         }//endl panels loop
      // }
    }
    else {
      dfdz_end = _hdfdz;
      phi0_end = _hphi0;
    }

    if (_debug > 10) {
      printf("[%s] strawhit type     X        Y        Z     index     \n", name.data());
      printf("[%s] ----------------------------------------------------\n", name.data());
      printf("[%s]    seeding   %9.3f %9.3f %9.3f   %i  \n", name.data(),p2.x(), p2.y(), p2.z(), SeedIndex.panel);//FIXME!
      printf("[%s]   candidate  %9.3f %9.3f %9.3f   %i  \n", name.data(),p1.x(), p1.y(), p1.z(), lastIndex.panel);//FIXME!
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
             name.data(),SeedIndex.panel,Mode,UseMPVDfDz,NPoints,Chi2,dfdz_end);//FIXME!
    }

    // can execution really come here with Mode == 0 ? - YES! not sure, why
    //      if (( NPoints >  _goodPointsTrkCandidate) ||
    //           ((NPoints == _goodPointsTrkCandidate) && (Chi2< _chi2TrkCandidate))) {
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

      float dz = p1.z() - p2.z();

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

      for (int f=SeedIndex.face; f>=0;  --f){
        facez     = &Helix._oTracker[f];
        int  firstPanel(0);
        if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
        for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
          panelz = &facez->panelZs[p];//Helix._oTracker[p];
          int  nhitsPerPanel  = panelz->nChHits();
          int  seedPanelIndex(0);
          if (nhitsPerPanel == 0)                          continue;
          if ( (f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0)) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

          for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
            int   index = panelz->idChBegin + i;
            hit = &Helix._chHitsToProcess[index];

            //            int index = facez->evalUniqueHitIndex(f,p,i);//p*FaceZ_t::kNMaxHitsPerPanel + i;
            if (Helix._hitsUsed[index] != 1)               continue;

            //  if (j < Helix.maxIndex()) {
              // Helix._diag.dist[j] = hit->_drFromPred;
              // Helix._diag.dz  [j] = hit->_dzFromSeed;
            //   ++j;
            // }
            // else {
            //   printf("ERROR in CalHelixFinderAlg::findTrack : index out limits. IGNORE; \n");
            // }
          }//end loop over the hits within a panel
        }//end panel loop
      }//end face loop
    }

  }

  //-----------------------------------------------------------------------------
  // helix parameters are defined at Z=p2.z, Phi0 corresponds to p2
  //-----------------------------------------------------------------------------
  bool CalHelixFinderAlg::calculateTrackParameters(const XYZVectorF&   p1       ,
                                                   const XYZVectorF&   p2       ,
                                                   const XYZVectorF&   p3       ,
                                                   XYZVectorF&         Center   ,
                                                   float&         Radius   ,
                                                   float&         Phi0     ,
                                                   float&         DfDz32) {
    //evaluate the area covered by the Triplet
    float dist2ij = (p1 - p2).Mag2();
    float dist2ik = (p1 - p3).Mag2();
    float dist2jk = (p2 - p3).Mag2();
    float area2   = (dist2ij*dist2jk + dist2ik*dist2jk + dist2ij*dist2ik) - 0.5*(dist2ij*dist2ij + dist2jk*dist2jk + dist2ik*dist2ik);
    if(area2 < _minarea2)              return false;


    Center.SetZ(p2.z());

    float x_m, y_m, x_n, y_n;
    // coordinates of the mean point between p1 and p3
    x_m = (p3.x() + p1.x())/2.;
    y_m = (p3.y() + p1.y())/2.;
    // coordinates of the mean point between p2 and p3
    x_n = (p3.x() + p2.x())/2.;
    y_n = (p3.y() + p2.y())/2.;
    //------------------------------------------------------------//

    //calculate now the term of the line ortoghonal to the mid point of
    //the cord which links p1 and p3
    float m = -1.*(p3.x() - p1.x())/(p3.y() - p1.y());
    float c = y_m - x_m*m;
    //the eq. is: y = x*m + c

    //calculate now the term of the line ortoghonal to the mid point of
    //the cord which links p2 and p3
    float k = -1.*(p3.x() - p2.x())/(p3.y() - p2.y());
    float t = y_n - x_n*k;
    //the eq. is: y = x*k + t

    if (std::isfinite(m) && std::isfinite(k) && std::fabs(k) > 1e-9) {
      if ( (m/k>0) && ( (m/k) - int(m/k) > _slopeRatioLimit) ) {//invert p3 with p1 and recalculate: x_n, y_n, k, t
        x_n = (p1.x() + p2.x())/2.;
        y_n = (p1.y() + p2.y())/2.;
        k   = -1.*(p1.x() - p2.x())/(p1.y() - p2.y());
        t   = y_n - x_n*k;
      }
    }

    if (!std::isfinite(m) || !std::isfinite(k) || std::fabs(m-k) < 1e-6) return false;

    // calculate Center.x and Center.y
    float x0 = (t - c)/(m - k);//(c - t) * (k*m)/(m-k);
    Center.SetX(x0);
    float y0 = m*x0 + c;   //(c - t) * m / (m - k) + t;
    Center.SetY(y0);

//-----------------------------------------------------------------------------
// calculate the radius,phi0, tanLambda assuming that the helix also crosses
// the point (0,0). Note that the Z-position of the stopping target is not used
//-----------------------------------------------------------------------------
    float dx3  = p3.x() - x0;
    float dy3  = p3.y() - y0;
    float dz32 = p3.z() - p2.z();

    Radius      = std::sqrt(dx3*dx3+dy3*dy3);

    float dx2  = (p2.x() - x0);
    float dy2  = (p2.y() - y0);

    Phi0        = polyAtan2(dy2,dx2);
//-----------------------------------------------------------------------------
// this assumes that the helix is right-handed, *FIXME*
// make sure that we are lookign for a particle which makes close to expected
// number of turns
//-----------------------------------------------------------------------------
    float dphi32 = polyAtan2(dy3,dx3) - Phi0;
    if (dphi32*_dfdzsign < 0.) dphi32 += 2.*M_PI;

    //    float exp_dphi = _mpDfDz*dz32;

    //check id DfDz is within the range
    if ( (fabs(DfDz32) < _minDfDz) || (fabs(DfDz32) > _maxDfDz)) DfDz32 = _mpDfDz;

    DfDz32 = dphi32/dz32;

    float   diff      = fabs(DfDz32 - _mpDfDz);
    float   diff_plus = fabs( (dphi32 + 2.*M_PI)/dz32 -_mpDfDz );
    while ( diff_plus < diff ){
      dphi32  = dphi32 + 2.*M_PI;
      DfDz32      = dphi32/dz32;
      diff      = fabs(DfDz32 - _mpDfDz);
      diff_plus = fabs( (dphi32 + 2.*M_PI)/dz32 -_mpDfDz );
    }

    float   diff_minus = fabs( (dphi32 - 2.*M_PI)/dz32 -_mpDfDz );
    while ( diff_minus < diff ){
      dphi32   = dphi32 - 2.*M_PI;
      DfDz32       = dphi32/dz32;
      diff       = fabs(DfDz32 - _mpDfDz);
      diff_minus = fabs( (dphi32 - 2.*M_PI)/dz32 -_mpDfDz );
    }

    //check id DfDz is within the range
    if ( (fabs(DfDz32) < _minDfDz) || (fabs(DfDz32) > _maxDfDz)) DfDz32 = _mpDfDz;

    if (_debug > 5) {
//-----------------------------------------------------------------------------
// in debug mode also want to print the helix parameters, calculate them
// phi00 - phi angle of (x0,y0) point - center of the helix
//-----------------------------------------------------------------------------
      float d0     = sqrt(x0*x0+y0*y0)-Radius;
      float phi00  = polyAtan2(y0,x0);                    // sign taken into account
      float tandip = DfDz32*_dfdzsign*Radius;             // signs of DfDz32 and _dfdzsign should be the same
      float dphi   = phi00-Phi0;
      if (dphi < 0) dphi += 2*M_PI;                        // *FIXME* right-handed helix

      float z0     = p2.z()-dphi*dz32/dphi32;

      printf("[CalHelixFinderAlg:calculateTrackParameters] X0: %9.3f Y0: %9.3f phi0: %8.5f p1.z = %9.3f p2.z = %9.3f p3.z = %9.3f dphi32 = %8.5f dfdz = %8.5f\n",
             Center.x(),Center.y(),Phi0,p1.z(),p2.z(),p3.z(),dphi32,DfDz32);
      printf("[CalHelixFinderAlg:calculateTrackParameters] z0 = %9.3f d0 = %8.4f  phi00 = %8.5f omega = %8.5f tandip = %8.4f\n",z0,d0,phi00,1/Radius,tandip);
    }

    return true;
  }

//-----------------------------------------------------------------------------
// the function is currently not called
//-----------------------------------------------------------------------------
  void  CalHelixFinderAlg::calculateDfDz(float phi0, float phi1, float z0, float z1, float& DfDz) {
    float   deltaPhi  = TVector2::Phi_mpi_pi(phi1-phi0);
    DfDz               = deltaPhi/(z1-z0);

    // 2018-01-02: don't do that!
    // float   diff      = fabs(DfDz - _mpDfDz);
    // float   diff_plus = fabs((deltaPhi + 2.*M_PI)/(z1-z0) -_mpDfDz);
    // while (diff_plus < diff) {
    //   deltaPhi  = deltaPhi + 2.*M_PI;
    //   DfDz      = deltaPhi/(z1-z0);
    //   diff      = fabs(DfDz - _mpDfDz);
    //   diff_plus = fabs( (deltaPhi + 2.*M_PI)/(z1-z0) -_mpDfDz );
    // }

    // float   diff_minus = fabs((deltaPhi - 2.*M_PI)/(z1-z0) -_mpDfDz);
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
  void  CalHelixFinderAlg::calculateDphiDz_2(CalHelixFinderData& Helix, HitInfo_t SeedIndex,
                                             int NHits, float X0, float Y0, float& DphiDz) {
    LsqSums2       sphiz;

    float         phi(0), phi0(0);//, phiCl(0);

    ComboHit*      hit(0);
    FaceZ_t*       facez(0);
    PanelZ_t*      panelz(0);
    XYZVectorF*        pos(0);

    int            counter(0);

    bool           isFirst(true);

    for (int f=SeedIndex.face; f>=0;  --f){
      facez     = &Helix._oTracker[f];
      int  firstPanel(0);
      if (f == SeedIndex.face) firstPanel = SeedIndex.panel;
      for (int p=firstPanel; p<FaceZ_t::kNPanels; ++p){
        if (counter > NHits  )                               break;
        panelz = &facez->panelZs[p];//Helix._oTracker[p];
        int  nhitsPerPanel  = panelz->nChHits();
        int  seedPanelIndex(0);
        if (nhitsPerPanel == 0)                           continue;
        if ( (f==SeedIndex.face) && (p==SeedIndex.panel) && (SeedIndex.panelHitIndex >=0) ) seedPanelIndex = SeedIndex.panelHitIndex - panelz->idChBegin;

        for (int i=seedPanelIndex; i<nhitsPerPanel; ++i){
          int   index = panelz->idChBegin + i;
          hit = &Helix._chHitsToProcess[index];
          //          int index = facez->evalUniqueHitIndex(f,p,i);//p*FaceZ_t::kNMaxHitsPerPanel + i;
          if (Helix._hitsUsed[index] != 1)                     continue;

          pos = &(hit->_pos);
          phi = polyAtan2(pos->y()-Y0,pos->x()-X0);
          if (isFirst) {
            phi0 = phi;
            isFirst = false;
          }

          if (phi-phi0 >  M_PI) phi -= 2*M_PI;
          if (phi-phi0 < -M_PI) phi += 2*M_PI;

          sphiz.addPoint(pos->z(),phi);
          ++counter;

          if (_debug > 10) printf("[CalHelixFinderAlg::calculateDphiDz_2:LOOP] panel,id,phi,z=%3i %3i %8.5f %9.3f\n",p,i,pos->z(),phi);
        }
      }//end panel loop
    }//end face loop
//-----------------------------------------------------------------------------
// define straight line phi = phi0+dPhi/Dz*z , where phi0 = phi(z=0)
//-----------------------------------------------------------------------------
    phi0   = sphiz.yMean();
    DphiDz = sphiz.dydx();

    if (_debug > 5) printf("[CalHelixFinderAlg::calculateDphiDz_2:END] phi0,DphiDz = %9.5f %9.5f \n",phi0,DphiDz);
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
    int n = Helix._nFiltPoints;

    printf("[CalHelixFinderAlg::printXYZP]-----------------------------------------------------------------------------------------\n");
    printf("[CalHelixFinderAlg::printXYZP]     f     p     i    shId    index     X         Y         Z      Z-Face    _debug: %5i nhits: %5i\n",_debug,n);
    printf("[CalHelixFinderAlg::printXYZP]-----------------------------------------------------------------------------------------\n");

    FaceZ_t*  facez(0);
    PanelZ_t* panelz(0);

    for (int f=StrawId::_ntotalfaces-1; f>=0;  --f){
      facez     = &Helix._oTracker[f];
      for (int p=0; p<FaceZ_t::kNPanels; ++p){
        panelz = &facez->panelZs[p];//Helix._oTracker[p];
        int  nhitsPerPanel  = panelz->nChHits();


        for (int i=0; i<nhitsPerPanel; ++i){
          int             index = panelz->idChBegin + i;
          mu2e::ComboHit* pt    = &Helix._chHitsToProcess[index];//panelz->_chHitsToProcess.at(i);

          printf("[CalHelixFinderAlg::printXYZP] %5i %5i %5i   %08x   %2i %9.3f %9.3f %9.3f  %9.3f \n",
                 f, p, i, pt->strawId().straw(), index, pt->_pos.x(), pt->_pos.y(), pt->_pos.z(), Helix._zFace[f]);
        }
      }
    }
  }

  bool CalHelixFinderAlg::isFaceUsed(CalHelixFinderData& Helix, FaceZ_t* facez){
    int c = 0;
    for (int p = 0; p < FaceZ_t::kNPanels; p++){
      PanelZ_t* panelz = &facez->panelZs[p];
      int  nhits = panelz->nChHits();
      for (int i=0; i<nhits; ++i){
        int index = panelz->idChBegin + i;
        if (Helix._hitsUsed[index] == 1){
          c++;
          break;
        }
      }
      if (c >= 1) break;
    }

    if (c > 0) return true;
    else       return false;
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void CalHelixFinderAlg::resolve2PiAmbiguity(ComboHit* Hit, XYZVectorF& Center, float &Phi_ref, float &DPhi){
    float dx      = (Hit->_pos.x() - Center.x());
    float dy      = (Hit->_pos.y() - Center.y());
    float phi     = polyAtan2(dy, dx);
    if (phi < 0) phi = phi + 2*M_PI;
    DPhi    = Phi_ref - phi;
    // resolve 2PI ambiguity
    while (DPhi > M_PI) {
      phi += 2*M_PI;
      DPhi = Phi_ref - phi;
    }
    while (DPhi < -M_PI) {
      phi -= 2*M_PI;
      DPhi = Phi_ref - phi;
    }
    // store the corrected value of phi
    Hit->_hphi = phi;

  }
  // void CalHelixFinderAlg::resolve2PiAmbiguity(CalHelixFinderData& Helix,const XYZVectorF& Center, float DfDz, float Phi0){

  //   const XYZVectorF*  pos;
  //   float         z, phi, phi_ref, dphi, dx, dy;

  //   FaceZ_t*        facez(0);
  //   PanelZ_t*       panelz(0);
  //   mu2e::ComboHit* hit(0);

  //   for (int f=0; f<StrawId::_ntotalfaces;  ++f){
  //     facez     = &Helix._oTracker[f];
  //     for (int p=0; p<FaceZ_t::kNPanels; ++p){
  //         panelz = &facez->panelZs[p];//Helix._oTracker[p];
  //         int  nhits          = panelz->nChHits();
  //         for (int i=0; i<nhits; ++i){
  //           hit = &panelz->_chHitsToProcess.at(i);
  //           pos = &hit->_pos;
  //           z   = pos->z();

  //           dx  = (pos->x() - Center.x());
  //           dy  = (pos->y() - Center.y());
  //           phi = polyAtan2(dy, dx);
  //           if (phi < 0) phi = phi + 2*M_PI;//phi = TVector2::Phi_0_2pi(phi);
  //           // predicted value of phi
  //           phi_ref = z*DfDz + Phi0;
  //           // signed residual
  //           dphi    = phi_ref - phi;
  //           // resolve the 2PI ambiguity
  //           while (dphi > M_PI) {
  //             phi += 2*M_PI;
  //             dphi = phi_ref - phi;
  //           }
  //           while (dphi < -M_PI) {
  //             phi -= 2*M_PI;
  //             dphi = phi_ref - phi;
  //           }
  //           // store the corrected value of phi
  //           // _phiCorrected[i] = phi;
  //         }

  //     }
  //   }
  //                 // don't know
  //   _phiCorrectedDefined = 0;
  // }



//---------------------------------------------------------------------------
// reset track paramters
// indices on the xyzp vector of: the straw hit seeding the search,
// the second strawhit used for recalculating the dfdz value
//---------------------------------------------------------------------------
  void CalHelixFinderAlg::resetTrackParamters() {

    fUseDefaultDfDz = 0;
    _hphi0                  = -9999.;
//-----------------------------------------------------------------------------
// quality paramters used for doing comparison between several track candidates
//-----------------------------------------------------------------------------
    _hdfdz                  = _mpDfDz;
  }

}
