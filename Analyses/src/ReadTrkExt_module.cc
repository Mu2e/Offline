//
//  Track Extrapolation Reader module
//  All in Detector (=Tracker) coordinate. Points other than tracker coordinates are properly transformed. 
//  B field basically requires mu2e coordinate
//
//  $Id: ReadTrkExt_module.cc,v 1.3 2012/08/31 22:36:33 brownd Exp $
//  $Author: brownd $
//  $Date: 2012/08/31 22:36:33 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/WorldG4.hh"

using namespace CLHEP;
#include "CLHEP/Vector/ThreeVector.h"

#include "TH1F.h"
#include "TTree.h"

#include "KalmanTests/inc/KalRepCollection.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/PointTrajectory.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "TrkBase/TrkParticle.hh"

#include "TrkExt/inc/TrkExtMCHits.hh"
#include "TrkExt/inc/TrkExtDetectors.hh"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"


using namespace std;

namespace mu2e {

  const unsigned int MAXNHOT = 100;
  const unsigned int MAXNTRK = 100;
  const unsigned int MAXNSIM = 5000;
  const unsigned int MAXNBACK = 5000;
  const unsigned int MAXNPAHITS = 50;

  class ReadTrkExt : public art::EDAnalyzer {

  public:
    explicit ReadTrkExt(fhicl::ParameterSet const& pset);
    ~ReadTrkExt() { }
    void beginJob();
    void beginRun(art::Run const &run);
    void beginSubRun(art::SubRun const& lblock );
    void analyze(const art::Event& event);
    void endJob();

  private:

    std::string _fitterModuleLabel;
    std::string _trkextModuleLabel;
    std::string _g4ModuleLabel;
    std::string _makerModuleLabel;
    bool _mcFlag;
    double _paClusterAcceptance;
    int _debugLevel;

    bool _flagPAHits;
    bool _flagSim;
    bool _flagPatRec;
    bool _flagHits;

    string _trkPatRecInstanceName;

    TTree * _hTrk;

    int _evtid, _trkid, _nhits, _nhots, _ntrks, _nsim, _nback, _nfwrd, _npahits, _npaclust, _nbackpaclust, _nbackstclust, _exitcode;
    float _hotx[MAXNHOT], _hoty[MAXNHOT], _hotz[MAXNHOT], _hott0[MAXNHOT];
    float _trkl0, _trkl1, _trkl[MAXNTRK], _trkx[MAXNTRK], _trky[MAXNTRK], _trkz[MAXNTRK];
    float _trkd0[MAXNTRK], _trkz0[MAXNTRK], _trkphi[MAXNTRK], _trkomega[MAXNTRK], _trktanDip[MAXNTRK];
    float _trkpx[MAXNTRK], _trkpy[MAXNTRK], _trkpz[MAXNTRK], _trkp[MAXNTRK];
    float _simx[MAXNSIM], _simy[MAXNSIM], _simz[MAXNSIM], _simp0, _simt0;
    float _vdx[5], _vdy[5], _vdz[5], _vdpx[5], _vdpy[5], _vdpz[5], _vdp[5];
    float _back_x[MAXNBACK], _back_y[MAXNBACK], _back_z[MAXNBACK], _back_px[MAXNBACK], _back_py[MAXNBACK], _back_pz[MAXNBACK], _back_p[MAXNBACK], _back_rho[MAXNBACK], _back_s, _back_t;
    float _back_ex[MAXNBACK], _back_ey[MAXNBACK], _back_ez[MAXNBACK], _back_epx[MAXNBACK], _back_epy[MAXNBACK], _back_epz[MAXNBACK], _back_ep[MAXNBACK], _back_er[MAXNBACK]; 
    int _back_vid[MAXNBACK];
    float _back_paclust_z[MAXNPAHITS], _back_paclust_dp[MAXNPAHITS], _back_paclust_dptot;
    float _back_stclust_z[MAXNPAHITS], _back_stclust_dp[MAXNPAHITS], _back_stclust_dptot;
    float _pa_z[MAXNPAHITS];
    int _paclust_nhits[MAXNPAHITS];
    float _paclust_px[MAXNPAHITS], _paclust_py[MAXNPAHITS], _paclust_pz[MAXNPAHITS], _paclust_p[MAXNPAHITS], _paclust_edep[MAXNPAHITS], _paclust_z[MAXNPAHITS]; 
    
    Hep3Vector _origin;
    Hep3Vector _mu2eOriginInWorld;

    void readHits ( const art::Event& event, TrkHotList const * hits) ;
    void readTrkPatRec(const KalRep & kalrep);
    int  readSimulation(const art::Event& event, const mu2e::TrkStrawHit* trkstrawHit) ;
    void doPAHitBooking ( vector<TrkExtTrajPoint> & hits) ;
    void doSTHitBooking ( vector<TrkExtTrajPoint> & hits) ;
    void doExtDataBooking (const TrkExtTraj & traj, vector<TrkExtTrajPoint> & pahits, vector<TrkExtTrajPoint> & sthits) ;

  };

  ReadTrkExt::ReadTrkExt(fhicl::ParameterSet const& pset):
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _trkextModuleLabel(pset.get<string>("trkextModuleLabel")),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _mcFlag(pset.get<bool>("mcFlag", false)),
    _paClusterAcceptance(pset.get<double>("paClusterAcceptance", 3)),
    _debugLevel(pset.get<int>("debugLevel", 1)),
    _hTrk(0)
  {

    _nhits = 0;
    _nhots = 0;
    _ntrks = 0;
    _nsim = 0;
    _nback = 0;
    switch (_debugLevel) {
      case 1:
        _flagPAHits = false;
        _flagSim = true;
        _flagPatRec = true;
        _flagHits = true;
        break;
      case 2:
        _flagPAHits = true;
        _flagSim = true;
        _flagPatRec = true;
        _flagHits = true;
        break;
      case 0:
      default:
        _flagPAHits = false;
        _flagSim = false;
        _flagPatRec = false;
        _flagHits = false;
    }

    _trkPatRecInstanceName = TrkFitDirection(TrkFitDirection::downstream).name() + TrkParticle(TrkParticle::e_minus).name();

  }

  void ReadTrkExt::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;

    // _hTrk
    _hTrk = tfs->make<TTree>("hTrk", "Track and recon info");
    _hTrk->Branch("evtid", &_evtid, "evtid/I");
    _hTrk->Branch("trkid", &_trkid, "trkid/I");
    _hTrk->Branch("nhits", &_nhits, "nhits/I");
    _hTrk->Branch("nhots", &_nhots, "nhots/I");
    _hTrk->Branch("ntrks", &_ntrks, "ntrks/I");
    _hTrk->Branch("exitcode", &_exitcode, "exitcode/I");
    if (_flagHits) {
      _hTrk->Branch("hotx", _hotx, "hotx[nhots]/F");
      _hTrk->Branch("hoty", _hoty, "hoty[nhots]/F");
      _hTrk->Branch("hotz", _hotz, "hotz[nhots]/F");
      _hTrk->Branch("hott0", _hott0, "hott0[nhots]/F");
    }
    if (_flagPatRec) {
      _hTrk->Branch("trkl0", &_trkl0, "trkl0/F");
      _hTrk->Branch("trkl1", &_trkl1, "trkl1/F");
      _hTrk->Branch("trkl", _trkl, "trkl[ntrks]/F");
      _hTrk->Branch("trkx", _trkx, "trkx[ntrks]/F");
      _hTrk->Branch("trky", _trky, "trky[ntrks]/F");
      _hTrk->Branch("trkz", _trkz, "trkz[ntrks]/F");
      _hTrk->Branch("trkd0", _trkd0, "trkd0[ntrks]/F");
      _hTrk->Branch("trkz0", _trkz0, "trkz0[ntrks]/F");
      _hTrk->Branch("trkphi", _trkphi, "trkphi0[ntrks]/F");
      _hTrk->Branch("trkomega", _trkomega, "trkomega[ntrks]/F");
      _hTrk->Branch("trktanDip", _trktanDip, "trktanDip[ntrks]/F");
      _hTrk->Branch("trkp", _trkp, "trkp[ntrks]/F");
      _hTrk->Branch("trkpx", _trkpx, "trkpx[ntrks]/F");
      _hTrk->Branch("trkpy", _trkpy, "trkpy[ntrks]/F");
      _hTrk->Branch("trkpz", _trkpz, "trkpz[ntrks]/F");
    }
    _hTrk->Branch("nback", &_nback, "nback/I");
    _hTrk->Branch("back_x", _back_x, "back_x[nback]/F");
    _hTrk->Branch("back_y", _back_y, "back_y[nback]/F");
    _hTrk->Branch("back_z", _back_z, "back_z[nback]/F");
    _hTrk->Branch("back_p", _back_p, "back_p[nback]/F");
    _hTrk->Branch("back_rho", _back_rho, "back_rho[nback]/F");
    _hTrk->Branch("back_px", _back_px, "back_px[nback]/F");
    _hTrk->Branch("back_py", _back_py, "back_py[nback]/F");
    _hTrk->Branch("back_pz", _back_pz, "back_pz[nback]/F");
    _hTrk->Branch("back_ex", _back_ex, "back_ex[nback]/F");
    _hTrk->Branch("back_ey", _back_ey, "back_ey[nback]/F");
    _hTrk->Branch("back_ez", _back_ez, "back_ez[nback]/F");
    _hTrk->Branch("back_epx", _back_epx, "back_epx[nback]/F");
    _hTrk->Branch("back_epy", _back_epy, "back_epy[nback]/F");
    _hTrk->Branch("back_epz", _back_epz, "back_epz[nback]/F");
    _hTrk->Branch("back_ep", _back_ep, "back_ep[nback]/F");
    _hTrk->Branch("back_er", _back_er, "back_er[nback]/F");
    _hTrk->Branch("back_s", &_back_s, "back_s/F");
    _hTrk->Branch("back_t", &_back_t, "back_t/F");
    _hTrk->Branch("back_vid", _back_vid, "back_vid[nback]/I");
    _hTrk->Branch("nbackpaclust", &_nbackpaclust, "nbackpaclust/I");
    _hTrk->Branch("back_paclust_z", _back_paclust_z, "back_paclust_z[nbackpaclust]/F"); 
    _hTrk->Branch("back_paclust_dp", _back_paclust_dp, "back_paclust_dp[nbackpaclust]/F"); 
    _hTrk->Branch("back_paclust_dptot", &_back_paclust_dptot, "back_paclust_dptot/F"); 
    _hTrk->Branch("nbackstclust", &_nbackstclust, "nbackstclust/I");
    _hTrk->Branch("back_stclust_z", _back_stclust_z, "back_stclust_z[nbackstclust]/F"); 
    _hTrk->Branch("back_stclust_dp", _back_stclust_dp, "back_stclust_dp[nbackstclust]/F"); 
    _hTrk->Branch("back_stclust_dptot", &_back_stclust_dptot, "back_stclust_dptot/F"); 
    if (_mcFlag) {
      if (_flagSim) {
        _hTrk->Branch("nsim", &_nsim, "nsim/I");
        _hTrk->Branch("simx", _simx, "simx[nsim]/F");
        _hTrk->Branch("simy", _simy, "simy[nsim]/F");
        _hTrk->Branch("simz", _simz, "simz[nsim]/F");
        _hTrk->Branch("simp0", &_simp0, "simp0/F");
        _hTrk->Branch("simt0", &_simt0, "simt0/F");
        _hTrk->Branch("vdx", _vdx, "vdx[5]/F");
        _hTrk->Branch("vdy", _vdy, "vdy[5]/F");
        _hTrk->Branch("vdz", _vdz, "vdz[5]/F");
        _hTrk->Branch("vdpx", _vdpx, "vdpx[5]/F");
        _hTrk->Branch("vdpy", _vdpy, "vdpy[5]/F");
        _hTrk->Branch("vdpz", _vdpz, "vdpz[5]/F");
        _hTrk->Branch("vdp", _vdp, "vdp[5]/F");
      }
      if (_flagPAHits) {
        _hTrk->Branch("npahits", &_npahits, "npahits/I");
        _hTrk->Branch("pa_z", _pa_z, "pa_z[npahits]/F");
        _hTrk->Branch("npaclust", &_npaclust, "npaclust/I");
        _hTrk->Branch("paclust_nhits", _paclust_nhits, "paclust_nhits[npaclust]/I");
        _hTrk->Branch("paclust_px", _paclust_px, "paclust_px[npaclust]/F");
        _hTrk->Branch("paclust_py", _paclust_py, "paclust_py[npaclust]/F");
        _hTrk->Branch("paclust_pz", _paclust_pz, "paclust_pz[npaclust]/F");
        _hTrk->Branch("paclust_p", _paclust_p, "paclust_p[npaclust]/F");
        _hTrk->Branch("paclust_edep", _paclust_edep, "paclust_edep[npaclust]/F");
        _hTrk->Branch("paclust_z", _paclust_z , "paclust_z[npaclust]/F");
      }
    }

    if (_paClusterAcceptance <=0) _paClusterAcceptance = 3.;

  }

  void ReadTrkExt::beginRun(art::Run const& run){
    GeomHandle<DetectorSystem> det;
    _origin = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );   // add this to transfer detector coord. to mu2e coord.
    GeomHandle<WorldG4> g4world;
    _mu2eOriginInWorld = g4world->mu2eOriginInWorld(); // add this to transfer mu2e coord. to g4 coord. 
    cerr << "ReadTrkExt_module: Detector coord origin in mu2e coord = " << _origin << endl;
    cerr << "ReadTrkExt_module: Mu2e coord origin in G4 coord = " << _mu2eOriginInWorld << endl;


  }

  void ReadTrkExt::beginSubRun(art::SubRun const& lblock ) {
  }

  void ReadTrkExt::endJob(){
  }


////////// Analysis ///////////

  void ReadTrkExt::analyze(const art::Event& event) {
    _evtid = event.id().event();
    //cerr << "ReadTrkExt: From analyze: " << _evtid << endl;


    art::Handle<KalRepCollection> trksHandle;
    event.getByLabel(_fitterModuleLabel, _trkPatRecInstanceName, trksHandle);
    KalRepCollection const& trks = *trksHandle;

    art::Handle<TrkExtTrajCollection> trkextHandle;
    event.getByLabel(_trkextModuleLabel, trkextHandle);
    TrkExtTrajCollection const & trkexts = *trkextHandle;


    for ( size_t i=0; i< trks.size(); ++i ){
      _trkid = i;
      KalRep const & trk = *(trks.at(i));
      TrkHotList const* hots  = trk.hotList();
      TrkExtTraj const& trkext = trkexts[i];

      _nhits = _nhots = hots->nHit();
      _ntrks = MAXNTRK;
      _nsim = 0;
      

      if (_nhits >(signed int) MAXNHOT || _nhots > (signed int) MAXNHOT ) {
        cerr << "ReadTrkExt_module: Too many hits or hots (" <<_nhits << ", " << _nhots <<  ") per track " << _trkid << " in event " << _evtid <<". Ignored.\n";
        continue;
      }

      readHits (event, hots) ;

      readTrkPatRec (trk);

      vector <TrkExtTrajPoint> pahits;
      vector <TrkExtTrajPoint> sthits;
      doExtDataBooking (trkext, pahits, sthits);
      doPAHitBooking (pahits); 
      doSTHitBooking (sthits); 

      _hTrk->Fill();
      

    }  // end of trks loop

  }


/////////// Read Hit List //////////////

  void ReadTrkExt::readHits (const art::Event& event, TrkHotList const * hits) {
    int i = 0;
    for (TrkHotList::hot_iterator iter = hits->begin() ; iter != hits->end() ; ++iter) {
      const TrkHitOnTrk * hit  = iter.get();
      const mu2e::TrkStrawHit* strawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hit);
      if (strawHit) {
        const HepPoint &point = strawHit->hitTraj()->position(strawHit->hitLen());
        _hotx[i] = point.x();
        _hoty[i] = point.y();
        _hotz[i] = point.z();
        _hott0[i] = strawHit->hitT0()._t0 / CLHEP::ns;
        ++i;
        if (_mcFlag) {
          if (_nsim <=0) _nsim = readSimulation(event, strawHit);
        }
      }
    }
  }

/////////// Read TrkPatRec //////////////

  void ReadTrkExt::readTrkPatRec(KalRep const & krep) {

    double _trkl0 = krep.startValidRange();
    double _trkl1 = krep.endValidRange();
    double step = (_trkl1 - _trkl0)/double(_ntrks);

    for(int i = 0; i < _ntrks; ++i) {
      double fltlen = _trkl0 + i * step;
      _trkl[i] = fltlen;
      HepPoint p = krep.position(fltlen);
      _trkx[i] = p.x();
      _trky[i] = p.y();
      _trkz[i] = p.z();
      HelixParams par = krep.helix(fltlen);
      _trkd0[i] = par.d0();
      _trkz0[i] = par.z0();
      _trkphi[i] = par.phi0();
      _trkomega[i] = par.omega();
      _trktanDip[i] = par.tanDip();
      _trkp[i] = krep.momentum(fltlen).mag();
      _trkpx[i] = krep.momentum(fltlen).x();
      _trkpy[i] = krep.momentum(fltlen).y();
      _trkpz[i] = krep.momentum(fltlen).z();
      Hep3Vector pp (p.x(), p.y(), p.z());

    }
    return;
  }

/////////// Read Simulation //////////////

  int ReadTrkExt::readSimulation(const art::Event& event, const mu2e::TrkStrawHit* trkstrawHit) {
    unsigned int i, j;
    const StrawHit& strawHit = trkstrawHit->strawHit();

    art::Handle<StrawHitCollection> shcHandle;
    event.getByLabel(_makerModuleLabel, shcHandle);
    StrawHitCollection const & shc = *shcHandle;

    int hitid = -1;

    for (i = 0 ; i <shc.size() ; ++i) {
      StrawHit  const& sh(shc.at(i));
      if (sh == strawHit) {
        hitid = int(i);
        break;
      }
    }
    if (hitid <0) return 0;

    art::Handle<PtrStepPointMCVectorCollection> stepsCHandle;
    event.getByLabel(_makerModuleLabel, "StrawHitMCPtr", stepsCHandle);
    PtrStepPointMCVectorCollection const & stepsC = *stepsCHandle;

    if (stepsC.size() <=0) return 0;

    PtrStepPointMCVector const & steps (stepsC.at(hitid));

    if (steps.size() <=0) return 0;

    StepPointMC const & step = *steps.at(0);

    SimParticle const &sim = *step.simParticle();
    
    cet::map_vector_key simid =  sim.id();

    art::Handle<PointTrajectoryCollection> trajHandle;
    event.getByLabel(_g4ModuleLabel, trajHandle);
    const PointTrajectoryCollection & trajC = *trajHandle;

    const PointTrajectory* traj = trajC.getOrNull(simid);
    if (traj == NULL) return 0;

    const vector<Hep3Vector>& pvec = traj->points();
    for (i = 0 ; i < pvec.size() ; ++i) {
      if (i >=MAXNSIM) {
        cerr << "ReadTrkExt_module: Too many sim traj point. ignored\n";
        break;
      }
      Hep3Vector position = pvec[i];
      _simx[i] = position.x();
      _simy[i] = position.y();
      _simz[i] = position.z() + 1800;
    }

    int ret = int(i);

    _simt0 = sim.startGlobalTime();
    _simp0 = sim.startMomentum().vect().mag();


    // read VD
    vector<Hep3Vector> vdp[5];
    vector<Hep3Vector> vdx[5];
    for (i = 0 ; i <5 ; ++i) {
      vdp[i].clear();
      vdx[i].clear();
    }
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,"virtualdetector",hits);
    for ( i=0; i<hits->size(); ++i ){
      const StepPointMC& hit = (*hits)[i];
      if (hit.trackId() != simid) continue;

      switch (hit.volumeId()) {
        case VirtualDetectorId::ST_In :
          j = 0;
          break;
        case VirtualDetectorId::ST_Out :
          j = 1;
          break;
        case VirtualDetectorId::TT_FrontHollow :
        case VirtualDetectorId::TT_FrontPA :
          j = 2;
          break;
        case VirtualDetectorId::TT_MidInner :
        case VirtualDetectorId::TT_Mid :
          j = 3;
          break;
        case VirtualDetectorId::TT_Back :
          j = 4;
          break;
        default :
          j = 9999;
          break;
      }
      if (j >10) continue;
      vdp[j].push_back(hit.momentum());
      vdx[j].push_back(hit.position()-_origin);
    }
    for (i = 0 ; i <5 ; ++i) {
      if (vdp[i].size() >0) {
        Hep3Vector psum(0,0,0);
        Hep3Vector xsum(0,0,0);
        for (unsigned int k = 0 ; k < vdp[i].size() ; ++k) {
          psum += (vdp[i])[k];
          xsum += (vdx[i])[k];
        }

        psum /=  double(vdp[i].size());
        xsum /=  double(vdx[i].size());
        _vdx[i] = xsum.x();
        _vdy[i] = xsum.y();
        _vdz[i] = xsum.z();
        _vdpx[i] = psum.x();
        _vdpy[i] = psum.y();
        _vdpz[i] = psum.z();
        _vdp[i] = psum.mag();
      }
      else {
        _vdx[i] = -99999;
        _vdy[i] = -99999;
        _vdz[i] = -99999;
        _vdpx[i] = -99999;
        _vdpy[i] = -99999;
        _vdpz[i] = -99999;
        _vdp[i] = -99999;
      }
    }

    TrkExtMCHits pareader(_g4ModuleLabel, "protonabsorber", 50., _paClusterAcceptance);
    pareader.readMCHits(event);
    vector<StepPointMCCollection> & pahitcol = pareader.getClusters();
    _npahits = 0;
    _npaclust = pareader.getNClusters();

    for (i = 0 ; i < pahitcol.size() ; ++i) {
      if (i>=100) {
        cerr << "ReadTrkExt_module: Too many (" << _npaclust << ") pa clusters for evt=" << _evtid << ". Skipping remanings" << endl;
        break;
      }
      StepPointMCCollection & mccol = pahitcol[i];
      Hep3Vector psum(0);
      Hep3Vector xsum(0);
      double     edeptot(0);
      for (j = 0 ; j <mccol.size() ; ++j) {
        StepPointMC & mc = mccol[j];
        psum += mc.momentum();
        xsum += mc.position();
        edeptot += mc.totalEDep();
      }
      _paclust_nhits[i] = mccol.size();
      _paclust_px[i] = psum.x() / double(mccol.size());
      _paclust_py[i] = psum.y() / double(mccol.size());
      _paclust_pz[i] = psum.z() / double(mccol.size());
      _paclust_p[i]  = psum.mag() / double(mccol.size());
      _paclust_edep[i] = edeptot; 
      _paclust_z[i] = xsum.z() / double(mccol.size());
    }

    

    return ret;
    
  }


///////// booking ////////////
  void ReadTrkExt::doExtDataBooking (const TrkExtTraj & traj, vector<TrkExtTrajPoint> & pahits, vector<TrkExtTrajPoint> & sthits) {
    unsigned int i;
    double x, y, z, px, py, pz, p, rho;
    pahits.clear();
    sthits.clear();
    for ( i = 0 ; i < traj.size() ; ++i) {
      TrkExtTrajPoint  hit = traj[i];
      _back_x[i] = x = hit.x();
      _back_y[i] = y = hit.y();
      _back_z[i] = z = hit.z();
      _back_rho[i] = rho = hit.position().rho();
      _back_p[i] = p = hit.momentum().mag();
      _back_px[i] = px = hit.px();
      _back_py[i] = py = hit.py();
      _back_pz[i] = pz = hit.pz();
      _back_ex[i] = hit.ex(); 
      _back_ey[i] = hit.ey(); 
      _back_ez[i] = hit.ez(); 
      _back_epx[i] = hit.epx();
      _back_epy[i] = hit.epy();
      _back_epz[i] = hit.epz();
      _back_ep[i] = safeSqrt( hit.covpxpx()*px*px + hit.covpypy()*py*py + hit.covpzpz()*pz*pz + 2.*px*py*hit.covpxpy() + 2.*py*pz*hit.covpypz() + 2.*pz*px*hit.covpxpz() ) / p;
      _back_er[i] = safeSqrt( hit.covxx()*x*x + hit.covyy()*y*y + 2.*x*y*hit.covxy() ) / rho;
      _back_vid[i] = hit.volumeId();

      if (_back_vid[i] == TrkExtDetectorList::ProtonAbsorber) pahits.push_back(hit);
      else if (i>0 && _back_vid[i-1] == TrkExtDetectorList::ProtonAbsorber) pahits.push_back(hit);
      else {;}

      if (_back_vid[i] == TrkExtDetectorList::StoppingTarget) sthits.push_back(hit);
      else if (i>0 && _back_vid[i-1] == TrkExtDetectorList::StoppingTarget) sthits.push_back(hit);
      else {;}
    }
    if (traj.size() <=0) {
      _back_s = 0;
      _back_t = 0;
      _exitcode = traj.exitCode();
      _nback = 0;
    }
    else {
      _back_s = traj.flightLength();
      _back_t = traj.flightTime();
      _exitcode = traj.exitCode();
      _nback = int(i);
    }
    return; 
  }




///////// PA/ST Hit anal ////////////

  void ReadTrkExt::doPAHitBooking ( vector<TrkExtTrajPoint> & _pahitdata) {
    unsigned int i, j, n;
    n = _pahitdata.size();
    if (n <=0) {
      _nbackpaclust = 0;
      return ;
    }

    vector < vector<TrkExtTrajPoint> > clusters;
    clusters.clear();

    vector<TrkExtTrajPoint> cluster;
    cluster.push_back(_pahitdata.front());
    clusters.push_back(cluster);
    bool flag;
    for (i = 1 ; i < _pahitdata.size() ; ++i) {
      TrkExtTrajPoint &hit = _pahitdata[i];
      flag = false;
      for (j = 0 ; j < clusters.size() ; ++j) {
        vector<TrkExtTrajPoint> & cluster = clusters[j];
        if (hit.trajPointId() == cluster.back().trajPointId() +1) {
          cluster.push_back(hit);
          flag = true;
          break;
        }
      }
      if (!flag) {
        vector<TrkExtTrajPoint> cluster;
        cluster.push_back(hit);
        clusters.push_back(cluster);
      }
    }

    double dpsum(0);
    for (i = 0 ; i < clusters.size() ; ++i) {
      vector<TrkExtTrajPoint> & cluster = clusters[i];
      double zsum(0);
      for (j = 0 ; j < cluster.size() ; ++j) {
        zsum += (cluster[j]).z();
      }
      _back_paclust_z[i] = zsum / double(cluster.size());
      _back_paclust_dp[i] = cluster.back().momentum().mag() - cluster.front().momentum().mag();
      dpsum += _back_paclust_dp[i];
    }
    _nbackpaclust = clusters.size();
    _back_paclust_dptot = dpsum;

    return; 
  }



  void ReadTrkExt::doSTHitBooking ( vector<TrkExtTrajPoint> & _sthitdata) {
    unsigned int i, j, n;
    n = _sthitdata.size();
    if (n <=0) {
      _nbackstclust = 0;
      return;
    }

    vector < vector<TrkExtTrajPoint> > clusters;
    clusters.clear();

    vector<TrkExtTrajPoint> cluster;
    cluster.push_back(_sthitdata.front());
    clusters.push_back(cluster);
    bool flag;
    for (i = 1 ; i < _sthitdata.size() ; ++i) {
      TrkExtTrajPoint &hit = _sthitdata[i];
      flag = false;
      for (j = 0 ; j < clusters.size() ; ++j) {
        vector<TrkExtTrajPoint> & cluster = clusters[j];
        if (hit.trajPointId() == cluster.back().trajPointId() +1) {
          cluster.push_back(hit);
          flag = true;
          break;
        }
      }
      if (!flag) {
        vector<TrkExtTrajPoint> cluster;
        cluster.push_back(hit);
        clusters.push_back(cluster);
      }
    }

    double dpsum(0);
    for (i = 0 ; i < clusters.size() ; ++i) {
      vector<TrkExtTrajPoint> & cluster = clusters[i];
      double zsum(0);
      for (j = 0 ; j < cluster.size() ; ++j) {
        zsum += (cluster[j]).z();
      }
      _back_stclust_z[i] = zsum / double(cluster.size());
      _back_stclust_dp[i] = cluster.back().momentum().mag() - cluster.front().momentum().mag();
      dpsum += _back_stclust_dp[i];
    }
    _back_stclust_dptot = dpsum;
    _nbackstclust = clusters.size();

    return ; 
  }



} // end namespace mu2e

using mu2e::ReadTrkExt;
DEFINE_ART_MODULE(ReadTrkExt);
