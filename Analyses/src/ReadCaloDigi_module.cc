//
// An EDAnalyzer module that reads back the hits created by the Calorimeter Digitization chain
//
// Original author

// ROOT includes
#include "TH1F.h"
#include "TF1.h"
#include "TSpline.h"
#include "TFile.h"
// #include "CalPatRec/inc/THackData.hh"
#include "TROOT.h"
#include "TFolder.h"
#include "TTree.h"
#include "TH2F.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

//tracker includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkRep.hh"

#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

//#include "TrkReco/inc/TrkStrawHit.hh"
#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {
  struct DAQ_t {
    TH1F*       _hNCaloDigi       [10];
    TH1F*       _hNSamplesPerDigi [10];
    TH1F*       _hNSamplesPerEvent[10];

    TH2F*       _hNSamplesVsRadius[10];
    TH2F*       _hNHitsVsRadius   [10];
   };

  class ReadCaloDigi : public art::EDAnalyzer {

  public:

    explicit ReadCaloDigi(fhicl::ParameterSet const& pset);
    virtual ~ReadCaloDigi() { }

    virtual void beginRun(art::Run const& ) override;

    virtual void beginJob() override;
    virtual void endJob() override ;

    // This is called for each event.
    virtual void analyze(const art::Event& e) override;

  private:

    typedef std::vector< art::Handle<mu2e::StepPointMCCollection> > StepMCHandleVector;
    typedef art::Ptr<mu2e::CaloHit>                          CaloHitPtr;

    void             initVHits     ();

    int                        _diagLevel;

    std::string                _caloDigisModuleLabel;
    std::string                _caloCrystalModuleLabel;
    std::string                _stepPoints;
    std::string                _rostepPoints;
    std::string                _caloClusterModuleLabel;
    std::string                _vdStepPoints;
    std::string                _trackModuleLabel;

    std::string                _instanceName;
    TrkParticle                _tpart;
    TrkFitDirection            _fdir;

    double                     _mbtime;
    double                     _mbbuffer;
    double                     _blindTime;

    int                        _fillWaveforms;
    double                     _psdThreshold;
    double                     _caloRmin;

    TTree*                     _Ntup;
    DAQ_t                      _histDisk0;
    DAQ_t                      _histDisk1;
    TH2F*                      _hNSamplesVsCrysId[10];

    int                        _nProcess;

    int                        _evt,_run,_caloCrystals,_caloDisk0Crystals,_caloDisk1Crystals,_nHits,_nCluster, _nCryRO;

    int                        _nRecoDigi, _recoDigiSamples[16384];
    float                      _recoDigiAmp[16384], _recoDigiEnergy[16384], _recoDigiTime[16384];

    int                        _recoDigiId[16384];
    float                      _recoDigiT0[16384];
    float                      _recoDigiX[16384], _recoDigiY[16384], _recoDigiZ[16384];
    float                      _recoDigiPulse[10000][350];

    int                        _cryId[16384],_crySectionId[16384];

    int                        _cluNcrys[204828];

    float                      _caloVolume, _crystalVolume;

    float                      _cryEtot,_cryTime[16384],_cryEdep[16384], _cryAmp[16384], _cryDose[16384];
    float                      _cryMCTime    [16384];
    float                      _cryMCEdep    [16384];
    int                        _cryNParticles[16384];
    float                      _cryPsd       [16384];
    float                      _cryIsConv    [16384];

    float                      _cryPosX[16384],_cryPosY[16384],_cryPosZ[16384], _cryPosR[16384];

    float                      _cluEnergy[2048], _cluCrysE[2048], _cluMeanTime[2048],_cluTime[2048], _cluMCTime[2048], _cluMCMeanTime[2048];
    float                      _cluCogX[2048],_cluCogY[2048],_cluCogZ[2048], _cluCogR[2048];

    int                        _cluConv[2048];

    Int_t                     _vNHits;
    Float_t                   _vP[2048], _vPx[2048], _vPy[2048], _vPz[2048], _vPt[2048];
    Float_t                   _vE[2048], _vEKin[2048], _vM[2048];
    Float_t                   _vT[2048];
    Float_t                   _vX[2048], _vY[2048], _vZ[2048], _vR[2048];
    Float_t                   _vCosth[2048], _vRadius[2048];
    Int_t                     _vId[2048];
    Int_t                     _vPdgId[2048];

    Int_t                     _nTrkGood;

    //some histograms for DAQ analyses purposes

    const Calorimeter*        _calorimeter; // cached pointer to the calorimeter geometry

  };

  void     ReadCaloDigi::initVHits(){

    _vNHits  = 0;

    // _vP      = 0.;
    // _vPx     = 0.;
    // _vPy     = 0.;
    // _vPz     = 0.;
    // _vPt     = 0.;
    // _vPdgId  = 0.;
    // _vM      = 0.;
    // _vE      = 0.;
    // _vEKin   = 0.;
    // _vT      = 0.;
    // _vX      = 0.;
    // _vY      = 0.;
    // _vZ      = 0.;
    // _vR      = 0.;

    // _vCosth  = 0.;
    // _vRadius = 0.;
    // _vId     = 0.;
  }

  ReadCaloDigi::ReadCaloDigi(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel                     (pset.get<int>   ("diagLevel")),
    _caloDigisModuleLabel          (pset.get<string>("caloDigisModuleLabel")),
    _caloCrystalModuleLabel        (pset.get<string>("caloCrystalModuleLabel")),
    _stepPoints                    (pset.get<string>("calorimeterStepPoints")),
    _rostepPoints                  (pset.get<string>("calorimeterROStepPoints")),
    _caloClusterModuleLabel        (pset.get<string>("caloClusterModuleLabel")),
    _vdStepPoints                  (pset.get<string>("vdStepPoints")),
    _trackModuleLabel              (pset.get<string>("trackModuleLabel")),
    _tpart     ((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _mbbuffer                      (pset.get<double>             ("TimeFoldingBuffer")),  // ns
    _blindTime                     (pset.get<double>             ("blindTime" )),         // ns
    _fillWaveforms                 (pset.get<int>                ("fillWaveforms" )),
    _psdThreshold                  (pset.get<double>("psdThreshold")),
    _caloRmin                      (pset.get<double>("caloRmin" )),
    _Ntup(0),_nProcess(0)
  {
    _instanceName = _fdir.name() + _tpart.name();
  }

  void ReadCaloDigi::beginRun(art::Run const& ){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
  }
  void ReadCaloDigi::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");

    for (int i=0; i<10; ++i){
      _histDisk0._hNCaloDigi        [i] = tfdir.make<TH1F>(Form("hNHitsDisk0%i",i)      ,
                                                          Form("Disk0: N caloDigis @ %i MeV threshold",i+1), 1e4, 0, 1e4);
      _histDisk0._hNSamplesPerDigi  [i] = tfdir.make<TH1F>(Form("hNSampleHitsDisk0%i",i),
                                                          Form("Disk0: N of words per caloDigi @ %i MeV threshold",i+1),
                                                           1e4, 0, 1e4);
      _histDisk0._hNSamplesPerEvent [i] = tfdir.make<TH1F>(Form("hNSampleDisk0%i",i),
                                                          Form("Disk0: N of words per event @ %i MeV threshold",i+1),
                                                           5e4, 0, 5e4);
      _histDisk1._hNCaloDigi        [i] = tfdir.make<TH1F>(Form("hNHitsDisk1%i",i)      ,
                                                          Form("Disk1: N caloDigis @ %i MeV threshold",i+1), 1e4, 0, 1e4);
      _histDisk1._hNSamplesPerDigi  [i] = tfdir.make<TH1F>(Form("hNSampleHitsDisk1%i",i),
                                                          Form("Disk1: N of words per caloDigi @ %i MeV threshold",i+1),
                                                           1e4, 0, 1e4);
      _histDisk1._hNSamplesPerEvent [i] = tfdir.make<TH1F>(Form("hNSampleDisk1%i",i),
                                                          Form("Disk1: N of words per event @ %i MeV threshold",i+1),
                                                           5e4, 0, 5e4);
      _hNSamplesVsCrysId            [i] = tfdir.make<TH2F>(Form("hNSamplesVsCrysId%i",i),
                                                          Form("N of words per event vs crystal Id @ %i MeV threshold",i+1),
                                                          2000, 0, 2e3, 600, 0, 600);

      _histDisk0._hNSamplesVsRadius [i] = tfdir.make<TH2F>(Form("hNSamplesVsRadiusDisk0%i",i),
                                                          "N of words per event vs radius on disk 0",
                                                          80, 300, 700, 300, 0, 300);
      _histDisk0._hNHitsVsRadius    [i] = tfdir.make<TH2F>(Form("hNHitsVsRadiusDisk0%i",i),
                                                          "N of hits per event vs radius on disk 0",
                                                          80, 300, 700, 2000, 0, 2000);

      _histDisk1._hNSamplesVsRadius [i] = tfdir.make<TH2F>(Form("hNSamplesVsRadiusDisk1%i", i),
                                                          "N of words per event vs radius on disk 0",
                                                          80, 300, 700, 300, 0, 300);
      _histDisk1._hNHitsVsRadius    [i] = tfdir.make<TH2F>(Form("hNHitsVsRadiusDisk1%i", i),
                                                          "N of hits per event vs radius on disk 0",
                                                          80, 300, 700, 2000, 0, 2000);
    }

    _Ntup  = tfs->make<TTree>("Calo", "Calo");

    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");
    _Ntup->Branch("caloCrystals", &_caloCrystals ,"caloCrystals/I");
    _Ntup->Branch("caloDisk0Crystals", &_caloDisk0Crystals ,"caloDisk0Crystals/I");
    _Ntup->Branch("caloDisk1Crystals", &_caloDisk1Crystals ,"caloDisk1Crystals/I");

    _Ntup->Branch("nRecoDigi",     &_nRecoDigi ,       "nRecoDigi/I");
    _Ntup->Branch("recoDigiEnergy",        &_recoDigiEnergy ,       "recoDigiEnergy[nRecoDigi]/F");
    _Ntup->Branch("recoDigiAmp", &_recoDigiAmp, "recoDigiAmp[nRecoDigi]/F");
    _Ntup->Branch("recoDigiTime"   , &_recoDigiTime, "recoDigiTime[nRecoDigi]/F");
    _Ntup->Branch("recoDigiSamples", &_recoDigiSamples, "recoDigiSamples[nRecoDigi]/I");
    _Ntup->Branch("recoDigiId"     , &_recoDigiId     , "recoDigiId[nRecoDigi]/I");
    _Ntup->Branch("recoDigiT0"     , &_recoDigiT0     , "recoDigiT0[nRecoDigi]/F");
    _Ntup->Branch("recoDigiX"     , &_recoDigiX     , "recoDigiX[nRecoDigi]/F");
    _Ntup->Branch("recoDigiY"     , &_recoDigiY     , "recoDigiY[nRecoDigi]/F");
    _Ntup->Branch("recoDigiZ"     , &_recoDigiZ     , "recoDigiZ[nRecoDigi]/F");
    _Ntup->Branch("recoDigiPulse"  , &_recoDigiPulse  , "recoDigiPulse[nRecoDigi][350]/F");

    _Ntup->Branch("cryEtot",      &_cryEtot ,    "cryEtot/F");

    _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
    _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
    _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
    _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
    _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
    _Ntup->Branch("cryPosR",      &_cryPosR ,     "cryPosR[nCry]/F");
    _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
    _Ntup->Branch("cryAmp",       &_cryAmp ,      "cryAmp[nCry]/F");
    _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
    _Ntup->Branch("cryMCTime",    &_cryMCTime    ,"cryMCTime[nCry]/F");
    _Ntup->Branch("cryMCEdep",    &_cryMCEdep    ,"cryMCEdep[nCry]/F");
    _Ntup->Branch("cryNParticles",&_cryNParticles,"cryNParticles[nCry]/I");
    _Ntup->Branch("cryPsd",       &_cryPsd       ,"cryPsd[nCry]/F");
    _Ntup->Branch("cryIsConv",    &_cryIsConv    ,"cryIsConv[nCry]/I");
    _Ntup->Branch("nCluster",     &_nCluster ,    "nCluster/I");
    _Ntup->Branch("cluEnergy",    &_cluEnergy ,   "cluEnergy[nCluster]/F");
    _Ntup->Branch("cluCrysE",     &_cluCrysE ,    "cluCrysE[nCluster]/F");
    _Ntup->Branch("cluTime",      &_cluTime ,     "cluTime[nCluster]/F");
    _Ntup->Branch("cluMeanTime",      &_cluMeanTime ,     "cluMeanTime[nCluster]/F");
    _Ntup->Branch("cluMCTime",    &_cluMCTime ,   "cluMCTime[nCluster]/F");
    _Ntup->Branch("cluMCMeanTime",    &_cluMCMeanTime ,   "cluMCMeanTime[nCluster]/F");
    _Ntup->Branch("cluCogX",      &_cluCogX ,     "cluCogX[nCluster]/F");
    _Ntup->Branch("cluCogY",      &_cluCogY ,     "cluCogY[nCluster]/F");
    _Ntup->Branch("cluCogZ",      &_cluCogZ ,     "cluCogZ[nCluster]/F");
    _Ntup->Branch("cluCogR",      &_cluCogR ,     "cluCogR[nCluster]/F");
    _Ntup->Branch("cluNcrys",     &_cluNcrys ,    "cluNCrys[nCluster]/I");
    _Ntup->Branch("cluConv",      &_cluConv ,     "cluConv[nCluster]/I");

    _Ntup->Branch("vNHits",  &_vNHits ,  "vNHits/I");
    _Ntup->Branch("vId"   ,  &_vId    ,  "vId[vNHits]/I");
    _Ntup->Branch("vPdgId",  &_vPdgId ,  "vPdgId[vNHits]/I");
    _Ntup->Branch("vP" ,     &_vP     ,  "vP[vNHits]/F");
    _Ntup->Branch("vPx",     &_vPx    ,  "vPx[vNHits]/F");
    _Ntup->Branch("vPy",     &_vPy    ,  "vPy[vNHits]/F");
    _Ntup->Branch("vPz",     &_vPz    ,  "vPz[vNHits]/F");
    _Ntup->Branch("vE" ,     &_vE     ,  "vE[vNHits]/F");
    _Ntup->Branch("vEKin",   &_vEKin  ,  "vEKin[vNHits]/F");
    _Ntup->Branch("vM",      &_vM     ,  "vM[vNHits]/F");
    _Ntup->Branch("vT",      &_vT     ,  "vT[vNHits]/F");
    _Ntup->Branch("vX",      &_vX     ,  "vX[vNHits]/F");
    _Ntup->Branch("vY",      &_vY     ,  "vY[vNHits]/F");
    _Ntup->Branch("vZ",      &_vZ     ,  "vZ[vNHits]/F");
    _Ntup->Branch("vR",      &_vR     ,  "vR[vNHits]/F");
    _Ntup->Branch("vCosth",  &_vCosth ,  "vCosth[vNHits]/F");
    _Ntup->Branch("vRadius", &_vRadius,  "vRadius[vNHits]/F");

    _Ntup->Branch("nTrkGood",&_nTrkGood ,  "nTrkGood/I");

  }

  void ReadCaloDigi::endJob(){
  }

  void ReadCaloDigi::analyze(const art::Event& event) {
    // if (_nProcess == 0){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    // }

    ++_nProcess;

    //load the timeoffset
    _mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();

    //data about hits in the calorimeter crystals
    art::Handle<CaloHitMCCollection> caloDigiMCHandle;
    event.getByLabel(_caloDigisModuleLabel, caloDigiMCHandle);

    const     CaloHitMCCollection*       caloDigiMCCol(0);
    int       nCaloHitMC(0);
    if (caloDigiMCHandle.isValid()){
       caloDigiMCCol = caloDigiMCHandle.product();
       nCaloHitMC   = caloDigiMCCol->size();
    }

    art::Handle<CaloRecoDigiCollection>   recoCaloDigiHandle;
    event.getByLabel(_caloDigisModuleLabel, recoCaloDigiHandle);

    const     CaloRecoDigiCollection*       recoCaloDigiCol(0);
    int       nCaloRecoDigi(0);
    if (recoCaloDigiHandle.isValid()){
       recoCaloDigiCol = recoCaloDigiHandle.product();
       nCaloRecoDigi   = recoCaloDigiCol->size();
    }

    int         nCaloCrystals(0);
    const     CaloHitCollection*  CaloHits(0);
    art::Handle<CaloHitCollection> CaloHitsHandle;
    event.getByLabel(_caloCrystalModuleLabel, CaloHitsHandle);
    if (!CaloHitsHandle.isValid()){
      //      printf("[ReadCaloDigi::analyze] no CaloHitCollection: BAILS OUT \n");
      nCaloCrystals = 0;
    }else {
     CaloHits = CaloHitsHandle.product();
     nCaloCrystals   = CaloHits->size();
    }

    //data about clusters
    art::Handle<CaloClusterCollection> caloClustersHandle;
    const     CaloClusterCollection*  caloClusters(0);
    int       nCaloClusters(0);
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    if (caloClustersHandle.isValid()){
      caloClusters  = caloClustersHandle.product();
      nCaloClusters = caloClusters->size();
    }

    //Handle to VD steps
    art::ProductInstanceNameSelector selector_vdhits("virtualdetector");
    art::Handle<StepPointMCCollection> *vdStepsHandle;
    const StepPointMCCollection *vdHits;
    StepMCHandleVector vdStepsHandleVec = event.getMany<StepPointMCCollection>(selector_vdhits);

    //Handle to tracks collection
    art::Handle<mu2e::KalRepCollection> dem_handle;
    // art::Selector  s_dem(art::ProcessNameSelector("")                &&
    //                   art::ModuleLabelSelector(_trackModuleLabel)    );
    // art::Selector  s_dem   (art::ProductInstanceNameSelector("") &&
    //                      art::ProcessNameSelector("")         &&
    //                      art::ModuleLabelSelector(_trackModuleLabel)            );
    //    fEvent->get(selector,krepsHandle);
    event.getByLabel(_trackModuleLabel, _instanceName, dem_handle);
    //    event.get(s_dem,dem_handle);

    int      nTraks(0);
    const mu2e::KalRepCollection* list_of_ele_tracks = nullptr;

    if (dem_handle.isValid()) {
      list_of_ele_tracks = dem_handle.product();
      nTraks             = list_of_ele_tracks->size();
    }

    GlobalConstantsHandle<ParticleDataList> pdt;

    _evt          = event.id().event();
    _run          = event.run();
    _caloCrystals = _calorimeter->nCrystals();
    _caloDisk0Crystals = _calorimeter->disk(0).nCrystals();
    _caloDisk1Crystals = _calorimeter->disk(1).nCrystals();
    _nTrkGood          = 0;

    //--------------------------------------------------------------------------------
    if (nTraks>0){
      double  fMinFitCons      = 2.e-3;
      double  fMinNActive      = 25;
      double  fMaxT0Err        = 0.9;           // in ns
      double  fMaxFitMomErr    = 0.25;                  // in MeV
      double  fMinTanDip       = tan(M_PI/6.);  // 0.5773
      double  fMaxTanDip       = 1.0;
      double  fMinD1           = -80.;          // in mm
      double  fMaxD1           = 105.;
      double  fTrkMomCut       = 100.;            //MeV/c

      KalRep *trk;

      trk = (KalRep*) &list_of_ele_tracks->at(0);
      CLHEP::Hep3Vector mom = trk->momentum(0);

      BbrVectorErr      momerr = trk->momentumErr(0);

      CLHEP::Hep3Vector momdir = trk->momentum(0).unit();
      CLHEP::HepVector  momvec(3);
      for (int i=0; i<3; i++) momvec[i] = momdir[i];

      double fitmom_err = sqrt(momerr.covMatrix().similarity(momvec));

      if ( (trk->chisqConsistency().consistency() > fMinFitCons  ) &&
           (trk->nActive()                        > fMinNActive  ) &&
           (trk->t0().t0Err()                     < fMaxT0Err    ) &&
           (fitmom_err                             < fMaxFitMomErr) &&
           (trk->helix(0).tanDip()                > fMinTanDip   ) &&
           (trk->helix(0).tanDip()                < fMaxTanDip   ) &&
           (trk->helix(0).d0()                    < fMaxD1       ) &&
           (trk->helix(0).d0()                    > fMinD1       ) &&
           (mom.mag()                             > fTrkMomCut   )){

        _nTrkGood = 1;
      }
    }
    //--------------------------------------------------------------------------------

    int           vdVecSize = vdStepsHandleVec.size();
    initVHits();

    //    double        timeLastHit(1e10);
    //    double        timeLastHit(0);

    for (int j=0; j< vdVecSize; ++j){
      vdStepsHandle = & vdStepsHandleVec[j];
      if (vdStepsHandle->isValid())  {
        vdHits        = vdStepsHandle->operator->();

        for (size_t i=0; i<vdHits->size(); ++i) {
          StepPointMC hit = vdHits->at(i);

          if (hit.simParticle()->fromGenerator()) {
            int id = hit.volumeId();

            if (id == VirtualDetectorId::EMC_Disk_0_SurfIn  ||
                id == VirtualDetectorId::EMC_Disk_1_SurfIn  ||
                id == VirtualDetectorId::EMC_Disk_0_EdgeIn  ||
                id == VirtualDetectorId::EMC_Disk_1_EdgeIn    ) {

              art::Ptr<SimParticle> const& simptr = hit.simParticle();
              //2016-01-10 G. PEzzullo temporary comment for using
              // a custom made gen particle

              SimParticle const& sim  = *simptr;

              if ( sim.fromGenerator() ){

                GenParticle* gen = (GenParticle*) &(*sim.genParticle());
                if ( gen->generatorId().isConversion() ){
                  continue;
                }
              }
              if ( !sim.fromGenerator() )                 continue;

              double     hitTime = fmod(hit.time(), _mbtime);
              if (hitTime < _mbbuffer) {
                if (hitTime+_mbtime > _blindTime) {
                  hitTime = hitTime + _mbtime;
                }
              }
              else {
                if (hitTime > (_mbtime - _mbbuffer)) {
                  if (hitTime - _mbtime > _blindTime) {
                    hitTime =   hitTime - _mbtime;
                  }
                }
              }

              //              if (hitTime > timeLastHit)                  continue;
              //              if (hitTime < timeLastHit)                  continue;

              //              timeLastHit = hitTime;

              _vT    [_vNHits]  = hitTime;

              _vP    [_vNHits]  = hit.momentum().mag();
              _vPx   [_vNHits]  = hit.momentum().x();
              _vPy   [_vNHits]  = hit.momentum().y();
              _vPz   [_vNHits]  = hit.momentum().z();
              _vPt   [_vNHits]  = std::sqrt( std::pow(_vPx[_vNHits],2.)+std::pow(_vPy[_vNHits],2.) );
              _vPdgId[_vNHits]  = hit.simParticle()->pdgId();
              _vM    [_vNHits]  = pdt->particle(_vPdgId[_vNHits]).mass();
              _vE    [_vNHits]  = sqrt(_vP[_vNHits]*_vP[_vNHits] + _vM[_vNHits]*_vM[_vNHits]);
              _vEKin [_vNHits]  = _vE[_vNHits] - _vM[_vNHits];

              _vX    [_vNHits]  = hit.position().x()+3904.;
              _vY    [_vNHits]  = hit.position().y();
              _vZ    [_vNHits]  = hit.position().z();
              _vR    [_vNHits]  = sqrt(_vX[_vNHits]*_vX[_vNHits] + _vY[_vNHits]*_vY[_vNHits]);

              _vCosth [_vNHits] = _vPz[_vNHits]/_vP[_vNHits];
              _vRadius[_vNHits] = std::sqrt(_vX[_vNHits]*_vX[_vNHits]+_vY[_vNHits]*_vY[_vNHits]);
              _vId    [_vNHits] = id;

              ++_vNHits;
            }
          }
        }
      }
    }

    _nHits = 0;
    _cryEtot = 0.0;

    //some helper variables
    const CaloCluster                       *cluster;
    const std::vector<CaloHitPtr>    *crystals;
    const CaloHit                    *crystalHit;
    const CaloRecoDigi                      *recoDigi;
    //    const CaloDigi                          *caloDigi;
    const CaloHitMC                        *caloDigiMC(0);
    const SimParticle                       *sim;

    //fill DAQ histograms
    double     amplitude(0);
    int        roId(0), crystalID(0), diskId(0);

    int        nThresholds(10), nWords(0);
    int        disk0NCaloDigi       [10] = {0};
    int        disk1NCaloDigi       [10] = {0};
    int        disk0NSamplesPerEvent[10] = {0};
    int        disk1NSamplesPerEvent[10] = {0};
    double     thresholds           [10] = {1., 2., 3, 4., 5., 6, 7, 8, 9, 10};
    //    double     ADC2mV                    = 2000./pow(2.,12);
    double     radius(0);

    //2016-01-27 G. Pezzullo conversion from thresholds given in MeV to mV
    double     p1(21.65);//for the pulse wfInput =2
    //    double     p1(14.76);//for the pulse wfInput =3

    for(int i=0; i<nThresholds; ++i){
      thresholds [i] = thresholds[i]*p1;// + p0;
    }

    int        ncrystals(_caloCrystals);
    vector<int>        nWordsCrystals(ncrystals,0);

    std::vector<int>   pulse;
    CLHEP::Hep3Vector  crystalPos(0);
    int        nROs = _caloCrystals*2;

    double     nWorsRO   [10][4000];
    double     radiusRO  [10][4000];
    double     diskIdRO  [10][4000];
    int        nHitsRO   [10][4000];
    for (int i=0;i<nThresholds; ++i){
      for (int j=0; j<nROs; ++j){
        nWorsRO   [i][j] = 0;
        radiusRO  [i][j] = 0;
        diskIdRO  [i][j] = 0;
        nHitsRO   [i][j] = 0;
      }
    }

    _nRecoDigi = nCaloRecoDigi;

    for (int i=0; i< nCaloRecoDigi; ++i){
      recoDigi   = &recoCaloDigiCol->at(i);
      //amplitude  = recoDigi->amplitude()*ADC2mV;
      roId       = recoDigi->SiPMID();
      crystalID  = CaloSiPMId(roId).crystal().id();
      diskId     = _calorimeter->crystal(crystalID).diskID();

      crystalPos = _calorimeter->geomUtil().mu2eToDiskFF(diskId,_calorimeter->crystal(crystalID).position());
      radius     = sqrt(crystalPos.x()*crystalPos.x() + crystalPos.y()*crystalPos.y());

      if (radius < _caloRmin)               continue;

      _recoDigiEnergy[i] = recoDigi->energyDep();

      const CaloDigi&        caloDigi = *recoDigi->caloDigiPtr();

      pulse      = caloDigi.waveform();
      nWords     = pulse.size();
      //get the amplitude
      for (int j=0; j<nWords; ++j){
        double content = pulse.at(j);
        if (content > amplitude) amplitude = content;
      }
      _recoDigiAmp     [i] = amplitude;

      _recoDigiSamples [i] = nWords;
      _recoDigiId      [i] = roId;
      _recoDigiT0      [i] = caloDigi.t0();
      _recoDigiX       [i] = crystalPos.x();
      _recoDigiY       [i] = crystalPos.y();
      _recoDigiZ       [i] = diskId;
      _recoDigiTime    [i] = recoDigi->time();

      if (_fillWaveforms == 1){
        for (int j=0; j<nWords; ++j){
          _recoDigiPulse[i][j] = pulse.at(j);//*ADC2mV;
        }
      }

      nWordsCrystals[crystalID] += nWords;

      for (int k=0; k<nThresholds; ++k){
        if (amplitude > thresholds[k]){
          nWorsRO   [k][roId] += nWords;
          radiusRO  [k][roId] = radius;
          diskIdRO  [k][roId] = diskId;
          nHitsRO   [k][roId] += 1;

          if (diskId ==0){
            disk0NCaloDigi[k]        = disk0NCaloDigi[k] + 1;
            disk0NSamplesPerEvent[k] = disk0NSamplesPerEvent[k] + nWords ;
            _histDisk0._hNSamplesPerDigi [k]->Fill(nWords);
          }else if (diskId == 1){
            disk1NCaloDigi[k]        = disk1NCaloDigi[k] + 1;
            disk1NSamplesPerEvent[k] = disk1NSamplesPerEvent[k] + nWords ;
            _histDisk1._hNSamplesPerDigi [k]->Fill(nWords);
          }
        }
      }//end loop on the thresholds

    }

    for (int k=0; k<nThresholds; ++k){
      _histDisk0._hNCaloDigi       [k]->Fill(disk0NCaloDigi       [k]);
      _histDisk0._hNSamplesPerEvent[k]->Fill(disk0NSamplesPerEvent[k]);

      _histDisk1._hNCaloDigi       [k]->Fill(disk1NCaloDigi       [k]);
      _histDisk1._hNSamplesPerEvent[k]->Fill(disk1NSamplesPerEvent[k]);

      double   content(0);

      for (int i=0; i<nROs; ++i){
        content     =  nWorsRO [k][i];
        radius      =  radiusRO[k][i];
        if (content > 0){
          if (    diskIdRO[k][i] == 0){
            _histDisk0._hNSamplesVsRadius[k]->Fill(radius, content);
            _histDisk0._hNHitsVsRadius   [k]->Fill(radius, nHitsRO[k][i]);
          }else {
            _histDisk1._hNSamplesVsRadius[k]->Fill(radius, content);
            _histDisk1._hNHitsVsRadius   [k]->Fill(radius, nHitsRO[k][i]);
          }
        }
      }

      for (int j=0; j<ncrystals; ++j){
        if (nWordsCrystals[j] > 0){
          _hNSamplesVsCrysId[k]->Fill(j, nWordsCrystals[j]);
        }
      }

    }

    for (int ic=0; ic<nCaloCrystals;++ic) {

      CaloHit const& hit    = CaloHits->at(ic);
      CLHEP::Hep3Vector crystalPos = _calorimeter->geomUtil().mu2eToDiskFF(diskId,_calorimeter->crystal(crystalID).position());

      _cryEtot             += hit.energyDep();
      _cryTime[_nHits]      = hit.time();
      _cryEdep[_nHits]      = hit.energyDep();

      _cryPosX[_nHits]      = crystalPos.x();
      _cryPosY[_nHits]      = crystalPos.y();
      _cryPosZ[_nHits]      = crystalPos.z();
      _cryPosR[_nHits]      = sqrt( _cryPosX[_nHits]*_cryPosX[_nHits] + crystalPos.y()*crystalPos.y() );
      _cryId[_nHits]        = hit.crystalID();
      _crySectionId[_nHits] = _calorimeter->crystal(hit.crystalID()).diskID();

      int    indexMC, nParticles(0);
      int    isConversion(0);

      recoDigi         = hit.recoCaloDigis().at(0).operator ->();
      //_cryAmp [_nHits]      = recoDigi->amplitude();

      indexMC          = 0;//caloDigi.index();

      if (nCaloHitMC > 0) {
        caloDigiMC       = &caloDigiMCCol->at(indexMC);

        for (unsigned k=0; k<caloDigiMC->nParticles(); ++k){
          sim =   caloDigiMC->energyDeposit(k).sim().operator ->();
          //    if ( sim->fromGenerator() ){
          const CLHEP::Hep3Vector genPos = sim->startPosition();

          // if (!_calorimeter->isInsideCalorimeter(genPos)){
          //   ++nParticles;

          //   int        pdgId       = sim->pdgId();
          //   double     ceEnergy    = 104.9;
          //   double     startEnergy = sim->startMomentum().e();
          //   if ( (pdgId == 11) && (startEnergy>ceEnergy))
          //     {
          //       isConversion = 1;
          //     }
          // }
        }//end loop on the particles inside the crystalHit

        _cryMCTime    [_nHits] = caloDigiMC->time();
        _cryMCEdep    [_nHits] = caloDigiMC->totalEnergyDep();
      } else {
        _cryMCTime    [_nHits] = 0;
        _cryMCEdep    [_nHits] = 0;
      }
      _cryNParticles[_nHits] = nParticles;
      //_cryPsd       [_nHits] = recoDigi->psd();
      _cryIsConv    [_nHits] = isConversion;

      ++_nHits;
    }

    _nCluster = nCaloClusters;//caloClusters->size();

    for (int i=0; i<_nCluster; ++i){
      cluster  = &caloClusters->at(i);

      crystals = &cluster->caloHitsPtrVector();

      int    nCrystals = crystals->size();
      int    indexMC;
      int    isConversion(0);

      double   energyMax(0), eMeanTot(0), clusterTime(0), clusterMCMeanTime(0), clusterMeanTime(0), clusterMCTime(0), eDep, psd, crystalTime(0);

      for (int j=0; j<nCrystals; ++j){
              crystalHit         = crystals->at(j).operator ->();
              recoDigi         = crystalHit->recoCaloDigis().at(0).operator ->();

              indexMC          = 0;//caloDigi.index();

              eDep             = crystalHit->energyDep();
              psd              = 0;//recoDigi  ->psd();

              if (psd >= _psdThreshold){
                crystalTime        =   crystalHit->time();

                if (eDep> 10.){
                  eMeanTot          += eDep;
                  clusterMeanTime   += crystalTime*eDep;
                }

                if (eDep > energyMax){
                  clusterTime       = crystalTime;
                  energyMax         = eDep;

            if (nCaloHitMC > 0) {
              caloDigiMC        = &caloDigiMCCol->at(indexMC);
              clusterMCMeanTime = caloDigiMC->time();
              clusterMCTime     = caloDigiMC->time();
            }
            }
              }

        if (nCaloHitMC > 0) {

          for (unsigned k=0; k<caloDigiMC->nParticles(); ++k){
            sim =   caloDigiMC->energyDeposit(k).sim().operator ->();
            int        pdgId       = sim->pdgId();
            double     ceEnergy    = 104.9;
            double     startEnergy = sim->startMomentum().e();
            if ( (pdgId == PDGCode::e_minus) && (startEnergy>ceEnergy))
              {
                isConversion = 1;
              }
            // if ( sim->fromGenerator() ){
            //   GenParticle* gen = (GenParticle*) &(sim->genParticle());
            //   if ( gen->generatorId().isConversion() ){
            //         isConversion = 1;
            //   }
            // }
          }//end loop on the particles inside the crystalHit
        }
      }
      if (eMeanTot>0){
        clusterMeanTime /= eMeanTot;
      }

      _cluEnergy    [i]     = cluster->energyDep();
      _cluTime      [i]     = clusterTime;
      _cluMeanTime  [i]     = clusterMeanTime;
      _cluMCTime    [i]     = clusterMCTime;
      _cluMCMeanTime[i]     = clusterMCMeanTime;
      _cluCrysE     [i]     = energyMax;
      _cluNcrys     [i]     = nCrystals;
      _cluCogX      [i]     = cluster->cog3Vector().x();
      _cluCogY      [i]     = cluster->cog3Vector().y();
      _cluCogR      [i]     = sqrt(_cluCogX[i]*_cluCogX[i] + _cluCogY[i]*_cluCogY[i]);
      _cluCogZ      [i]     = cluster->diskID();//cog3Vector().z();
      _cluConv      [i]     = isConversion;

    }//end filling calo clusters info

    _Ntup->Fill();
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ReadCaloDigi)
