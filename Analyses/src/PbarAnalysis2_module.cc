//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
// Original author Bertrand Echenard
//

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
//tracker includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
// data
#include "Offline/RecoDataProducts/inc/TrackClusterMatch.hh"


#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"

#include "Offline/RecoDataProducts/inc/TrkCaloMatch.hh"
// data
#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"


#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>



namespace mu2e {


  class PbarAnalysis2 : public art::EDAnalyzer {

     public:

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;



       explicit PbarAnalysis2(fhicl::ParameterSet const& pset);
       virtual ~PbarAnalysis2() { }

       virtual void beginJob();
       virtual void endJob();

       // This is called for each event.
       virtual void analyze(const art::Event& e);





     private:


    int    findBestCluster(TrackClusterMatchCollection const& trackClusterMatches, int trkId, double maxChi2);
    int    findBestTrack  (TrackClusterMatchCollection const& trackClusterMatches, int cluId, double maxChi2);
    float  findBestChi2  (TrackClusterMatchCollection const& trackClusterMatches, int cluId, double maxChi2);


    int _diagLevel;
    bool _doGenerated;
    int _nProcess;

    std::string _g4ModuleLabel;
    std::string _psVacuumStepPoints;
    std::string _generatorModuleLabel;
    art::InputTag _simParticleTag;
    std::string _simParticleModuleLabel;
    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;
    std::string _trackerHitsModuleLabel;
    std::string _trackerHitsInstanceLabel;

    std::string _shLabel;
    std::string _shpLabel;
    std::string _shfLabel;
    std::string _bkfLabel;

    const StrawHitCollection*         _shcol;
    const StrawHitFlagCollection*     _shfcol;
    const StrawHitFlagCollection*     _bkfcol;
    StrawHitFlagCollection*           _flags;
    const StrawHitPositionCollection* _shpcol;

    std::string _caloCrystalModuleLabel;
    std::string _caloClusterModuleLabel;
    std::string _caloDigiTruthModuleLabel;
    std::string _caloClusterTruthModuleLabel;
    std::string _caloClusterAlgorithm;
    std::string _caloClusterSeeding;
    const std::string _producerName;
    std::string _virtualDetectorLabel;
    std::string _stepPointMCLabel;
    std::string _trkFitterModuleLabel;
    std::string _trackClusterMatchModuleLabel;

    double _maxChi2Match;
    bool _writeVertexFile;
    std::string _pbarVertexOutName;
    std::ofstream _pbarVertexOut;
    TTree* _Ntup;


       TH1F *_hcryE,*_hcryT,*_hcryX,*_hcryY,*_hcryZ;
       TH1F *_hcluE,*_hcluT,*_hcluX,*_hcluY,*_hcluZ,*_hcluE1Et,*_hcluE1E9,*_hcluE1E25;

       TH2F *_hxy;


       int   _evt,_run;

       int   _nGen,_genPdgId[16384],_genCrCode[16384];
       float _genmomX[16384],_genmomY[16384],_genmomZ[16384],_genStartX[16384],_genStartY[16384],_genStartZ[16384],_genStartT[16384];

       int   _nHits,_cryId[163840],_crySectionId[163840],_crySimIdx[163840],_crySimLen[163840];
       float _cryEtot,_cryTime[163840],_cryEdep[163840],_cryPosX[163840],_cryPosY[163840],_cryPosZ[163840],_cryLeak[163840];

       int   _nSim,_motId[500000],_motPdgId[500000],_motcrCode[500000],_motGenIdx[500000];
       float _motmom[500000],_motStartX[500000],_motStartY[500000],_motStartZ[500000],_motStartT[500000];
       float _motTime[500000],_motEdep[16348],_motPosX[500000],_motPosY[500000],_motPosZ[500000];

       int   _nCluster,_nCluSim,_cluNcrys[16384];
       float _cluEnergy[16384],_cluTime[16384],_cluCogX[16384],_cluCogY[16384],_cluCogZ[16384],_cluE1[16384],_cluE9[16384],_cluE25[16384],_cluSecMom[16384];
       int   _cluSplit[16384],_cluConv[16384],_cluSimIdx[16384],_cluSimLen[16384];
       std::vector<std::vector<int> > _cluList;

       int   _clusimId[16384],_clusimPdgId[16384],_clusimGenId[16384],_clusimGenPdg[16384],_clusimCrCode[16384];
       float _clusimMom[16384],_clusimMom2[16384],_clusimPosX[16384],_clusimPosY[16384],_clusimPosZ[16384],_clusimStartX[16384],
             _clusimStartY[16384],_clusimStartZ[16384],_clusimTime[16384],_clusimEdep[16384];

       int   _nVd,_vdId[16384],_vdPdgId[16384],_vdenIdx[16384];
       float _vdTime[16384],_vdPosX[16384],_vdPosY[16384],_vdPosZ[16384],_vdMom[16384],_vdMomX[16384],_vdMomY[16384],_vdMomZ[16384];

       int   _nTrk,_trkstat[8192],_trknHit[8192];
       float _trkDip[8192],_trkpt[8192],_trkcon[8192],_trkmomErr[8192];

       float _trkt0[8192],_trkMom[8192], _trkd0[8192],_trkz0[8192],_trkOmega[8192],_trkPhi0[8192],_trkt0Err[8192]; //rhb added


    float _hitTime[32768];
    float _hitEnergy[32768];
    float _hitZ[32768];
    int _hitIsDelta[32768];
    int _nTrackerHits;
    int _nGoodHits;
    int _nDeltas;
    float _bestCluster[8192];
    float _bestChi2[8192];
    int _nTrackClusterMatches;

    float _xStart;
    float _yStart;
    float _zStart;
    float _pxStart;
    float _pyStart;
    float _pzStart;

    float _pxIncomingProton;
    float _pyIncomingProton;
    float _pzIncomingProton;

    float _pxInitialAntiProton;
    float _pyInitialAntiProton;
    float _pzInitialAntiProton;

    float _endGlobalTime;


  };


  PbarAnalysis2::PbarAnalysis2(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _doGenerated(pset.get<bool>("doGenerated",false)),
    _nProcess(pset.get<int>("nProcess",0)),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _psVacuumStepPoints(pset.get<std::string>("psVacuumStepPoints","AntiProtonSteps")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel")),
    _simParticleTag(pset.get<std::string>("simParticleTag")),
    _simParticleModuleLabel(pset.get<std::string>("simParticleModuleLabel")),
    //    _trackerHitsModuleLabel(pset.get<std::string>("trackerHitsModuleLabel")),
    //_trackerHitsInstanceLabel(pset.get<std::string>("trackerHitsInstanceLabel")),
    _shLabel(pset.get<std::string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel    (pset.get<std::string>("StrawHitFlagCollectionLabel" ,"FlagStrawHits"  )),
    _bkfLabel    (pset.get<std::string>("StrawHitFlagCollectionLabel" ,"FlagBkgHits"  )),
    _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _caloDigiTruthModuleLabel(pset.get<std::string>("caloHitTruthModuleLabel")),
    _caloClusterTruthModuleLabel(pset.get<std::string>("caloClusterTruthModuleLabel")),
    _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
    _stepPointMCLabel(pset.get<std::string>("stepPointMCLabel")),
    _trkFitterModuleLabel(pset.get<std::string>("trkFitterModuleLabel")),
    _trackClusterMatchModuleLabel(pset.get<std::string>("trackClusterMatchModuleLabel")),
    _maxChi2Match(pset.get<double>("maxChi2Match")),
    _writeVertexFile(pset.get<bool>("writeVertexFile")),
    _pbarVertexOutName(pset.get<std::string>("pbarVertexOut")),
    _Ntup(0)
  {
  }

  void PbarAnalysis2::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;
       if (_writeVertexFile){
         _pbarVertexOut.open(_pbarVertexOutName.c_str());
       }
       _Ntup  = tfs->make<TTree>("PbarOut2", "PbarOut2");



       _Ntup->Branch("evt",          &_evt ,        "evt/I");
       _Ntup->Branch("run",          &_run ,        "run/I");
       _Ntup->Branch("cryEtot",      &_cryEtot ,    "cryEtot/F");

       _Ntup->Branch("nGen",         &_nGen ,        "nGen/I");
       _Ntup->Branch("genId",        &_genPdgId,     "genId[nGen]/I");
       _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode[nGen]/I");
       _Ntup->Branch("genMomX",      &_genmomX,      "genMomX[nGen]/F");
       _Ntup->Branch("genMomY",      &_genmomY,      "genMomY[nGen]/F");
       _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomZ[nGen]/F");
       _Ntup->Branch("genStartX",    &_genStartX,    "genStartX[nGen]/F");
       _Ntup->Branch("genStartY",    &_genStartY,    "genStartY[nGen]/F");
       _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ[nGen]/F");
       _Ntup->Branch("genStartT",    &_genStartT,    "genStartT[nGen]/F");

       _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
       _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
       _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
       _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
       _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
       _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
       _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
       _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
       _Ntup->Branch("crySimIdx",    &_crySimIdx ,   "crySimIdx[nCry]/I");
       _Ntup->Branch("crySimLen",    &_crySimLen ,   "crySimLen[nCry]/I");

       _Ntup->Branch("nSim",         &_nSim ,        "nSim/I");
       _Ntup->Branch("simId",        &_motId ,       "simId[nSim]/I");
       _Ntup->Branch("simPdgId",     &_motPdgId ,    "simPdgId[nSim]/I");
       _Ntup->Branch("simCrCode",    &_motcrCode ,   "simCrCode[nSim]/I");
       _Ntup->Branch("simMom",       &_motmom ,      "simMom[nSim]/F");
       _Ntup->Branch("simStartX",    &_motStartX ,   "simStartX[nSim]/F");
       _Ntup->Branch("simStartY",    &_motStartY ,   "simStartY[nSim]/F");
       _Ntup->Branch("simStartZ",    &_motStartZ ,   "simStartZ[nSim]/F");
       _Ntup->Branch("simStartT",    &_motStartT ,   "simStartT[nSim]/F");
       _Ntup->Branch("simPosX",      &_motPosX ,     "simPosX[nSim]/F");
       _Ntup->Branch("simPosY",      &_motPosY ,     "simPosY[nSim]/F");
       _Ntup->Branch("simPosZ",      &_motPosZ ,     "simPosZ[nSim]/F");
       _Ntup->Branch("simTime",      &_motTime ,     "simTime[nSim]/F");
       _Ntup->Branch("simEdep",      &_motEdep ,     "simEdep[nSim]/F");
       _Ntup->Branch("simGenIdx",    &_motGenIdx ,   "simGenIdx[nSim]/I");


       _Ntup->Branch("nCluster",     &_nCluster ,    "nCluster/I");
       _Ntup->Branch("cluEnergy",    &_cluEnergy ,   "cluEnergy[nCluster]/F");
       _Ntup->Branch("cluTime",      &_cluTime ,     "cluTime[nCluster]/F");
       _Ntup->Branch("cluCogX",      &_cluCogX ,     "cluCogX[nCluster]/F");
       _Ntup->Branch("cluCogY",      &_cluCogY ,     "cluCogY[nCluster]/F");
       _Ntup->Branch("cluCogZ",      &_cluCogZ ,     "cluCogZ[nCluster]/F");
       _Ntup->Branch("cluNcrys",     &_cluNcrys ,    "cluNcrys[nCluster]/I");
       _Ntup->Branch("cluE1",        &_cluE1 ,       "cluE1[nCluster]/F");
       _Ntup->Branch("cluE9",        &_cluE9 ,       "cluE9[nCluster]/F");
       _Ntup->Branch("cluE25",       &_cluE25 ,      "cluE25[nCluster]/F");
       _Ntup->Branch("cluSecMom",    &_cluSecMom ,   "cluSecMom[nCluster]/F");
       _Ntup->Branch("cluSplit",     &_cluSplit ,    "cluSplit[nCluster]/I");
       _Ntup->Branch("cluConv",      &_cluConv ,     "cluConv[nCluster]/I");
       _Ntup->Branch("cluSimIdx",    &_cluSimIdx ,   "cluSimIdx[nCluster]/I");
       _Ntup->Branch("cluSimLen",    &_cluSimLen ,   "cluSimLen[nCluster]/I");
       _Ntup->Branch("cluList",      &_cluList);

       _Ntup->Branch("nCluSim",      &_nCluSim ,     "nCluSim/I");
       _Ntup->Branch("clusimId",     &_clusimId ,    "clusimId[nCluSim]/I");
       _Ntup->Branch("clusimPdgId",  &_clusimPdgId , "clusimPdgId[nCluSim]/I");
       _Ntup->Branch("clusimGenId",  &_clusimGenId , "clusimGenId[nCluSim]/I");
       _Ntup->Branch("clusimGenPdg", &_clusimGenPdg, "clusimGenPdg[nCluSim]/I");
       _Ntup->Branch("clusimCrCode", &_clusimCrCode ,"clusimCrCode[nCluSim]/I");
       _Ntup->Branch("clusimMom",    &_clusimMom ,   "clusimMom[nCluSim]/F");
       _Ntup->Branch("clusimMom2",   &_clusimMom2 ,  "clusimMom2[nCluSim]/F");
       _Ntup->Branch("clusimPosX",   &_clusimPosX ,  "clusimPosX[nCluSim]/F");
       _Ntup->Branch("clusimPosY",   &_clusimPosY ,  "clusimPosY[nCluSim]/F");
       _Ntup->Branch("clusimPosZ",   &_clusimPosZ ,  "clusimPosZ[nCluSim]/F");
       _Ntup->Branch("clusimStartX", &_clusimStartX ,"clusimStartX[nCluSim]/F");
       _Ntup->Branch("clusimStartY", &_clusimStartY ,"clusimStartY[nCluSim]/F");
       _Ntup->Branch("clusimStartZ", &_clusimStartZ ,"clusimStartZ[nCluSim]/F");
       _Ntup->Branch("clusimTime",   &_clusimTime ,  "clusimTime[nCluSim]/F");
       _Ntup->Branch("clusimEdep",   &_clusimEdep ,  "clusimEdep[nCluSim]/F");

       _Ntup->Branch("nVd",      &_nVd ,     "nVd/I");
       _Ntup->Branch("vdId",     &_vdId ,    "vdId[nVd]/I");
       _Ntup->Branch("vdPdgId",  &_vdPdgId , "vdPdgId[nVd]/I");
       _Ntup->Branch("vdMom",    &_vdMom ,   "vdMom[nVd]/F");
       _Ntup->Branch("vdMomX",   &_vdMomX ,  "vdMomX[nVd]/F");
       _Ntup->Branch("vdMomY",   &_vdMomY ,  "vdMomY[nVd]/F");
       _Ntup->Branch("vdMomZ",   &_vdMomZ ,  "vdMomZ[nVd]/F");
       _Ntup->Branch("vdPosX",   &_vdPosX ,  "vdPosX[nVd]/F");
       _Ntup->Branch("vdPosY",   &_vdPosY ,  "vdPosY[nVd]/F");
       _Ntup->Branch("vdPosZ",   &_vdPosZ ,  "vdPosZ[nVd]/F");
       _Ntup->Branch("vdTime",   &_vdTime ,  "vdTime[nVd]/F");
       _Ntup->Branch("vdGenIdx", &_vdenIdx , "vdGenIdx[nVd]/I");

       _Ntup->Branch("nTrk",         &_nTrk ,        "nTrk/I");
       _Ntup->Branch("trkDip",       &_trkDip ,      "trkDip[nTrk]/F");
       _Ntup->Branch("trkpt",        &_trkpt ,       "trkpt[nTrk]/F");
       _Ntup->Branch("trkstat",      &_trkstat ,     "trkstat[nTrk]/I");
       _Ntup->Branch("trkcon",       &_trkcon ,      "trkcon[nTrk]/F");
       _Ntup->Branch("trknHit",      &_trknHit ,     "trknHit[nTrk]/I");
       _Ntup->Branch("trkmomErr",    &_trkmomErr ,   "trkmomErr[nTrk]/F");
       _Ntup->Branch("trkt0",        &_trkt0 ,       "trkt0[nTrk]/F");
       _Ntup->Branch("trkt0Err",     &_trkt0Err,     "trkt0Err[nTrk]/F");

       _Ntup->Branch("trkMom",       &_trkMom ,      "trkMom[nTrk]/F");
       _Ntup->Branch("trkd0",        &_trkd0 ,       "trkd0[nTrk]/F");
       _Ntup->Branch("trkz0",        &_trkz0 ,       "trkz0[nTrk]/F");
       _Ntup->Branch("trkOmega",     &_trkOmega ,    "trkOmega[nTrk]/F");
       _Ntup->Branch("trkPhi0",      &_trkPhi0 ,     "trkPhi0[nTrk]/F");
       _Ntup->Branch("bestCluster",  &_bestCluster , "bestCluster[nTrk]/I");
       _Ntup->Branch("bestChi2",     &_bestChi2 ,    "bestChi2[nTrk]/F");

       _Ntup->Branch("nTrackerHits",         &_nTrackerHits,         "nTrackerHits/I");
       _Ntup->Branch("nTrackClusterMatches", &_nTrackClusterMatches, "nTrackClusterMatches/I");
       _Ntup->Branch("nGoodHits",            &_nGoodHits,            "nGoodHits/I");
       _Ntup->Branch("nDeltas",              &_nDeltas,              "nDeltas/I");
       _Ntup->Branch("hitTime",              &_hitTime,              "hitTime[nGoodHits]/F");
       _Ntup->Branch("hitEnergy",            &_hitEnergy,            "hitEnergy[nGoodHits]/F");
       _Ntup->Branch("hitZ",                 &_hitZ,                 "hitZ[nGoodHits]/F");
       _Ntup->Branch("hitIsDelta",           &_hitIsDelta,           "hitIsDelta[nGoodHits]/I");

       _Ntup->Branch("pxPbar",                &_pxStart,               "pxPbar/F");
       _Ntup->Branch("pyPbar",                &_pyStart,               "pyPbar/F");
       _Ntup->Branch("pzPbar",                &_pzStart,               "pzPbar/F");

       _Ntup->Branch("xPbar",                &_xStart,               "xPbar/F");
       _Ntup->Branch("yPbar",                &_yStart,               "yPbar/F");
       _Ntup->Branch("zPbar",                &_zStart,               "zPbar/F");

       _Ntup->Branch("pxIncomingProton",                &_pxIncomingProton,               "pxIncomingProton/F");
       _Ntup->Branch("pyIncomingProton",                &_pyIncomingProton,               "pyIncomingProton/F");
       _Ntup->Branch("pzIncomingProton",                &_pzIncomingProton,               "pzIncomingProton/F");

       _Ntup->Branch("pxInitialAntiProton",                &_pxInitialAntiProton,               "pxInitialAntiProton/F");
       _Ntup->Branch("pyInitialAntiProton",                &_pyInitialAntiProton,               "pyInitialAntiProton/F");
       _Ntup->Branch("pzInitialAntiProton",                &_pzInitialAntiProton,               "pzInitialAntiProton/F");


       _hcryE     = tfs->make<TH1F>("cryEdep",  "Energy deposited / crystal", 100,    0., 50.   );
       _hcryT     = tfs->make<TH1F>("cryTime",  "Time of crystal hit",        100,    0., 2000. );
       _hcryX     = tfs->make<TH1F>("cryX",     "X coord of crystal hit",     100,  300., 700.  );
       _hcryY     = tfs->make<TH1F>("cryY",     "Y coord of crystal hit",     100,  300., 700.  );
       _hcryZ     = tfs->make<TH1F>("cryZ",     "Z coord of crystal hit",     100,11000., 13000.);
       _hcluE     = tfs->make<TH1F>("cluEdep",  "Energy deposited / clustal", 150,    0., 150.  );
       _hcluT     = tfs->make<TH1F>("cluTime",  "Time of clustal hit",        100,    0., 2000. );
       _hcluX     = tfs->make<TH1F>("cluX",     "X coord of clustal hit",     100,  300., 700.  );
       _hcluY     = tfs->make<TH1F>("cluY",     "Y coord of clustal hit",     100,  300., 700.  );
       _hcluZ     = tfs->make<TH1F>("cluZ",     "Z coord of clustal hit",     100,11000., 13000.);
       _hcluE1Et  = tfs->make<TH1F>("cluE1Et",  "E1/Etot",                    100,    0., 1.1   );
       _hcluE1E9  = tfs->make<TH1F>("cluE1E9",  "E1/E9",                      100,    0., 1.1   );
       _hcluE1E25 = tfs->make<TH1F>("cluE1E25", "E1/E25",                     100,    0., 1.1   );
       _hxy       = tfs->make<TH2F>("cryxy", "cryxy",                     350,-700,700,350,-700,700  );

  }



  void PbarAnalysis2::endJob(){
  }




  void PbarAnalysis2::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from PbarAnalysis2 =  "<<_nProcess <<std::endl;

      double _mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();


      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloHitCollection> CaloHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, CaloHitsHandle);
      CaloHitCollection const& CaloHits(*CaloHitsHandle);

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);

      //Virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);

      //Calo digi truth assignment
      art::Handle<CaloHitMCTruthAssn> caloDigiTruthHandle;
      event.getByLabel(_caloDigiTruthModuleLabel, caloDigiTruthHandle);
      const CaloHitMCTruthAssn& caloDigiTruth(*caloDigiTruthHandle);

       //Calo cluster truth assignment
      art::Handle<CaloClusterMCTruthAssn> caloClusterTruthHandle;
      event.getByLabel(_caloClusterTruthModuleLabel, caloClusterTruthHandle);
      const CaloClusterMCTruthAssn& caloClusterTruth(*caloClusterTruthHandle);

     // Get tracks
      art::Handle<KalRepPtrCollection> trksHandle;
      event.getByLabel(_trkFitterModuleLabel,trksHandle);
      KalRepPtrCollection const& trksPtrColl(*trksHandle);

     // Get track calorimeter matches
      //     art::Handle<TrkCaloMatchCollection>  trkCaloMatchHandle;
      //event.getByLabel(_trkCaloMatchModuleLabel, trkCaloMatchHandle);
      //TrkCaloMatchCollection const& trkCaloMatches(*trkCaloMatchHandle);
      art::Handle<TrackClusterMatchCollection>  trackClusterMatchHandle;
      event.getByLabel(_trackClusterMatchModuleLabel, trackClusterMatchHandle);
      TrackClusterMatchCollection const& trackClusterMatches(*trackClusterMatchHandle);


      std::map<art::Ptr<SimParticle>, const StepPointMC*> vdMap;
      if (vdhits.isValid())
      {
         for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
         {
            const StepPointMC& hit = *iter;
            if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;
            vdMap[hit.simParticle()] = &hit;
         }
      }

      art::Handle<SimParticleCollection> simParticleHandle;
      event.getByLabel(_simParticleModuleLabel,simParticleHandle);
      bool haveSimPart = simParticleHandle.isValid();
      SimParticleCollection const& simParticles(*simParticleHandle);
      if (haveSimPart) haveSimPart = !(simParticles.empty());

      if (_diagLevel > 1){
        std::cout << "haveSimPart = " << haveSimPart << std::endl;
      }
      if (_writeVertexFile){
        if (haveSimPart){
          for (auto iSim : simParticles){
            //          mu2e::SimParticle& obj = iSim;
            const auto& p = iSim.second;
            //
            // check this is the initial pbar and add what I need to the ntuple.
            // very unfortunately implementation dependent: how many stages
            //
            // set these now, write them and initial proton momentum vector later
            if (p.id().asInt() == 2000001 && p.pdgId() == PDGCode::anti_proton){
              //              if (p.pdgId() == PDGCode::anti_proton){
              _xStart = p.startPosition().x();
              _yStart = p.startPosition().y();
              _zStart = p.startPosition().z();
              _pxStart = p.startMomentum().x();
              _pyStart = p.startMomentum().y();
              _pzStart = p.startMomentum().z();
              _endGlobalTime = p.endGlobalTime();
            }

            if (_diagLevel > 2 && p.id().asInt() == 1000001 && p.pdgId() == PDGCode::anti_proton){
              std::cout<<" id = "<<p.id()
                       <<", parent="<<p.parentId()
                       <<", pdgId="<<p.pdgId()
                       <<", start="<<p.startPosition()
                       <<", pstart="<<p.startMomentum()
                       <<", end="<<p.endPosition()
                       <<", pend="<<p.endMomentum()
                       <<", nSteps = "<<p.nSteps()
                       <<", preLastStepKE = "<<p.preLastStepKineticEnergy()
                       <<", creationCode = "<<p.creationCode()
                       <<", stoppingCode = "<<p.stoppingCode()
                       <<", startG4Status = "<<p.startG4Status()
                       <<", endG4Status = "<<p.endG4Status()
                       <<", _endGlobalTime = " << p.endGlobalTime() << std::endl;
            }
          }
         }
      }
    art::Handle<GenParticleCollection> genParticleHandle;
    event.getByLabel(_generatorModuleLabel, genParticleHandle);
    bool haveGenPart = genParticleHandle.isValid();
    GenParticleCollection const& genParticles(*genParticleHandle);
    if ( haveGenPart ) haveGenPart = !(genParticles.empty());
    CLHEP::HepLorentzVector initialProtonFourMomentum(0.,0.,0.,0.);
    CLHEP::HepLorentzVector initialAntiProtonFourMomentum(0.,0.,0.,0.);
    //
    // look at the genParticles and see that the geantino with the initial proton is still there.
    if (haveGenPart){
     for (auto iGen : genParticles ) {
        if (_diagLevel > 0){
          std::cout << "  in gen particles: particle id " << iGen.pdgId() << " particle momentum " << iGen.momentum()
                    << " position " << iGen.position()
                    << "generator id " << iGen.generatorId()<< std::endl;
        }

         //
        // if reading from file, the gen particle has the incoming pbar info. Use the fake z to distinguish.
        if (iGen.pdgId() == PDGCode::anti_proton && !_writeVertexFile && iGen.position().z() > 1000. ){
          _xStart = iGen.position().x();
          _yStart = iGen.position().y();
          _zStart = iGen.position().z();
          _pxStart = iGen.momentum().x();
          _pyStart = iGen.momentum().y();
          _pzStart = iGen.momentum().z();
          if (_diagLevel > 2){
            std::cout << "initial pbar pos at stopping target = " << _xStart << " " << _yStart << " " << _zStart << std::endl;
            std::cout << "initial pbar mom at stopping target = " << _pxStart << " " << _pyStart << " " << _pzStart << std::endl;
          }
        }
        if (iGen.pdgId() == PDGCode::anti_proton && iGen.position().z() < 1000.){
          initialAntiProtonFourMomentum = iGen.momentum();
          _pxInitialAntiProton = iGen.momentum().x();
          _pyInitialAntiProton = iGen.momentum().y();
          _pzInitialAntiProton = iGen.momentum().z();
          if (_diagLevel > 2){
            std::cout << "pos and momentum of initial pbar = " << iGen.momentum()
                      << " " << iGen.position() << std::endl;

          }
        }
        //
        // the predecessor put the proton back into the event as a genParticle after reading from file
        //
        // protons are a delta fcn at t=0 now.  Can spread these later in later stage job
        // huh?  how can this ever worked, pdgid was zero. also why is genid = 38?  4/13/2018  was checking pdgcode pdg code for proton??
        if (iGen.pdgId() == 0){
          initialProtonFourMomentum = iGen.momentum();
          _pxIncomingProton = iGen.momentum().x();
          _pyIncomingProton = iGen.momentum().y();
          _pzIncomingProton = iGen.momentum().z();
          if (_diagLevel > 2){
            std::cout << "initial beam proton momentum = " << _pxIncomingProton << " " << _pyIncomingProton << " " << _pzIncomingProton << std::endl;
          }
        }

      }

     if (_writeVertexFile){
       _pbarVertexOut <<  _xStart << " "
                      <<  _yStart << " "
                      <<  _zStart << " "
                      << _pxStart << " "
                      << _pyStart << " "
                      << _pzStart << " "
                      << _endGlobalTime << " "
                      << _pxInitialAntiProton << " "
                      << _pyInitialAntiProton << " "
                      << _pzInitialAntiProton << " "
                      << _pxIncomingProton << " "
                      << _pyIncomingProton << " "
                      << _pzIncomingProton << std::endl;
     }
    }





       //--------------------------  Do generated particles --------------------------------


       _evt = event.id().event();
       _run = event.run();

       if (_diagLevel == 3){std::cout << "processing event in calo_example " << _nProcess << " run and event  = " << _run << " " << _evt << std::endl;}



       if (_doGenerated)
       {
           //Get generated particles
           art::Handle<GenParticleCollection> gensHandle;
           event.getByLabel(_generatorModuleLabel, gensHandle);
           GenParticleCollection const& genParticles(*gensHandle);

           _nGen = genParticles.size();

           if (_diagLevel > 2){
             std::cout << "number of generated particles = " << _nGen << std::endl;
           }
           for (unsigned int i=0; i < genParticles.size(); ++i)
           {
               GenParticle const* gen = &genParticles[i];
               _genPdgId[i]   = gen->pdgId();
               _genCrCode[i]  = gen->generatorId().id();
               _genmomX[i]    = gen->momentum().vect().x();
               _genmomY[i]    = gen->momentum().vect().y();
               _genmomZ[i]    = gen->momentum().vect().z();
               _genStartX[i]  = gen->position().x();
               _genStartY[i]  = gen->position().y();
               _genStartZ[i]  = gen->position().z();
               _genStartT[i]  = gen->time();
               if (_diagLevel > 2){
                 std::cout << "pdgId = " << gen->pdgId() << " "
                           << "generator Id = " << gen->generatorId().id() << "\n "
                           << "momentum     = " << gen->momentum() << "\n "
                           << "position     = " << gen->position() << "\n "
                           << "time         = " << gen->time()
                           << std::endl;
               }
           }
       } else {_nGen=0;}


       //--------------------------  Do calorimeter hits --------------------------------

       _nHits = _nSim = 0;
       _cryEtot = 0.0;

       for (unsigned int ic=0; ic<CaloHits.size();++ic)
       {
           const CaloHit &hit     = CaloHits.at(ic);
           int diskId                    = cal.crystal(hit.crystalID()).diskID();
           CLHEP::Hep3Vector crystalPos  = cal.crystal(hit.crystalID()).localPosition();  //in disk FF frame

           auto itMC = caloDigiTruth.begin();
           while (itMC != caloDigiTruth.end()) {if (itMC->first.get() == &hit) break; ++itMC;}
           unsigned nCrySims = (itMC != caloDigiTruth.end()) ? itMC->second->nParticles() : 0;

           _cryEtot             += hit.energyDep();
           _cryTime[_nHits]      = hit.time();
           _cryEdep[_nHits]      = hit.energyDep();
           _cryPosX[_nHits]      = crystalPos.x();
           _cryPosY[_nHits]      = crystalPos.y();
           _cryPosZ[_nHits]      = crystalPos.z();
           _cryId[_nHits]        = hit.crystalID();
           _crySectionId[_nHits] = diskId;

           _crySimIdx[_nHits]    = _nSim;
           _crySimLen[_nHits]    = nCrySims;

           for (unsigned i=0;i< nCrySims;++i)
           {
               const auto& eDepMC = itMC->second->energyDeposit(i);

               auto parent(eDepMC.sim());
               while ( parent->hasParent()) parent = parent->parent();

               _motId[_nSim]      = eDepMC.sim()->id().asInt();
               _motPdgId[_nSim]   = eDepMC.sim()->pdgId();
               _motmom[_nSim]     = eDepMC.momentumIn();
               _motcrCode[_nSim]  = eDepMC.sim()->creationCode();
                      _motTime[_nSim]    = eDepMC.time();
               _motEdep[_nSim]    = eDepMC.energyDep();

               _motStartX[_nSim]  = parent->startPosition().x();
               _motStartY[_nSim]  = parent->startPosition().y();
               _motStartZ[_nSim]  = parent->startPosition().z();

               ++_nSim;
            }


           _hcryE->Fill(hit.energyDep());
           _hcryT->Fill(hit.time());
           _hcryX->Fill(crystalPos.x());
           _hcryY->Fill(crystalPos.y());
           _hcryZ->Fill(crystalPos.z());
           _hxy->Fill(crystalPos.x(),crystalPos.y(),hit.energyDep());
           ++_nHits;
       }


       // do tracks, old school
      //--------------------------  Do tracks  --------------------------------

       _nTrk = 0;

       //       std::cout << "size of trks collection = " << trksPtrColl.size() << std::endl;
       for ( size_t itrk=0; itrk< trksPtrColl.size(); ++itrk )
       {
         KalRep const* krep = trksPtrColl.at(itrk).get();

         // Arc length from center of tracker to the most upstream point on the fitted track.
         double s0 = krep->startValidRange();

         // Momentum and position at s0.
         CLHEP::Hep3Vector p0     = krep->momentum(s0);
         HepPoint          pos0   = krep->position(s0);

         // Some other properties at s0.
         double loclen(0.);
         const TrkSimpTraj* ltraj = krep->localTrajectory(s0,loclen);
         HepVector momvec(3);
         momvec = p0.unit();
         BbrVectorErr momCov = krep->momentumErr(s0);
         double fitMomErr    = sqrt(momCov.covMatrix().similarity(momvec));
         double tanDip       = ltraj->parameters()->parameter()[HelixTraj::tanDipIndex];
         double omega        = ltraj->parameters()->parameter()[HelixTraj::omegaIndex];
         double d0           = ltraj->parameters()->parameter()[HelixTraj::d0Index];
         double z0           = ltraj->parameters()->parameter()[HelixTraj::z0Index];
         double trkPhi0      = ltraj->parameters()->parameter()[HelixTraj::phi0Index];
         double fitCon       = krep->chisqConsistency().significanceLevel();
         double trkt0Err     = krep->t0().t0Err();
         double fitmompt = p0.mag()*(1.0-p0.cosTheta()*p0.cosTheta());

          _trkDip[_nTrk] = tanDip;
          _trkpt[_nTrk]  = fitmompt;
          _trkstat[_nTrk]  = krep->fitStatus().success();
          _trkcon[_nTrk]  = fitCon;
          _trknHit[_nTrk]  = krep->nActive();
          _trkmomErr[_nTrk]  = fitMomErr ;
          _trkt0[_nTrk] = krep->t0().t0();
          _trkMom[_nTrk] = p0.mag();
          _trkd0[_nTrk] = d0;
          _trkz0[_nTrk] = z0;
          _trkOmega[_nTrk] = omega;
          _trkPhi0[_nTrk] = trkPhi0;
          _trkt0Err[_nTrk] = trkt0Err;
          ++_nTrk;

       }




       //--------------------------  Do tracks  --------------------------------
       // we just want the cluster that matches the track.
       // For information about the track at the entrance of the calorimeter, look at the
       // Analysis/src/ReadTrackCaloMatching_module.cc mdoule

       //giani promises this index goes with the track index

       if (_diagLevel > 0)     {std::cout << "size of match collection = " << trackClusterMatches.size() << std::endl;}
       _nTrackClusterMatches = trackClusterMatches.size();
       for (unsigned int i=0;i< trackClusterMatches.size(); ++i)
       {
         //KalRepPtr const& tkrPtr = trksPtrColl.at(i);
         _bestCluster[i] = findBestCluster(trackClusterMatches,i,_maxChi2Match);
         _bestChi2   [i] = findBestChi2   (trackClusterMatches,i,_maxChi2Match);
        }



       //--------------------------  Do clusters --------------------------------
       _nCluster = _nCluSim = 0;
       _cluList.clear();
       for (unsigned int ic=0; ic<caloClusters.size();++ic)
       {
          const CaloCluster& cluster = caloClusters.at(ic);
          std::vector<int> cryList;
          for (auto cryPtr : cluster.caloHitsPtrVector()) cryList.push_back(int(cryPtr.get()- &CaloHits.at(0)));

          //Find the caloDigiMC in the truth map
          auto itMC = caloClusterTruth.begin();
          while (itMC != caloClusterTruth.end()) {if (itMC->first.get() == &cluster) break; ++itMC;}
          const auto eDepMCs = (itMC != caloClusterTruth.end()) ? itMC->second->energyDeposits() : std::vector<CaloEDepMC>{};

          bool isConversion(false);
          if (itMC != caloClusterTruth.end())
          {
             for (auto& edep : eDepMCs)
             {
                auto parent(edep.sim());
                while (parent->hasParent()) parent = parent->parent();
                if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) isConversion=true;
             }
          }
          ClusterUtils cluUtil(cal, cluster);

           _cluEnergy[_nCluster] = cluster.energyDep();
           _cluTime[_nCluster]   = cluster.time();
           _cluNcrys[_nCluster]  = cluster.size();
           _cluCogX[_nCluster]   = cluster.cog3Vector().x(); //in disk FF frame
           _cluCogY[_nCluster]   = cluster.cog3Vector().y();

           //FF as in CaloExample isn't really useful; replace with diskId
           _cluCogZ[_nCluster]   = cluster.diskID();
           _cluE1[_nCluster]     = cluUtil.e1();
           _cluE9[_nCluster]     = cluUtil.e9();
           _cluE25[_nCluster]    = cluUtil.e25();
           _cluSecMom[_nCluster] = cluUtil.secondMoment();
           _cluSplit[_nCluster]  = cluster.isSplit();
           _cluConv[_nCluster]   = isConversion;
           _cluList.push_back(cryList);

           _cluSimIdx[_nCluster] = _nCluSim;
           _cluSimLen[_nCluster] = eDepMCs.size();

           for (unsigned i=0;i< eDepMCs.size();++i)
           {
               const auto& eDepMC = eDepMCs[i];
               art::Ptr<SimParticle> sim = eDepMC.sim();

               art::Ptr<SimParticle> smother(sim);
               while (smother->hasParent() && !smother->genParticle() ) smother = smother->parent();
               int genId=-1;
               if (smother->genParticle()) genId = smother->genParticle()->generatorId().id();
               int genPdg=-1;
               if (smother->genParticle()) genPdg = smother->genParticle()->pdgId();


               double simMom(-1);
               CLHEP::Hep3Vector simPos(0,0,0);
               auto vdMapEntry = vdMap.find(sim);
               if (vdMapEntry != vdMap.end())
               {
                  simMom = vdMapEntry->second->momentum().mag();
                  CLHEP::Hep3Vector simPos = cal.geomUtil().mu2eToDiskFF(cluster.diskID(), vdMapEntry->second->position());
               }

               _clusimId[_nCluSim]     = sim->id().asInt();
               _clusimPdgId[_nCluSim]  = sim->pdgId();
               _clusimGenId[_nCluSim]  = genId;
               _clusimGenPdg[_nCluSim] = genPdg;
               _clusimCrCode[_nCluSim] = sim->creationCode();
               _clusimTime[_nCluSim]   = eDepMC.time();
               _clusimEdep[_nCluSim]   = eDepMC.energyDep();
               _clusimMom[_nCluSim]    = eDepMC.momentumIn();
               _clusimMom2[_nCluSim]   = simMom;
               _clusimPosX[_nCluSim]   = simPos.x(); // in disk FF frame
               _clusimPosY[_nCluSim]   = simPos.y();
               _clusimPosZ[_nCluSim]   = simPos.z();
               _clusimStartX[_nCluSim] = sim->startPosition().x(); // in disk FF frame
               _clusimStartY[_nCluSim] = sim->startPosition().y();
               _clusimStartZ[_nCluSim] = sim->startPosition().z();

               ++_nCluSim;
            }

           _hcluE->Fill(cluster.energyDep());
           _hcluT->Fill(cluster.time());
           _hcluX->Fill(cluster.cog3Vector().x());
           _hcluY->Fill(cluster.cog3Vector().y());
           _hcluZ->Fill(cluster.cog3Vector().z());


           ++_nCluster;
       }




       //--------------------------  Do virtual detectors --------------------------------
       //73/74/77/78 front back inner outer edges disk 0
       //75/76/79/80 front back inner outer edges disk 1
       _nVd = 0;
       if (vdhits.isValid())
       {
           for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
             {
               const StepPointMC& hit = *iter;

               if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;

               double hitTime         = fmod(hit.time(),_mbtime);

               CLHEP::Hep3Vector VDPos = cal.geomUtil().mu2eToTracker(hit.position());

               _vdId[_nVd]    = hit.volumeId();
               _vdPdgId[_nVd] = hit.simParticle()->pdgId();
               _vdTime[_nVd]  = hitTime;
               _vdPosX[_nVd]  = VDPos.x(); //tracker frame
               _vdPosY[_nVd]  = VDPos.y();
               _vdPosZ[_nVd]  = VDPos.z();
               _vdMomX[_nVd]  = hit.momentum().x();
               _vdMomY[_nVd]  = hit.momentum().y();
               _vdMomZ[_nVd]  = hit.momentum().z();
               _vdenIdx[_nVd] = hit.simParticle()->generatorIndex();
               ++_nVd;
             }
         }



       //--------------------------  Do tracker hits  --------------------------------

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (event.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in CaloExample::findData: StrawHitCollection with label=%s not found.\n",
             _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (event.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in CaloExample::findData: StrawHitPositionCollection with label=%s not found.\n",
             _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (event.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in CaloExample::findData: StrawHitFlagCollection with label=%s not found.\n",
             _shfLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> bkflagH;
    if (event.getByLabel(_bkfLabel,bkflagH)) {
      _bkfcol = bkflagH.product();
    }
    else {
      _bkfcol = 0;
      printf(" >>> ERROR in CaloExample::findData: StrawHitFlagCollection with label=%s not found.\n",
             _bkfLabel.data());
    }


    _nTrackerHits = _shcol->size();
    _nGoodHits = 0;
    _nDeltas = 0;

    //   std::cout << "size of flag collection = " << _shfcol->size() << std::endl;
    //std::cout << "_ntrackerHits = " << _nTrackerHits << std::endl;


    for (int ithHit =0; ithHit<_nTrackerHits; ++ithHit)
      {
        //        std::cout << "inside hit loop " << _shcol->at(ithHit).time() << std::endl;

        if (_shcol->at(ithHit).time() > 650. && _nGoodHits< 32767)
          {
            _hitTime[_nGoodHits] = _shcol->at(ithHit).time();
            _hitEnergy[_nGoodHits] = _shcol->at(ithHit).energyDep();
            _hitZ[_nGoodHits] = _shpcol->at(ithHit).pos().z();
            //            bool deltaRay = _bkfcol->at(ithHit).hasAllProperties(StrawHitFlag::delta);
            //                        bool deltaRay = _shpcol->at(ithHit).flag().hasAllProperties(StrawHitFlag::delta);
            //            if (deltaRay){
            //              _hitIsDelta[_nGoodHits] = 1;
              //              ++_nDeltas;
            //            } else{_hitIsDelta[_nGoodHits] = 0;}
            //if (deltaRay){std::cout << "found a delta" << std::endl;}
            //std::cout << "time, energy,z " << _shcol->at(ithHit).time() << " " << _shcol->at(ithHit).energyDep() << " " <<
            //                      _shpcol->at(ithHit).pos().z()
            //           << std::endl;
            ++_nGoodHits;
          }
      }

    //std::cout << " number of good hits = " << _nGoodHits << std::endl;




        _Ntup->Fill();
  }

   int PbarAnalysis2::findBestCluster(TrackClusterMatchCollection const& trackClusterMatches, int trkId, double maxChi2)
   {
     double chi2Best(maxChi2);
     int cluBest(-1);
      if (_diagLevel>0){std::cout << "\n \n " <<std::endl;}
      for (auto const& trkClusterMatch: trackClusterMatches)
      {
               if (trkClusterMatch.iex()==trkId && trkClusterMatch.chi2() < chi2Best)
         {
           if (_diagLevel > 0) {
             std::cout << " chi2Best, chi2, cluBest = " << chi2Best << " " << trkClusterMatch.chi2()  << " " << cluBest << std::endl;
             std::cout << " and the du,dv = " << trkClusterMatch.du() << " " << trkClusterMatch.dv() << std::endl;
           }
           chi2Best= trkClusterMatch.chi2();
           cluBest = trkClusterMatch.icl();
         }
      }
      if (_diagLevel>0){std::cout << "best chi2 = " << chi2Best <<std::endl;}

      return cluBest;
   }


   float PbarAnalysis2::findBestChi2(TrackClusterMatchCollection const& trackClusterMatches, int trkId, double maxChi2)
   {
     double chi2Best(maxChi2);
     int cluBest(-1);
      if (_diagLevel>0){std::cout << "\n \n " <<std::endl;}
      for (auto const& trkClusterMatch: trackClusterMatches)
      {
               if (trkClusterMatch.iex()==trkId && trkClusterMatch.chi2() < chi2Best)
         {
           if (_diagLevel > 0) {
             std::cout << " chi2Best, chi2, cluBest = " << chi2Best << " " << trkClusterMatch.chi2()  << " " << cluBest << std::endl;
             std::cout << " and the du,dv = " << trkClusterMatch.du() << " " << trkClusterMatch.dv() << std::endl;
           }
           chi2Best= trkClusterMatch.chi2();
           cluBest = trkClusterMatch.icl();
         }
      }
      if (_diagLevel>0){std::cout << "best chi2 = " << chi2Best <<std::endl;}

      return static_cast<float>(chi2Best);
   }

   int PbarAnalysis2::findBestTrack(TrackClusterMatchCollection const& trackClusterMatches, int cluId, double maxChi2)
   {
      double chi2Best(maxChi2), trkBest(-1);
      for (auto const& trkClusterMatch: trackClusterMatches)
      {
         if (trkClusterMatch.icl()==cluId && trkClusterMatch.chi2() < chi2Best)
         {
           chi2Best= trkClusterMatch.chi2();
           trkBest = trkClusterMatch.iex();
         }
      }
      return trkBest;
   }


}

DEFINE_ART_MODULE(mu2e::PbarAnalysis2)
