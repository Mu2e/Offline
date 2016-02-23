//
// An EDAnalyzer module that writes out calorimeter and tracker information so that it can be used to study combined track/calo information.  Largely based on Echenard's CaloExample_module.cc


//
// Original author Robert Bernstein
//

#include "CLHEP/Units/SystemOfUnits.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "art/Framework/Core/EDAnalyzer.h"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CalorimeterPhysicalConstants.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"


#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CaloCluster/inc/CaloContentMC.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
//
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"

#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"


// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"

#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>


#include "GeneralUtilities/inc/sqrtOrThrow.hh"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using CLHEP::keV;



namespace mu2e {


  class KineticFracAnalysis : public art::EDAnalyzer {

  public:

    typedef art::Ptr<StepPointMC> StepPtr;
    typedef std::vector<StepPtr>  StepPtrs;
    typedef std::map<int,StepPtrs > HitMap;



    explicit KineticFracAnalysis(fhicl::ParameterSet const& pset);

    virtual ~KineticFracAnalysis() { }

    virtual void beginJob();
    virtual void endJob();

    // This is called for each event.
    virtual void analyze(const art::Event& e);





  private:

    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
    typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;
    typedef art::Ptr<SimParticle> SimParticlePtr;


    int _diagLevel;
    int _nProcess;
    std::string _g4ModuleLabel;
    std::string _generatorModuleLabel;
    art::InputTag _simParticleTag;
    std::string _caloReadoutModuleLabel;
    std::string _caloCrystalModuleLabel;
    std::string _caloHitMCCrystalPtrLabel;
    std::string _caloClusterModuleLabel;
    std::string      _trkCaloMatchModuleLabel;
    std::string      _trkIntersectModuleLabel;
    std::string      _trkFitterModuleLabel;
    std::string _virtualDetectorLabel;
    std::string _stepPointMCLabel;
    std::string _trkPatRecModuleLabel;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    SimParticleTimeOffset _toff;  // time offset smearing
    std::string _shLabel;
    std::string _shpLabel; 
    std::string _shfLabel;
    std::string _bkfLabel;
    std::string _trkfitInstanceName; 
    std::string _caloClusterAlgorithm;
    std::string _caloClusterSeeding;
    const std::string _producerName;
    std::string _instanceName;

 
    const StrawHitCollection*         _shcol;
    const StrawHitFlagCollection*     _shfcol;
    const StrawHitFlagCollection*     _bkfcol;
    const StrawHitPositionCollection* _shpcol;



    TH1F *_hcryE,*_hcryT,*_hcryX,*_hcryY,*_hcryZ;
    TH1F *_hcluE,*_hcluT,*_hcluX,*_hcluY,*_hcluZ,*_hcluE1Et,*_hcluE1E9,*_hcluE1E25;


    TH1F *_hfitstatus
      ,*_htrkcon
      ,*_htrkcon1
      ,*_htrkcon2
      ,*_htrkmomErr
      ,*_htrkmomErr1
      ,*_htrkmomErr2
      ,*_htrkt0
      ,*_hTimeDiffBest1st
      ,*_hTimeDiffBest2nd
      ,*_htrkt0Err
      ,*_htrkt0Err1
      ,*_htrkt0Err2
      ,*_hfitpar_d0
      ,*_hfitpar_z0
      ,*_hfitpar_td
      ,*_hfitpar_om
      ,*_hfitpar_d0omega
      ,*_hPhi0
      ,*_hnactive
      ,*_hfitmom
      ,*_hevtwt
      ,*_hfitmomCutSetC
      ,*_hfitmomFinal
      ,*_hClusterEnergy
      ,*_hNumberOfClusters
      ,*_hClusterEnergyInWindow
      ,*_hNumberOfClustersInWindow
      ,*_hTimeOfCluster
      ,*_hHighestEnergyCluster
      ,*_hNumberOfTracks
      ,*_hMomentumOfTrackForCal
      ,*_hTotalEnergyInWindow
      ,*_hDistX
      ,*_hDistY
      ,*_htime
      ,*_hShortestDistance, *_hShortestFirst
      ,*_hClusterRadius
      ,*_hEnergyAtLargeRadius
      ,*_hNumberAtLargeRadius
      ,*_hTotalEnergyAtLargeRadius
      ,*_hPdgIdOfAllClustersAtLargeRadius
      ,*_hCloseDeltasMean
      ,*_hCloseDeltasSigma
      ,*_hCloseDeltasAmplitude
      ,*_hTimeOfDeltaWrtT0Zoom
      ,*_hEnergyAtSmallRadius
      ,*_hNumberAtSmallRadius
      ,*_hTotalEnergyAtSmallRadius
      ,*_hPdgIdOfAllClustersAtSmallRadius;

    TH1F 
    *_hPdgId
      ,*_hCluConv
      ,*_hEoverP1st
      ,*_hEoverP2nd
      ,*_hEOverPZoom
      ,*_hDistToTrack
      ,*_hNumberOfHits
      ,*_hNumberOfHitsAfter650Nsec
      ,*_hTimeDiff
      ,*_hNHitsInTime
      ,*_hZOfHit
      ,*_hEnergyDeposit
      ,*_hEnergyDepositFront
      ,*_hNumberOfDeltas
      ,*_hTimeOfDeltaWrtT0
      ,*_hCloseDeltas
      ,*_hMuonSampleMomentum
      ,*_hMuonFrac
      ,*_hPionFrac
      ,*_hMuonFracCutMomentum
      ,*_hPionFracCutMomentum
      ,*_hWeirdFracCutMomentum
      ,*_hWeirdFracCluster
      ,*_hMuonFracCluster
      ,*_hPionFracCluster
      ,*_hNumberOfCrystalsMedFrac
      ,*_hNumberOfCrystalsLowFrac
      ,*_hNumberOfCrystalsHighFrac
      ,*_hSingleCrystalFrac
      ,*_hDoubleCrystalFrac
      ,*_hTripleCrystalFrac
      ,*_hCascade
      ,*_hShortestDistanceClusterTime;

    TProfile 
       *_hErrVsMomMuon
      ,*_hErrVsD0Muon
      ,*_hErrVsOmegaMuon
      ,*_hErrVsTanDipMuon
      ,*_hErrVsMomPion
      ,*_hErrVsD0Pion
      ,*_hErrVsOmegaPion
      ,*_hErrVsTanDipPion;
       TH2F *_hTanDipVsMomentumPion;
       TH2F *_hTanDipVsMomentumMuon;

      TProfile
       *_hErrVsMomHiFrac
      ,*_hErrVsD0HiFrac
      ,*_hErrVsOmegaHiFrac
      ,*_hErrVsTanDipHiFrac;

    TH2F *_hTanDipVsMomentumHiFrac, *_hTanDipVsMomentumLoFrac;


    TProfile *_hErrVsMomLoFrac, *_hErrVsD0LoFrac, *_hErrVsOmegaLoFrac, 
      *_hErrVsTanDipLoFrac;

    TH2F *_hFracVsTime, *_hRadiusVsTime, *_hFracVsRadius, *_hFracVsLocation, 
      *_hMuonFracVsLocation, *_hFracVsP;

    int _timeDiff;

    TTree* _Ntup;


    int   _evt,_run;

    int   _nGen,_genPdgId[16384],_genCrCode[16384];
    float _genmomX[16384],_genmomY[16384],_genmomZ[16384],_genStartX[16384],_genStartY[16384],_genStartZ[16384],_genStartT[16384];

    int   _nHits,_cryId[16384],_crySectionId[16384],_crySimIdx[16384],_crySimLen[16384];
    float _cryEtot,_cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384],_cryLeak[16384];

    int   _nSim,_motId[16384],_motPdgId[16384],_motcrCode[16384],_motGenIdx[16384];
    float _motmom[16384],_motStartX[16384],_motStartY[16384],_motStartZ[16384],_motStartT[16384];
    float _motTime[16384],_motEdep[16348],_motPosX[16384],_motPosY[16384],_motPosZ[16384];

    int   _nCluster,_nCluSim,_cluNcrys[16384];
    float _cluEnergy[16384],_cluTime[16384],_cluCogX[16384],_cluCogY[16384],_cluCogZ[16384];
    int   _cluConv[16384],_cluSimIdx[16384],_cluSimLen[16384];
    std::vector<std::vector<int> > _cluList;

    int   _clusimId[16384],_clusimPdgId[16384],_clusimGenIdx[16384];
    float _clusimMom[16384],_clusimPosX[16384],_clusimPosY[16384],_clusimPosZ[16384],_clusimTime[16384],_clusimEdep[16384];

    int   _nVd,_vdId[16384],_vdPdgId[16384],_vdenIdx[16384];
    float _vdTime[16384],_vdPosX[16384],_vdPosY[16384],_vdPosZ[16384],_vdMom[16384];

    int   _nTrkOk,_nTrk,_trkOk[8192],_trkstat[8192],_trknHit[8192];
    float _trkDip[8192],_trkpt[8192],_trkcon[8192],_trkmomErr[8192];

    float _trkt0[8192],_trkMom[8192], _trkd0[8192],_trkz0[8192],_trkOmega[8192],_trkPhi0[8192],_trkt0Err[8192]; //rhb added


    float _hitTime[32768];
    float _hitEnergy[32768];
    float _hitZ[32768];
    int _hitIsDelta[32768];
    int _nTrackerHits;
    int _nGoodHits;
    int _nDeltas;
       


          int    _nMatch, _mTrkId[1024],_mCluId[1024];
          float  _mChi2[1024],_mChi2Pos[1024],_mChi2Time[1024];
    int    _nTrkHel,_trknHitHel[1024],_trkStatHel[1024],_trkCluIdxHel[1024];
          float  _trkx[1024],_trky[1024],_trkz[1024],_trkFFx[1024],_trkFFy[1024],_trkFFz[1024],_trke[1024],_trkt[1024];
          float  _trkpx[1024],_trkpy[1024],_trkpz[1024],_trkprob[1024];
          float  _trkphi0[1024],_trkomega[1024],_trkcdip[1024],_trkdlen[1024];


    int  _goodKfrac;
    //    CLHEP::Hep3Vector _firstDiskLoc;
    //CLHEP::Hep3Vector _secondDiskLoc;
    double _firstDiskZ;
    double _secondDiskZ;
    double _eCritPos;
    double _eCritNeg;
    double _density;
    double _tol = 1.0*CLHEP::mm; // for comparisons, slight differences in computations at 1e-9 level
  };


  KineticFracAnalysis::KineticFracAnalysis(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _simParticleTag(pset.get<string>("simParticleTag")),
    _caloReadoutModuleLabel(pset.get<string>("caloReadoutModuleLabel")),
    _caloCrystalModuleLabel(pset.get<string>("caloCrystalModuleLabel")),
    _caloHitMCCrystalPtrLabel(pset.get<string>("calorimeterHitMCCrystalPtr")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _trkCaloMatchModuleLabel(pset.get<std::string>("trkCaloMatchModuleLabel")),
    _trkIntersectModuleLabel(pset.get<std::string>("trkIntersectModuleLabel")),
    _trkFitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
    _virtualDetectorLabel(pset.get<string>("virtualDetectorName")),
    _stepPointMCLabel(pset.get<string>("stepPointMCLabel")),
    _trkPatRecModuleLabel(pset.get<string>("trkPatRecModuleLabel")),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel(pset.get<string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel    (pset.get<std::string>("StrawHitFlagCollectionLabel" ,"FlagStrawHits"  )),
    _bkfLabel    (pset.get<std::string>("StrawHitBkgFlagCollectionLabel" ,"FlagBkgHits"  )),
    _timeDiff(pset.get<int>("timeDifference",15)),
    _Ntup(0)
 {
    _instanceName = _fdir.name() + _tpart.name();
            _trkfitInstanceName = _fdir.name() + _tpart.name();
  
 }

  /*            art::EDAnalyzer(pset),
            _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
            _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
            _trkCaloMatchModuleLabel(pset.get<std::string>("trkCaloMatchModuleLabel")),
            _trkIntersectModuleLabel(pset.get<std::string>("trkIntersectModuleLabel")),
            _trkFitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
            _tpart((TrkParticle::type)(pset.get<int>("fitparticle"))),
            _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection"))),
            _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
            _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
            _Ntup(0)
          {
             _trkfitInstanceName = _fdir.name() + _tpart.name();
          }
  */
  void KineticFracAnalysis::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("KineticFrac", "KineticFrac");



    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");
    _Ntup->Branch("cryEtot",      &_cryEtot ,    "cryEtot/F");
 
     _Ntup->Branch("nMatch",    &_nMatch ,   "nMatch/I");
       _Ntup->Branch("mTrkId",    &_mTrkId,    "mTrkId[nMatch]/I");
       _Ntup->Branch("mCluId",    &_mCluId,    "mCluId[nMatch]/I");
       _Ntup->Branch("mChi2",     &_mChi2,     "mChi2[nMatch]/F");
       _Ntup->Branch("mChi2Pos",  &_mChi2Pos,  "mChi2Pos[nMatch]/F");
       _Ntup->Branch("mChi2Time", &_mChi2Time, "mChi2Time[nMatch]/F");



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
    _Ntup->Branch("cryDose",      &_cryDose ,     "cryDose[nCry]/F");
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
    _Ntup->Branch("cluConv",      &_cluConv ,     "cluConv[nCluster]/I");
    _Ntup->Branch("cluSimIdx",    &_cluSimIdx ,   "cluSimIdx[nCluster]/I");
    _Ntup->Branch("cluSimLen",    &_cluSimLen ,   "cluSimLen[nCluster]/I");
    _Ntup->Branch("cluList",      &_cluList);

    _Ntup->Branch("nCluSim",      &_nCluSim ,     "nCluSim/I");
    _Ntup->Branch("clusimId",     &_clusimId ,    "clusimId[nCluSim]/I");
    _Ntup->Branch("clusimPdgId",  &_clusimPdgId , "clusimPdgId[nCluSim]/I");
    _Ntup->Branch("clusimGenIdx", &_clusimGenIdx ,"clusimGenIdx[nCluSim]/I");
    _Ntup->Branch("clusimMom",    &_clusimMom ,   "clusimMom[nCluSim]/F");
    _Ntup->Branch("clusimPosX",   &_clusimPosX ,  "clusimPosX[nCluSim]/F");
    _Ntup->Branch("clusimPosY",   &_clusimPosY ,  "clusimPosY[nCluSim]/F");
    _Ntup->Branch("clusimPosZ",   &_clusimPosZ ,  "clusimPosZ[nCluSim]/F");
    _Ntup->Branch("clusimTime",   &_clusimTime ,  "clusimTime[nCluSim]/F");
    _Ntup->Branch("clusimEdep",   &_clusimEdep ,  "clusimEdep[nCluSim]/F");

    _Ntup->Branch("nVd",      &_nVd ,     "nVd/I");
    _Ntup->Branch("vdId",     &_vdId ,    "vdId[nVd]/I");
    _Ntup->Branch("vdPdgId",  &_vdPdgId , "vdPdgId[nVd]/I");
    _Ntup->Branch("vdMom",    &_vdMom ,   "vdMom[nVd]/F");
    _Ntup->Branch("vdPosX",   &_vdPosX ,  "vdPosX[nVd]/F");
    _Ntup->Branch("vdPosY",   &_vdPosY ,  "vdPosY[nVd]/F");
    _Ntup->Branch("vdPosZ",   &_vdPosZ ,  "vdPosZ[nVd]/F");
    _Ntup->Branch("vdTime",   &_vdTime ,  "vdTime[nVd]/F");
    _Ntup->Branch("vdGenIdx", &_vdenIdx , "vdGenIdx[nVd]/I");

    _Ntup->Branch("nTrkOk",       &_nTrkOk ,      "nTrkOk/I");
    _Ntup->Branch("nTrk",         &_nTrk ,        "nTrk/I");
    _Ntup->Branch("trkDip",       &_trkDip ,      "trkDip[nTrk]/F");
    _Ntup->Branch("trkOk",        &_trkOk ,       "trkOk[nTrk]/I");
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

    _hcryE = tfs->make<TH1F>("cryEdep","Energy deposited / crystal",100,0.,50.);
    _hcryT = tfs->make<TH1F>("cryTime","Time of crystal hit",100,0.,2000.);
    _hcryX = tfs->make<TH1F>("cryX","X coord of crystal hit",100,300.,700.);
    _hcryY = tfs->make<TH1F>("cryY","Y coord of crystal hit",100,300.,700.);
    _hcryZ = tfs->make<TH1F>("cryZ","Z coord of crystal hit",100,11000.,13000.);

    _hcluE = tfs->make<TH1F>("cluEdep","Energy deposited / clustal",150,0.,150.);
    _hcluT = tfs->make<TH1F>("cluTime","Time of clustal hit",100,0.,2000.);
    _hcluX = tfs->make<TH1F>("cluX","X coord of clustal hit",100,300.,700.);
    _hcluY = tfs->make<TH1F>("cluY","Y coord of clustal hit",100,300.,700.);
    _hcluZ = tfs->make<TH1F>("cluZ","Z coord of clustal hit",100,11000.,13000.);
    _hcluE1Et = tfs->make<TH1F>("cluE1Et","E1/Etot",100,0,1.1);
    _hcluE1E9 = tfs->make<TH1F>("cluE1E9","E1/E9",100,0,1.1);
    _hcluE1E25 = tfs->make<TH1F>("cluE1E25","E1/E25",100,0,1.1);

    _hPdgId = tfs->make<TH1F>("_hpdgId","pdg Id for cluster associated with fitted track",100,-50.,50.);
    _hfitstatus    = tfs->make<TH1F>("_hfitstatus","fitstatus",20,0.,20.);
    _hnactive  = tfs->make<TH1F>("_hnactive","trknHit[ithtrk]",100,0.,100.);
    _htrkcon       = tfs->make<TH1F>("_htrkcon","trkcon",100,0.,1.);
    _htrkcon1       = tfs->make<TH1F>("_htrkcon1","trkcon2",100,0.,1.);
    _htrkcon2       = tfs->make<TH1F>("_htrkcon2","trkcon2",100,0.,1.);
    _htrkmomErr    = tfs->make<TH1F>("_htrkmomErr","trkmomErr",100,0.,2.);
    _htrkmomErr1    = tfs->make<TH1F>("_htrkmomErr1","trkmomErr1",100,0.,2.);
    _htrkmomErr2    = tfs->make<TH1F>("_htrkmomErr2","trkmomErr2",100,0.,2.);
    _htrkt0           = tfs->make<TH1F>("_htrkt0","trkt0",170,0.,1700.);
    _htrkt0Err        = tfs->make<TH1F>("_htrkt0Err","trkt0Err",100,0.,10.);
    _htrkt0Err1        = tfs->make<TH1F>("_htrkt0Err1","trkt0Err1",100,0.,10.);
    _htrkt0Err2        = tfs->make<TH1F>("_htrkt0Err2","trkt0Err2",100,0.,10.);
    _hfitpar_om    = tfs->make<TH1F>("_hfitpar_om","om",100,0.,0.02);
    _hfitpar_d0    = tfs->make<TH1F>("_hfitpar_d0","d0",100,-200.,500.);
    _hfitpar_z0    = tfs->make<TH1F>("_hfitpar_z0","z0",1000,-1000.,1000.);
    _hfitpar_d0omega    = tfs->make<TH1F>("_hfitpar_d0omega","d0 + 2/om",120,-200.,1000.);
    _hfitpar_td    = tfs->make<TH1F>("_hfitpar_td","tandip",100,0.,3.);
    _hfitmom       = tfs->make<TH1F>("_hfitmom","fitmom",300,0.,300.);
    _hevtwt        = tfs->make<TH1F>("_hevtwt","event weight",300,0.,.035);
    _hfitmomCutSetC = tfs->make<TH1F>("_hfitmomCutSetC","fitmom",300,0.,150.);
    _hfitmomFinal   = tfs->make<TH1F>("_hfitmomFinal","fitmom",10,103.5,105.);
    _hNumberOfTracks = tfs->make<TH1F>("_hNumberOfTracks", "Number of Tracks", 5,0.,5.);
    _hTimeOfCluster = tfs->make<TH1F>("_hTimeOfCluster", "Time Of Cluster", 1700,0.,1700.);
    _hNumberOfClusters          = tfs->make<TH1F>("_hNumberOfClusters","Number Of Clusters", 100,0.,100.);
    _hNumberOfClustersInWindow  = tfs->make<TH1F>("_hNumberOfClustersInWindow","Number Of Clusters In Window", 20,0.,20.);
    _hClusterEnergy             = tfs->make<TH1F>("_hClusterEnergy","Energy Of Clusters", 60,0.,300.);
    _hClusterEnergyInWindow     = tfs->make<TH1F>("_hClusterEnergyInWindow","Energy Of Clusters In Window", 60,0.,300.);
    _hHighestEnergyCluster      = tfs->make<TH1F>("_hHighestEnergyCluster", "Highest Energy Cluster", 60,0.,300.);
    _hMomentumOfTrackForCal     = tfs->make<TH1F>("_hMomentumOfTrackForCal","Momentum of Track Defining Cluster Window", 300,0.,150.);
    _hTotalEnergyInWindow = tfs->make<TH1F>("_hTotalEnergyInWindow","Total Energy In Window", 300,0.,300.);
    _hDistX = tfs->make<TH1F>("_hDistX","Distance from Track to Highest Energy Cluster, X", 100,-1000.,1000.);
    _hDistY = tfs->make<TH1F>("_hDistY","Distance from Track to Highest Energy Cluster, Y", 100,-1000.,1000.);
    _htime = tfs->make<TH1F>("_htime","Time of Cluster - Time of Track at Calorimeter CoGZ",100, -100.,100.);
    _hPhi0 = tfs->make<TH1F>("_hPhi0","Phi0",64, -3.2,3.2);
    _hShortestDistance = tfs->make<TH1F>("_hShortestDistance","Shortest Distance Cluster",100,0.,500.);
    _hShortestFirst    = tfs->make<TH1F>("_hShortestFirst","Shortest Distance First Disk Cluster",100,0.,500.);
    _hShortestDistanceClusterTime = tfs->make<TH1F>("_hShortestDistanceClusterTime","Shortest Distance Cluster Time - Track at Disk",50,-25.,25.);
    _hClusterRadius = tfs->make<TH1F>("_hClusterRadius","Radius of Clusters", 100,0.,1000.);
    _hEnergyAtLargeRadius = tfs->make<TH1F>("_hEnergyAtLargeRadius","Cluster Energy at Large Radius (get rid of DIO and track)",60,0.,300.);
    _hTotalEnergyAtLargeRadius = tfs->make<TH1F>("_hTotalEnergyAtLargeRadius","Summed Energy at Large Radius (get rid of DIO and track)",60,0.,300.);
    _hNumberAtLargeRadius = tfs->make<TH1F>("_hNumberAtLargeRadius","Number at Large Radius (get rid of DIO and track)",20,0.,20.);
    _hPdgIdOfAllClustersAtLargeRadius = tfs->make<TH1F>("_hPdgIdOfAllClustersAtLargeRadius","PDGId of Clusters at Large Radius",5000,-2500.,2500.);

    _hEnergyAtSmallRadius = tfs->make<TH1F>("_hEnergyAtSmallRadius","Cluster Energy at Small Radius (get rid of DIO and track)",60,0.,300.);
    _hTotalEnergyAtSmallRadius = tfs->make<TH1F>("_hTotalEnergyAtSmallRadius","Summed Energy at Small Radius (get rid of DIO and track)",60,0.,300.);
    _hNumberAtSmallRadius = tfs->make<TH1F>("_hNumberAtSmallRadius","Number at Small Radius (get rid of DIO and track)",20,0.,20.);
    _hPdgIdOfAllClustersAtSmallRadius = tfs->make<TH1F>("_hPdgIdOfAllClustersAtSmallRadius","PDGId of Clusters at Small Radius",5000,-2500.,2500.);
    
    _hCluConv = tfs->make<TH1F>("_hCluConv","CluConv", 10,-5.,5.);
    _hEoverP1st = tfs->make<TH1F>("_hEoverP1st","E/p for shortest distance cluster disk 1", 150,0.,1.5);
    _hEoverP2nd = tfs->make<TH1F>("_hEoverP2nd","E/p for shortest distance cluster disk 2", 150,0.,1.5);
    _hEOverPZoom = tfs->make<TH1F>("_hEOverPZoom","E/p for shortest distance cluster", 50,0.8,1.3);
    _hDistToTrack = tfs->make<TH1F>("_hDistToTrack","Distance of Cluster to Track at first disk",100,0.,1000.);

    _hNumberOfHits = tfs->make<TH1F>("_hNumberOfHits","Total Number of Hits",100,0.,5000.);
    _hNumberOfHitsAfter650Nsec = tfs->make<TH1F>("_hNumberOfHitsAfter650Nsec","Total Number of Hits After 650 nsec",100,0.,5000.);
    _hTimeDiffBest1st = tfs->make<TH1F>("_hTimeDiffBest1st","time diff to t0, 1st disk",20, -20.,20.);
    _hTimeDiffBest2nd = tfs->make<TH1F>("_hTimeDiffBest2nd","time diff to t0, 2nd disk",20, -20.,20.);

    _hNHitsInTime = tfs->make<TH1F>("_hNHitsInTime","number of hits in time", 1000,0.,1000.);
    _hZOfHit = tfs->make<TH1F>("_hZOfHit","Z of intime hit", 150,-1500,1500.);
    _hEnergyDeposit = tfs->make<TH1F>("_hEnergyDeposit","Enegy Deposit of Hits",500, 0.,.025);
    _hEnergyDepositFront = tfs->make<TH1F>("_hEnergyDepositFront","Energy Deposit of Hits, Front 10",500, 0.,.025);
    _hNumberOfDeltas = tfs->make<TH1F>("_hNumberOfDeltas","Number of Deltas in Time",1000,0.,1000.);
    _hTimeOfDeltaWrtT0 = tfs->make<TH1F>("_hTimeOfDeltaWrtT0","Time of Delta wrt T0",100,-50.,50.);
    _hTimeOfDeltaWrtT0Zoom = tfs->make<TH1F>("_hTimeOfDeltaWrtT0Zoom","Time of Delta wrt T0 Zoom",100,0.,50.);
    _hCloseDeltas = tfs->make<TH1F>("_hCloseDeltas","CloseInTime Deltas", 50,0.,1.);
    _hCloseDeltasMean = tfs->make<TH1F>("_hCloseDeltasMean","CloseInTime Deltas Mean", 10,0.,50.);
    _hCloseDeltasSigma = tfs->make<TH1F>("_hCloseDeltasSigma","CloseInTime Deltas Sigma", 50,0.,50.);
    _hCloseDeltasAmplitude = tfs->make<TH1F>("_hCloseDeltasAmplitude","CloseInTime Deltas Amplitude", 100,0.,1000.);
    _hMuonSampleMomentum = tfs->make<TH1F>("_hMuonSampleMomentum","Muons with Cut Set C , 8 Active, T>500 nsec",100,0.,200.);
    _hMuonFracCutMomentum = tfs->make<TH1F>("_hMuonFracCutMomentum","Muons with Cut Set C , kFrac>0.7 and < 1 , 8 Active, T>500 nsec",100,0.,200.);
    _hPionFracCutMomentum = tfs->make<TH1F>("_hPionFracCutMomentum","Muons with Cut Set C , kFrac<0.85, 8 Active, T>500 nsec",100,0.,200.);
    _hWeirdFracCutMomentum = tfs->make<TH1F>("_hWeirdFracCutMomentum","Muons with Cut Set C , kFrac>0.85 to 0.9, 8 Active, T>500 nsec",100,0.,200.);
    _hFracVsTime = tfs->make<TH2F>("_hFracVsTime","  E/p vs time diff", 50, -25.,25.,50,0.,1.25);
    _hRadiusVsTime = tfs->make<TH2F>("_hRadiusVsTime","  track to cluster dist vs time diff", 50,0.,250.,50,-25.,25.);
    _hFracVsRadius = tfs->make<TH2F>("_hFracVsRadius","  E/p vs track to cluster dist", 50,0.,250.,50,0.,1.25);
    _hFracVsLocation = tfs->make<TH2F>("_hFracVsLocation","  E/p vs r2 for cluster", 50,0.,1000,50,0.,1.25);
    _hMuonFracVsLocation = tfs->make<TH2F>("_hMuonFracVsLocation","  MuonFrac vs r2 for cluster", 50,0.,1000,50,0.75,1.25);

    _hFracVsP = tfs->make<TH2F>("_hFracVsP","E/p vs p",30,60.,120.,100,0.,1.);
    _hMuonFrac =  tfs->make<TH1F>("_hMuonFrac","ClusterEnergy/KE assuming Muon",100,0.5,1.5);
    _hPionFrac =  tfs->make<TH1F>("_hPionFrac","ClusterEnergy/KE assuming Pion, Muon Frac lt 0.85",100,0.5,1.5);
    _hWeirdFracCluster = tfs->make<TH1F>("_hWeirdFracCluster","Cluster Energy for frac 0.85 to 0.9",100,0.,100.);
    _hMuonFracCluster = tfs->make<TH1F>("_hMuonFracCluster","Cluster Energy for frac gt 0.85",100,0.,100.);
    _hPionFracCluster = tfs->make<TH1F>("_hPionFracCluster","Cluster Energy for frac lt 0.85",100,0.,100.);
    _hNumberOfCrystalsHighFrac = tfs->make<TH1F>("_hNumberOfCrystalsHighFrac","Number of Crystals, Frac 0.96 to 0.98",10,0.,10.);
    _hNumberOfCrystalsLowFrac = tfs->make<TH1F>("_hNumberOfCrystalsLowFrac","Number of Crystals, Frac 0.86 to 0.88",10,0.,10.);
    _hNumberOfCrystalsMedFrac = tfs->make<TH1F>("_hNumberOfCrystalsMedFrac","Number of Crystals, Frac gt 0.85 and < 0.9",10,0.,10.);
    _hSingleCrystalFrac = tfs->make<TH1F>("_hSingleCrystalFrac","Frac for Clusters With Single Crystals",100,0.5,1.5);
    _hDoubleCrystalFrac = tfs->make<TH1F>("_hDoubleCrystalFrac","Frac for Clusters With Two Crystals",100,0.5,1.5);
    _hTripleCrystalFrac = tfs->make<TH1F>("_hTripleCrystalFrac","Frac for Clusters With Three Crystals",100,0.5,1.5);
    _hCascade = tfs->make<TH1F>("_hCascade","Building up Track Cuts",10, 0.,10.);
    // TF1* func = new TF1("timeFit",timeFit,0.,50.,3);

    _hErrVsMomHiFrac = tfs->make<TProfile>("_hErrVsMomHiFrac","Err vs Momentum HiFrac Profile",300,0.,300.);
    _hErrVsD0HiFrac = tfs->make<TProfile>("_hErrVsD0HiFrac","Err vs D0 HiFrac Profile",300,-100,500.);
    _hErrVsOmegaHiFrac = tfs->make<TProfile>("_hErrVsOmegaHiFrac","Err vs Omega HiFrac Profile",100,0.,.01);
    _hErrVsTanDipHiFrac = tfs->make<TProfile>("_hErrVsTanDipHiFrac","Err vs TanDip HiFrac Profile",300,0.,3.);
    _hTanDipVsMomentumHiFrac = tfs->make<TH2F>("_hTanDipVsMomentumHiFrac","Tan Dip vs Momentum HiFrac",300,0.,300.,100,0.,3.);

    _hErrVsMomLoFrac = tfs->make<TProfile>("_hErrVsMomLoFrac","Err vs Momentum LoFrac Profile",300,0.,300.);
    _hErrVsD0LoFrac = tfs->make<TProfile>("_hErrVsD0LoFrac","Err vs D0 LoFrac Profile",300,-100,500.);
    _hErrVsOmegaLoFrac = tfs->make<TProfile>("_hErrVsOmegaLoFrac","Err vs Omega LoFrac Profile",100,0.,.01);
    _hErrVsTanDipLoFrac = tfs->make<TProfile>("_hErrVsTanDipLoFrac","Err vs TanDip LoFrac Profile",300,0.,3.);
    _hTanDipVsMomentumLoFrac = new TH2F("_hTanDipVsMomentumLoFrac","Tan Dip vs Momentum LoFrac",300,0.,300.,100,0.,3.);


    _hErrVsMomMuon = tfs->make<TProfile>("_hErrVsMomMuon","Err vs Momentum Muon Profile",300,0.,300.);
    _hErrVsD0Muon = tfs->make<TProfile>("_hErrVsD0Muon","Err vs D0 Muon Profile",300,-100,500.);
    _hErrVsOmegaMuon = tfs->make<TProfile>("_hErrVsOmegaMuon","Err vs Omega Muon Profile",100,0.,.01);
    _hErrVsTanDipMuon = tfs->make<TProfile>("_hErrVsTanDipMuon","Err vs TanDip Muon Profile",300,0.,3.);
    _hTanDipVsMomentumMuon = new TH2F("_hTanDipVsMomentumMuon","Tan Dip vs Momentum Muon",300,0.,300.,100,0.,3.);



    _hErrVsMomPion = tfs->make<TProfile>("_hErrVsMomPion","Err vs Momentum Pion Profile",300,0.,300.);
    _hErrVsD0Pion = tfs->make<TProfile>("_hErrVsD0Pion","Err vs D0 Pion Profile",300,-100,500.);
    _hErrVsOmegaPion = tfs->make<TProfile>("_hErrVsOmegaPion","Err vs Omega Pion Profile",100,0.,.01);
    _hErrVsTanDipPion = tfs->make<TProfile>("_hErrVsTanDipPion","Err vs TanDip Pion Profile",300,0.,3.);
    _hTanDipVsMomentumPion = tfs->make<TH2F>("_hTanDipVsMomentumPion","Tan Dip vs Momentum Pion",300,0.,300.,100,0.,3.);


    //set up some constants


  }



  void KineticFracAnalysis::endJob(){

  }




  void KineticFracAnalysis::analyze(const art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from KineticFracAnalysis =  "<<_nProcess << " with instance name " << _instanceName <<std::endl;

    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);

    ConditionsHandle<CalorimeterPhysicalConstants> calPhys("ignored");
    _density = calPhys->density();
    //    std::cout << "density is " << _density/(CLHEP::g/CLHEP::cm3) << std::endl;
    _eCritPos = calPhys->criticalEnergyPos();
    _eCritNeg = calPhys->criticalEnergyNeg();

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    _firstDiskZ  = cal.toTrackerFrame(cal.section(0).frontFaceCenter()).z();
    _secondDiskZ = cal.toTrackerFrame(cal.section(1).frontFaceCenter()).z();

    //    std::cout << "checking disk locations " << _firstDiskZ << " " <<  _secondDiskZ << std::endl;

    //    std::cout << " and other end " << cal.toTrackerFrame(cal.section(0).frontFaceCenter()) << " " << cal.toTrackerFrame(cal.section(1).frontFaceCenter())<<std::endl;
    //Get handle to the tracker
    if( ! geom->hasElement<TTracker>() ) return;
    //      TTracker const & tracker = *(GeomHandle<TTracker>());
       art::Handle<TrkCaloMatchCollection>  trkCaloMatchHandle;
        event.getByLabel(_trkCaloMatchModuleLabel, trkCaloMatchHandle);
        TrkCaloMatchCollection const& trkCaloMatches(*trkCaloMatchHandle);

        art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle;
        event.getByLabel(_trkIntersectModuleLabel, trjIntersectHandle);
        TrkCaloIntersectCollection const& trkIntersect(*trjIntersectHandle);

    //Get generated particles
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);
    GenParticleCollection const& genParticles(*gensHandle);

    //Get calorimeter readout hits (2 readout / crystal as of today)
    art::Handle<CaloHitCollection> caloHitsHandle;
    event.getByLabel(_caloReadoutModuleLabel, caloHitsHandle);
    CaloHitCollection const& caloHits(*caloHitsHandle);

    //Get calorimeter readout hits MC level - energy/time/type
    art::Handle<CaloHitMCTruthCollection> caloHitMCTruthHandle;
    event.getByLabel(_caloReadoutModuleLabel, caloHitMCTruthHandle);
    CaloHitMCTruthCollection const& caloHitsMCTruth(*caloHitMCTruthHandle);

    //Get simParticles and stepPointMC summary for crystal readout hits
    art::Handle<CaloHitSimPartMCCollection> caloHitSimMCHandle;
    event.getByLabel(_caloReadoutModuleLabel, caloHitSimMCHandle);
    CaloHitSimPartMCCollection const& caloHitSimPartMC(*caloHitSimMCHandle);

    //Get calo crystal hits (average from readouts)
    art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
    event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
    CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);

    //Get calo cluster
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    CaloClusterCollection const& caloClusters(*caloClustersHandle);

    //Get virtual detector hits
    art::Handle<StepPointMCCollection> vdhits;
    event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);

    //Get tracks
    art::Handle<KalRepCollection> krepsHandle;
    event.getByLabel(_trkPatRecModuleLabel,_instanceName,krepsHandle);
    KalRepCollection const& kreps = *krepsHandle;


    //
    // handle to PDG
    GlobalConstantsHandle<ParticleDataTable> pdt;
    ParticleDataTable const & pdt_ = *pdt;
    
    //Utility to match  cloHits with MCtruth, simParticles and StepPoints
    CaloHitMCNavigator caloHitNavigator(caloHits, caloHitsMCTruth, caloHitSimPartMC);


    const double CrDensity = 4.9*(CLHEP::g/CLHEP::cm3);
    const double CrMass    = CrDensity*cal.caloGeomInfo().crystalVolume();

    double numberOfTracks = 0;
    double numberOfClusters = 0;




    //--------------------------  Do generated particles --------------------------------


    _evt = event.id().event();
    _run = event.run();

    if (_diagLevel == 3){std::cout << "processing event in Kinetic_Frac " << _nProcess << " run and event  = " << _run << " " << _evt << " with instance name = " << _instanceName << std::endl;}


    _nGen = genParticles.size();
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
      }



    //--------------------------  Do calorimeter hits --------------------------------

    _nHits = _nSim = 0;
    _cryEtot = 0.0;

    for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
      {
	CaloCrystalHit const& hit    = caloCrystalHits.at(ic);
	int sectionId                  = cal.crystal(hit.id()).sectionId();
	CLHEP::Hep3Vector crystalPos   = cal.crystal(hit.id()).localPositionFF();  //in disk FF frame
	CaloHit const& caloHit       = *(hit.readouts().at(0));


	CaloHitSimPartMC const& hitSim = caloHitNavigator.sim(caloHit);
	int nPartInside                = hitSim.simParticles().size();


	_hcryE->Fill(hit.energyDep());
	_hcryT->Fill(hit.time());
	_hcryX->Fill(crystalPos.x());
	_hcryY->Fill(crystalPos.y());
	_hcryZ->Fill(crystalPos.z());



	_cryEtot             += hit.energyDep();
	_cryTime[_nHits]      = hit.time();
	_cryEdep[_nHits]      = hit.energyDep();
	_cryDose[_nHits]      = hit.energyDep() / CrMass / (CLHEP::joule/CLHEP::kg); //dose
	_cryPosX[_nHits]      = crystalPos.x();
	_cryPosY[_nHits]      = crystalPos.y();
	_cryPosZ[_nHits]      = crystalPos.z();
	//	   std::cout << "z Position of crystal " << crystalPos.z() << std::endl;


	_cryId[_nHits]        = hit.id();
	_crySectionId[_nHits] = sectionId;
	_crySimIdx[_nHits]    = _nSim;
	_crySimLen[_nHits]    = nPartInside;


	for (int ip=0; ip<nPartInside;++ip)
	  {

	    art::Ptr<SimParticle> const& mother = hitSim.simParticles().at(ip);

	    art::Ptr<SimParticle> grandMother = mother;
	    while (grandMother->hasParent()) grandMother = grandMother->parent();
	    GenParticle const* generated = grandMother->genParticle() ? grandMother->genParticle().get() : 0;
	    CLHEP::Hep3Vector hitSimPos = cal.toSectionFrameFF(sectionId,hitSim.position().at(ip)); //in disk FF frame

	    _motId[_nSim]      = mother->id().asInt();
	    _motPdgId[_nSim]   = mother->pdgId();
	    _motmom[_nSim]     = hitSim.momentum().at(ip);
	    _motcrCode[_nSim]  = mother->creationCode();
	    _motStartX[_nSim]  = mother->startPosition().x(); //in Mu2e frame
	    _motStartY[_nSim]  = mother->startPosition().y();
	    _motStartZ[_nSim]  = mother->startPosition().z();
	    _motStartT[_nSim]  = mother->startGlobalTime();
	    _motPosX[_nSim]    = hitSimPos.x(); // in disk FF frame
	    _motPosY[_nSim]    = hitSimPos.y();
	    _motPosZ[_nSim]    = hitSimPos.z();
	    _motTime[_nSim]    = hitSim.time().at(ip);
	    _motEdep[_nSim]    = hitSim.eDep().at(ip);
	     
	    _motGenIdx[_nSim]  = -1;
	    if (generated) _motGenIdx[_nSim] = generated - &(genParticles.at(0));
	    ++_nSim;

	  }

	++_nHits;
      }



    //--------------------------  Do clusters --------------------------------

    _nCluster = _nCluSim = 0;
    _cluList.clear();
    for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
      {

	CaloContentMC clutil(caloHitNavigator,*clusterIt);
	std::vector<art::Ptr<SimParticle> > const& sim1 = clutil.simPart();

	std::vector<int> _list;
	for (int i=0;i<clusterIt->size();++i)
          {
            int idx = int(clusterIt->caloCrystalHitsPtrVector().at(i).get()- &caloCrystalHits.at(0));
            _list.push_back(idx);
          }
	_cluEnergy[_nCluster] = clusterIt->energyDep();
	_cluTime[_nCluster]   = clusterIt->time();
	_cluNcrys[_nCluster]  = clusterIt->size();
	_cluCogX[_nCluster]   = clusterIt->cog3Vector().x();
	_cluCogY[_nCluster]   = clusterIt->cog3Vector().y();
	_cluCogZ[_nCluster]   = clusterIt->cog3Vector().z();
	if (clusterIt->sectionId() == 0)
	  {_cluCogZ[_nCluster]   = _firstDiskZ;} 
	else 
	  {_cluCogZ[_nCluster]   = _secondDiskZ;}

	_cluConv[_nCluster]   = clutil.hasConversion();
	_cluSimIdx[_nCluster] = _nCluSim;
	_cluSimLen[_nCluster] = sim1.size();
	_cluList.push_back(_list);

	_hcluE->Fill(clusterIt->energyDep());
	_hcluT->Fill(clusterIt->time());
	_hcluX->Fill(clusterIt->cog3Vector().x());
	_hcluY->Fill(clusterIt->cog3Vector().y());
	_hcluZ->Fill(clusterIt->cog3Vector().z());
	_hcluE1Et->Fill(clusterIt->e1()/clusterIt->energyDep());
	_hcluE1E9->Fill(clusterIt->e1()/clusterIt->e9());
	_hcluE1E25->Fill(clusterIt->e1()/clusterIt->e25());
	//note: a signal cluster is defined as _cluConv==1 and _cluSimLen==1 (if it has a signal and there is only one inside)

	for (unsigned int ip=0; ip<sim1.size();++ip)
          {
	    art::Ptr<SimParticle> smother = sim1[ip];
	    while (smother->hasParent()) smother = smother->parent();
	    int genIdx=-1;
	    if (smother->genParticle()) genIdx = smother->genParticle()->generatorId().id();
	    CLHEP::Hep3Vector cluSimPos = cal.toSectionFrameFF(clusterIt->sectionId(),clutil.position().at(ip));

	    _clusimId[_nCluSim]     = sim1[ip]->id().asInt();
	    _clusimPdgId[_nCluSim]  = sim1[ip]->pdgId();
	    _clusimGenIdx[_nCluSim] = genIdx;
	    _clusimMom[_nCluSim]    = clutil.momentum().at(ip);
	    _clusimPosX[_nCluSim]   = cluSimPos.x(); // in disk FF frame
	    _clusimPosY[_nCluSim]   = cluSimPos.y();
	    _clusimPosZ[_nCluSim]   = cluSimPos.z();  
	    _clusimTime[_nCluSim]   = clutil.time().at(ip);
	    _clusimEdep[_nCluSim]   = clutil.edepTot().at(ip);
	     
	    ++_nCluSim;
          }

	++_nCluster;
      }
    numberOfClusters = _nCluster;


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

	    double hitTimeUnfolded = _toff.timeWithOffsetsApplied(hit);
	    double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	    CLHEP::Hep3Vector VDPos = cal.toTrackerFrame(hit.position());

	    _vdId[_nVd]    = hit.volumeId();
	    _vdPdgId[_nVd] = hit.simParticle()->pdgId();
	    _vdTime[_nVd]  = hitTime;//hit.time();
	    _vdPosX[_nVd]  = VDPos.x(); //tracker frame
	    _vdPosY[_nVd]  = VDPos.y();
	    _vdPosZ[_nVd]  = VDPos.z();

	    _vdMom[_nVd]   = hit.momentum().mag();
	    _vdenIdx[_nVd] = hit.simParticle()->generatorIndex();
	    ++_nVd;
	  }
      }

        _nMatch=0;
        for (auto const& trkCaloMatch: trkCaloMatches)
        {
              _mTrkId[_nMatch]     = trkCaloMatch.trkId();
              _mCluId[_nMatch]     = trkCaloMatch.cluId();
              _mChi2[_nMatch]      = trkCaloMatch.chi2();
              _mChi2Pos[_nMatch]   = trkCaloMatch.chi2Pos();
              _mChi2Time[_nMatch]  = trkCaloMatch.chi2Time();

              ++_nMatch;
        }

        _nTrk=0;
        for (auto const& intersect: trkIntersect)
        {
              KalRep const& trk = *(intersect.trk());

              double pathLength             = intersect.pathLengthEntrance();
              CLHEP::Hep3Vector trkMomentum = trk.momentum(pathLength);
              double trkTime                = trk.arrivalTime(pathLength);
              double tprob                  = trk.chisqConsistency().significanceLevel();
              int    tNhits                 = trk.nActive();
              int    tStatus                = trk.fitStatus().success();

	      std::cout << "time, prob, nhits, status = " << trkTime << " " << tprob << " " << " " << tNhits << " " << tStatus << std::endl;

              HepPoint point = trk.position(pathLength);
              CLHEP::Hep3Vector posTrkInTracker(point.x(),point.y(),point.z());
              CLHEP::Hep3Vector posTrkInSectionFF = cal.toSectionFrameFF(intersect.sectionId(),cal.fromTrackerFrame(posTrkInTracker));

              HelixTraj trkHel(trk.helix(pathLength).params(),trk.helix(pathLength).covariance());

	      std::cout << "trk.d0 = " << trkHel.d0() << std::endl;
           }



    //--------------------------  Do tracks  --------------------------------
    _nTrkOk = 0;
    _nTrk = 0;

    for ( KalRepCollection::const_iterator i=kreps.begin(); i != kreps.end(); ++i ){

      KalRep const& krep = **i;

      // Arc length from center of tracker to the most upstream point on the fitted track.
      double s0 = krep.startValidRange();

      // Momentum and position at s0.
      CLHEP::Hep3Vector p0     = krep.momentum(s0);
      HepPoint          pos0   = krep.position(s0);

      // Some other properties at s0.
      double loclen(0.);
      const TrkSimpTraj* ltraj = krep.localTrajectory(s0,loclen);
      HepVector momvec(3);
      momvec = p0.unit();
      BbrVectorErr momCov = krep.momentumErr(s0);
      double fitMomErr    = sqrt(momCov.covMatrix().similarity(momvec));
      double tanDip       = ltraj->parameters()->parameter()[HelixTraj::tanDipIndex];
      double omega        = ltraj->parameters()->parameter()[HelixTraj::omegaIndex];
      double d0           = ltraj->parameters()->parameter()[HelixTraj::d0Index];
      double z0           = ltraj->parameters()->parameter()[HelixTraj::z0Index];
      double trkPhi0      = ltraj->parameters()->parameter()[HelixTraj::phi0Index];
      double fitCon       = krep.chisqConsistency().significanceLevel();
      double trkt0Err     = krep.t0().t0Err();
      double fitmompt = p0.mag()*(1.0-p0.cosTheta()*p0.cosTheta());

      // Does this fit pass cut set C?
      bool cutC = ( krep.fitStatus().success() >0) &&
	( fitCon             > 2.e-3  ) &&
	( krep.nActive()    >= 20    ) &&
	( fitMomErr          < 0.25   ) &&
	( krep.t0().t0() > 700 && krep.t0().t0() < 1695);

      _trkDip[_nTrk] = tanDip;
      _trkOk[_nTrk]  = cutC? 1 : 0;
      _trkpt[_nTrk]  = fitmompt;
      _trkstat[_nTrk]  = krep.fitStatus().success();
      _trkcon[_nTrk]  = fitCon;
      _trknHit[_nTrk]  = krep.nActive();
      _trkmomErr[_nTrk]  = fitMomErr ;
      _trkt0[_nTrk] = krep.t0().t0();
      _trkMom[_nTrk] = p0.mag();
      _trkd0[_nTrk] = d0;
      _trkz0[_nTrk] = z0;
      _trkOmega[_nTrk] = omega;
      _trkPhi0[_nTrk] = trkPhi0;
      _trkt0Err[_nTrk] = trkt0Err;
      ++_nTrk;

      //              std:: cout << "KineticFracAnalysis t0 = " << krep.t0().t0() << std::endl;
      //      std:: cout << "KineticFracAnalysis fitstatus = " << krep.fitStatus().success() << " " << _fitStatus[_nTrk] << std::endl;

      if (cutC)  ++_nTrkOk;
    }
    numberOfTracks = _nTrk;
    //--------------------------  Do tracker hits  --------------------------------

    //   const Tracker& tracker = getTrackerOrThrow();

    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if (event.getByLabel(_shLabel,strawhitsH)) {
      _shcol = strawhitsH.product();
    }
    else {
      _shcol  = 0;
      printf(" >>> ERROR in KineticFracAnalysis::findData: StrawHitCollection with label=%s not found.\n",
	     _shLabel.data());
    }

    art::Handle<mu2e::StrawHitPositionCollection> shposH;
    if (event.getByLabel(_shpLabel,shposH)) {
      _shpcol = shposH.product();
    }
    else {
      _shpcol = 0;
      printf(" >>> ERROR in KineticFracAnalysis::findData: StrawHitPositionCollection with label=%s not found.\n",
	     _shpLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if (event.getByLabel(_shfLabel,shflagH)) {
      _shfcol = shflagH.product();
    }
    else {
      _shfcol = 0;
      printf(" >>> ERROR in KineticFracAnalysis::findData: StrawHitFlagCollection with label=%s not found.\n",
	     _shfLabel.data());
    }

    art::Handle<mu2e::StrawHitFlagCollection> bkflagH;
    if (event.getByLabel(_bkfLabel,bkflagH)) {
      _bkfcol = bkflagH.product();
    }
    else {
      _bkfcol = 0;
      printf(" >>> ERROR in KineticFracAnalysis::findData: StrawHitFlagCollection with label=%s not found.\n",
	     _bkfLabel.data());
    }


    _nTrackerHits = _shcol->size();
    _nGoodHits = 0;
    _nDeltas = 0;

    //   std::cout << "size of flag collection = " << _shfcol->size() << std::endl;
    //std::cout << "_ntrackerHits = " << _nTrackerHits << std::endl;


    for (int ithHit =0; ithHit<_nTrackerHits; ++ithHit)
      {
	//	std::cout << "inside hit loop " << _shcol->at(ithHit).time() << std::endl;

	if (_shcol->at(ithHit).time() > 650. && _nGoodHits< 32767)
	  {
	    _hitTime[_nGoodHits] = _shcol->at(ithHit).time();
	    _hitEnergy[_nGoodHits] = _shcol->at(ithHit).energyDep();
	    _hitZ[_nGoodHits] = _shpcol->at(ithHit).pos().z();
	    bool deltaRay = _bkfcol->at(ithHit).hasAllProperties(StrawHitFlag::delta);
	    //	    	    bool deltaRay = _shpcol->at(ithHit).flag().hasAllProperties(StrawHitFlag::delta);
	    if (deltaRay){
	      _hitIsDelta[_nGoodHits] = 1;
	      ++_nDeltas;
	    } else{_hitIsDelta[_nGoodHits] = 0;}
	    //if (deltaRay){std::cout << "found a delta" << std::endl;}
	    //std::cout << "time, energy,z " << _shcol->at(ithHit).time() << " " << _shcol->at(ithHit).energyDep() << " " <<
	    //		      _shpcol->at(ithHit).pos().z() 
	    //           << std::endl;
	    ++_nGoodHits;
	  }
      }

    //std::cout << " number of good hits = " << _nGoodHits << std::endl;

    //    _Ntup->Fill();

    // extrapolate the track, see if lines up with a cluster, and calculate the kinetic fraction.

    // set validity variable.
    //    std::cout << "from KineticFrac, upstream disk edge " << cal.section(0).zUpInTracker() << " " << cal.section(1).zUpInTracker() << std::endl;
    //  std::cout << "                 and downstream ends " << cal.section(0).zDownInTracker() << " " << cal.section(1).zDownInTracker() << std::endl;
    _goodKfrac = 0;

    if (numberOfTracks   == 0){_goodKfrac =1;}
    if (numberOfClusters == 0 && numberOfTracks > 0){_goodKfrac = 2;}
    if (numberOfClusters > 0 && numberOfTracks == 0){_goodKfrac = 3;}
    int signFlip = 1;
    //    std::cout << "number of tracks = " << numberOfTracks << std::endl;
    for (int ithtrk =  0; ithtrk < numberOfTracks;    ++ithtrk) {
      if (_trkstat[ithtrk] > 0){

	//	    std::cout << "inside event id, fitstatus " << entry << " " << trkstat[ithtrk] << std::endl;
	_hnactive->Fill(_trknHit[ithtrk]);
	//if (nActiveHits < 15){std::cout << "nActiveHits = " << nActiveHits << std::endl;}
	_htrkcon->Fill(_trkcon[ithtrk]);
	_htrkmomErr->Fill(_trkmomErr[ithtrk]);
	_htrkt0->Fill(_trkt0[ithtrk]);
	_htrkt0Err->Fill(_trkt0Err[ithtrk]);
	_hfitpar_d0->Fill(signFlip*_trkd0[ithtrk]);
	_hfitpar_z0->Fill(_trkz0[ithtrk]);
	_hfitpar_d0omega->Fill(signFlip*(_trkd0[ithtrk]+ 2./_trkOmega[ithtrk]));
	_hfitpar_om->Fill(signFlip*_trkOmega[ithtrk]);
	_hfitpar_td->Fill(_trkDip[ithtrk]);
	_hfitmom->Fill(_trkMom[ithtrk]);
	_hPhi0->Fill(_trkPhi0[ithtrk]);
	/*
	  _hErrVsMom->Fill(trkMom[ithtrk],trkmomErr[ithtrk]);
	  _hErrVsD0->Fill(signFlip*trkd0[ithtrk],trkmomErr[ithtrk]);
	  _hErrVsOmega->Fill(signFlip*trkOmega[ithtrk],trkmomErr[ithtrk]);
	  _hErrVsTanDip->Fill(trkDip[ithtrk],trkmomErr[ithtrk]);
	  _hTanDipVsMomentum->Fill(trkMom[ithtrk],trkDip[ithtrk]);
	*/

	//
	// coordinate change from center of tracker to Mu2e
	//_trkz0[ithtrk] = _trkz0[ithtrk] + 10175.;

	//
	//known offset from pezzullo, see his producer file
		_trkt0[ithtrk] -= 1.4;
		_trkt0[ithtrk] -= 1.4;// do twice for a minute while I debug
		
      }	
    }
    
    for (int ithtrk =  0; ithtrk < numberOfTracks;    ++ithtrk) {

      // first thing to do is extrapolate track to the disks.  Does it get through first disk and end up striking second disk?

      if (!_trkOk[ithtrk]) break;
      
      //
      // extrapolate to first disk.  just use the helix since the spread on the 
      //tracker is small compared to the location from the calorimeter.  On the 
      //other hand, knowing from the track where the shower hit might be helpful 
      //in making energy corrections for energy going down cracks.  
      //So eventually will have to look at both

      std::cout << "got a track" << ithtrk << std::endl;

      double betaZ = _trkDip[ithtrk]/sqrtOrThrow(1. + _trkDip[ithtrk]*_trkDip[ithtrk],0.00001);

      //
      // extapolate track to first disk.  If it is in the active area, say this is a first disk track; if it goes in the hole 
      // and extrapolates in the active area of the second disk, call it a second disk track. If neither, don't bother with it.
      double flightLength = (_firstDiskZ - _trkz0[ithtrk])/_trkDip[ithtrk];
      double trkX =  sin(_trkPhi0[ithtrk] + _trkOmega[ithtrk]*flightLength)/_trkOmega[ithtrk]  
	- (1./_trkOmega[ithtrk] + _trkd0[ithtrk])*sin(_trkPhi0[ithtrk] );
      double trkY = -cos(_trkPhi0[ithtrk] + _trkOmega[ithtrk]*flightLength)/_trkOmega[ithtrk]  
	+ (1./_trkOmega[ithtrk] + _trkd0[ithtrk])*cos(_trkPhi0[ithtrk] );

      double flightLengthSecond = (_secondDiskZ  - _trkz0[ithtrk])/_trkDip[ithtrk];
      double trkX2nd =  sin(_trkPhi0[ithtrk] + _trkOmega[ithtrk]*flightLengthSecond)/_trkOmega[ithtrk]  
	- (1./_trkOmega[ithtrk] + _trkd0[ithtrk])*sin(_trkPhi0[ithtrk] );
      double trkY2nd = -cos(_trkPhi0[ithtrk] + _trkOmega[ithtrk]*flightLengthSecond)/_trkOmega[ithtrk]  
	+ (1./_trkOmega[ithtrk] + _trkd0[ithtrk])*cos(_trkPhi0[ithtrk] );

      double radiusAtFirstDisk =  sqrt(trkX*trkX + trkY*trkY);
      double radiusAtSecondDisk = sqrt(trkX2nd*trkX2nd + trkY2nd*trkY2nd);
      /*
      HelixTraj trkHel1(ltraj->helix(flightLength) , ltraj->helix(flightLengt).covariance());
      HelixTraj trkHel2(ltraj->helix(flightLengthSecond) , ltraj->helix(flightLengthSecond).covariance());

      double phi01         = trkHel1.phi0();
      double omega1         = trkHel1.omega();
      double radius1        = 1.0/omega;
      double cosDip1        = trkHel1.cosDip();
      double sinDip1        = sqrt(1-cosDip*cosDip);
      double centerCircleX1 = (trkHel1.d0() + radius)*sin(phi0);
      double centerCircleY1 = (trkHel1.d0() + radius)*cos(phi0);	     
      double trkX1Full1 =  radius1*sin(phi01+omega1*cosDip*length) - centerCircleX1;
      double trkY1Full1 = -radius1*cos(phi01+omega1*cosDip*length) + centerCircleY1;

      
      double phi02         = trkHel2.phi0();
      double omega2         = trkHel2.omega();
      double radius2        = 1.0/omega;
      double cosDip2        = trkHel2.cosDip();
      double sinDip2        = sqrt(1-cosDip*cosDip);
      double centerCircleX2 = (trkHel2.d0() + radius)*sin(phi0);
      double centerCircleY2 = (trkHel2.d0() + radius)*cos(phi0);	     
      double trkXFull2 =  radius2*sin(phi02+omega2*cosDip*length) - centerCircleX2;
      double trkYFull2 = -radius2*cos(phi02+omega2*cosDip*length) + centerCircleY2;

      std::cout << "full and partial = " << trkx1 << " " << trkXFull1 << " " << trkYFull2 << std::endl;
      */
      std::cout << " radius at first disk =  " << radiusAtFirstDisk << " " << _firstDiskZ << std::endl;
      std::cout << " radius at second disk = " << radiusAtSecondDisk << " " << _secondDiskZ << std::endl;

      //	    std::cout << " a few checks " << zClus << " "  << radiusAtFirstDisk << " " << _firstDiskZ << std::endl;
	    //arbitrarily put 1 cm safety zone in place: which disk?


	    bool hitFirstDisk  = false; 
	    bool hitSecondDisk = false;

	    double safetyOffset = 0.;
	    //	    std::cout << "inner, outer radii = " << cal.section(0).innerEnvelopeR() << " " << cal.section(1).innerEnvelopeR() << " " << cal.section(0).outerEnvelopeR() << " " << cal.section(0).outerEnvelopeR() << std::endl;
	    if (radiusAtFirstDisk > cal.section(0).innerEnvelopeR() + safetyOffset && radiusAtFirstDisk < cal.section(0).outerEnvelopeR())
	      {
	      hitFirstDisk = true;
	      hitSecondDisk = false;
	      //	      std::cout << "hit first disk " << hitFirstDisk << " " << hitSecondDisk << std::endl;
	      }

	    if (radiusAtFirstDisk < cal.section(0).innerEnvelopeR() - safetyOffset && !(radiusAtFirstDisk > cal.section(0).outerEnvelopeR())
		&&
		(radiusAtSecondDisk > cal.section(1).innerEnvelopeR() + safetyOffset && radiusAtSecondDisk < cal.section(1).outerEnvelopeR()))
	      {
		hitFirstDisk = false;
		hitSecondDisk = true;
		//
		// replace by 2nd disk extrapolation
		std::cout << "hit second disk " << hitFirstDisk << " " << radiusAtFirstDisk << " " << radiusAtSecondDisk << " " << hitSecondDisk << std::endl;
		std::cout << "trkx,y" << " " << trkX2nd << " " << trkY2nd << " " << flightLengthSecond << std::endl;
		trkX = trkX2nd;
		trkY = trkY2nd;
		flightLength = flightLengthSecond;
	      }
	    //
	    //if you didn't hit either disk stop.

	    std::cout << "trkx, trky before at first disk = " << trkX << " " << trkY << std::endl;
	    if (!hitFirstDisk && !hitSecondDisk) break;

	    int ithClusMatch = -1;
	    bool foundMatch = false;
	    double fracBest = -1;
	    double timeBest = -100.;
	    double minDist = cal.section(1).outerEnvelopeR()*cal.section(1).outerEnvelopeR(); // big number
	    
      
	    for (int ithclus = 0; ithclus < numberOfClusters; ++ithclus){

	      //	std::cout << "got a cluster" << std::endl;
	      //if this particle is being fitted under an e+/e- hypothesis, want to make a correction;
	      //if under a muon hypothesis, correction is zero
	      Float_t showerCorr;
	      if (_tpart == TrkParticle::e_minus || _tpart == TrkParticle::e_plus)
		{ 

		  //	   std::cout << "in here" << std::endl;
		  //std::cout << "critical energy is " << _eCritPos << " or " << _eCritPos*CLHEP::MeV << std::endl;
		  //std::cout << "cluster energy is " << _cluEnergy[ithclus] << " or " << _cluEnergy[ithclus]*CLHEP::MeV << std::endl;
		  //std::cout << "some checks" << " " << _cluEnergy[ithclus]/11. << std::endl;
		  //std::cout << "some more checks" << " " << _cluEnergy[ithclus]/_eCritPos << std::endl;
		  //std::cout << "even more checks" << " " << _cluEnergy[ithclus]*CLHEP::MeV/(_eCritPos*CLHEP::MeV) << std::endl;
		  //std::cout << "radiation length = " << calPhys->radiationLength() << std::endl;
		  //std::cout << "and the correction to shower location is: " << showerCorr << std::endl;
		  //
		  //adjusts for shower development, formula out of Rossi
		  if (_tpart == TrkParticle::e_plus)
		    {
		      showerCorr = (  log(_cluEnergy[ithclus]*CLHEP::MeV/_eCritPos*CLHEP::MeV) - 0.1 )*calPhys->radiationLength()*CLHEP::mm;
		    }
		  else if (_tpart == TrkParticle::e_minus)
		    {
		      showerCorr = (  log(_cluEnergy[ithclus]*CLHEP::MeV/_eCritNeg*CLHEP::MeV) - 0.1 )*calPhys->radiationLength()*CLHEP::mm;
		    }
		  else {showerCorr = 0.;}
		  // showerCorr along axis; how far to project in z?
		  showerCorr *= betaZ;
		  // std::cout << "shower and showerCorr = " << _cluCogZ[ithclus] << " " << showerCorr << std::endl;
		  double zClus = _cluCogZ[ithclus] + showerCorr ; //need to add code for which disk this is, use right coordinates, etc.  placeholder
		  //
		  // recalculate x,y with shower correction
		  trkX =  sin(_trkPhi0[ithtrk] + _trkOmega[ithtrk]*(flightLength+showerCorr))/_trkOmega[ithtrk]  
		    - (1./_trkOmega[ithtrk] + _trkd0[ithtrk])*sin(_trkPhi0[ithtrk] );
		  trkY = -cos(_trkPhi0[ithtrk] + _trkOmega[ithtrk]*(flightLength+showerCorr))/_trkOmega[ithtrk]  
		    + (1./_trkOmega[ithtrk] + _trkd0[ithtrk])*cos(_trkPhi0[ithtrk] );

		  //std::cout << "trkx, y after shower " <<  trkX << " " << trkY << std::endl;
		  //
		  // check how close this cluster is to the track, but make sure it is (a) close enough in time and (b) in the correct disk!
		  double flightZ = (zClus - _trkz0[0]); 

		  double flightTime = 0;

		  if (_tpart == TrkParticle::e_plus ||_tpart == TrkParticle::e_minus)
		    {
		      double betaSpeed = _trkMom[ithtrk]/(sqrt(_trkMom[ithtrk]*_trkMom[ithtrk] + pow(pdt_.particle(TrkParticle::e_minus).ref().mass(),2)));
		      flightTime = flightZ/(betaZ*3.0e11)/1.0e-9/betaSpeed; //in nsec
		      //		      std::cout << "flight time = " << flightTime << std::endl;
		    }

		  if (_tpart == TrkParticle::mu_plus ||_tpart == TrkParticle::mu_minus)
		    {
		      double betaSpeed = _trkMom[ithtrk]/(sqrt(_trkMom[ithtrk]*_trkMom[ithtrk] + pow(pdt_.particle(TrkParticle::mu_minus).ref().mass(),2)));
		      flightTime = flightZ/(betaZ*3.0e11)/1.0e-9/betaSpeed; //in nsec
		    }

		  double timeAtDisk = _trkt0[ithtrk] + flightTime;
		  //		  std::cout << " time at disk and clutime = " << timeAtDisk << " " << _cluTime[ithclus] << " timediff = " << _timeDiff << std::endl;
		  if (abs(_cluTime[ithclus] - timeAtDisk) < _timeDiff){
		    //  std::cout << "dingding disk time " << hitFirstDisk << " " << hitSecondDisk << " " << _cluCogZ[ithclus] <<
		    //    " " << _firstDiskZ << " " << _secondDiskZ << std::endl;
		    if (_cluCogZ[ithclus] == _firstDiskZ) {std::cout << "first disk match "  << std::endl;}
		    if (_cluCogZ[ithclus] == _secondDiskZ) {std::cout << "second disk match " << std::endl;}
		    if (hitFirstDisk && (abs(_cluCogZ[ithclus] == _firstDiskZ) < _tol) )
		      {
			double dist = sqrt( (trkX-_cluCogX[ithclus])*(trkX-_cluCogX[ithclus]) + (trkY-_cluCogY[ithclus])*(trkY-_cluCogY[ithclus]) );
			//std::cout << "disk 1 dist = " << dist << std::endl;
			  if (dist < minDist)
			    {
			      ithClusMatch = ithclus;
			      minDist = dist;
			      foundMatch = true;
			      fracBest = _cluEnergy[ithclus]/_trkMom[ithtrk];
			      timeBest = _cluTime[ithclus] - (_trkt0[ithtrk]+flightTime);
			    }
		      }

		    if (hitSecondDisk &&(abs(_cluCogZ[ithclus] == _secondDiskZ) < _tol) )
		      {
			double dist = sqrt( (trkX-_cluCogX[ithclus])*(trkX-_cluCogX[ithclus]) + (trkY-_cluCogY[ithclus])*(trkY-_cluCogY[ithclus]) );
			//			std::cout << "disk 2 dist = " << dist << std::endl;
			  if (dist < minDist)
			    {
			      ithClusMatch = ithclus;
			      minDist = dist;
			      foundMatch = true;
			      fracBest = _cluEnergy[ithclus]/_trkMom[ithtrk];
			      timeBest = _cluTime[ithclus] - (_trkt0[ithtrk] + flightTime);
			    }
		      }
		  }

		}
	    }
	    if (foundMatch){
	      if (hitFirstDisk && fracBest > 0.8){
		_hShortestDistance->Fill(minDist);
		_htrkmomErr1->Fill(_trkmomErr[ithtrk]);
		_htrkt0Err1->Fill(_trkt0Err[ithtrk]);
		_htrkcon1->Fill(_trkcon[ithtrk]);	 
		_hEoverP1st->Fill(fracBest); _hTimeDiffBest1st->Fill(timeBest);
		std::cout << "and the 1st disk match says: " << foundMatch << " " << ithClusMatch << " " << minDist  << std::endl;
	      }
	      if (hitSecondDisk && fracBest > 0.8) {
		_hShortestFirst->Fill(minDist);
		_hEoverP2nd->Fill(fracBest);
		_htrkmomErr2->Fill(_trkmomErr[ithtrk]);
		_htrkt0Err2->Fill(_trkt0Err[ithtrk]);
		_htrkcon2->Fill(_trkcon[ithtrk]);	 
		_hTimeDiffBest2nd->Fill(timeBest);
		std::cout << "and the 2nd disk match says: " << foundMatch << " " << ithClusMatch << " " << minDist << std::endl;
	      }
	    }
    }
  }
}// end namespace mu2e

DEFINE_ART_MODULE(mu2e::KineticFracAnalysis);


