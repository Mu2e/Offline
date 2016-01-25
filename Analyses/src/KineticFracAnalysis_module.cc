//
// An EDAnalyzer module that writes out calorimeter and tracker information so that it can be used to study combined track/calo information.  Largely based on Echenard's CaloExample_module.cc


//
// Original author Robert Bernstein
//

#include "CLHEP/Units/SystemOfUnits.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CalorimeterPhysicalConstants.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"


#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CaloCluster/inc/CaloContentMC.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
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
    std::string _caloClusterAlgorithm;
    std::string _caloClusterSeeding;
    const std::string _producerName;
    std::string _virtualDetectorLabel;
    std::string _stepPointMCLabel;
    std::string _trkPatRecModuleLabel;
    std::string _instanceName;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    SimParticleTimeOffset _toff;  // time offset smearing

    std::string _shLabel;
    std::string _shpLabel; 
    std::string _shfLabel;
    std::string _bkfLabel;

    const StrawHitCollection*         _shcol;
    const StrawHitFlagCollection*     _shfcol;
    const StrawHitFlagCollection*     _bkfcol;
    const StrawHitPositionCollection* _shpcol;



    TH1F *_hcryE,*_hcryT,*_hcryX,*_hcryY,*_hcryZ;
    TH1F *_hcluE,*_hcluT,*_hcluX,*_hcluY,*_hcluZ,*_hcluE1Et,*_hcluE1E9,*_hcluE1E25;


    TH1F *_hfitstatus  ,*_htrkcon
      ,*_htrkmomErr, *_htrkt0
      ,*_htrkt0Err
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
      ,*_hShortestDistance
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
      ,*_hEoverP
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
       

    int  _goodKfrac;

    //
    //for BaF2
    double _eCrit = 1./ (0.7833/(610/(56+1.24)) + 0.2167/(610/(9.+1.24)));

    //    CLHEP::Hep3Vector _firstDiskLoc;
    //CLHEP::Hep3Vector _secondDiskLoc;
    double _firstDiskZ;
    double _secondDiskZ;

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
  }

  void KineticFracAnalysis::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("KineticFrac", "KineticFrac");



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
    _htrkmomErr    = tfs->make<TH1F>("_htrkmomErr","trkmomErr",100,0.,2.);
    _htrkt0           = tfs->make<TH1F>("_htrkt0","trkt0",170,0.,1700.);
    _htrkt0Err        = tfs->make<TH1F>("_htrkt0Err","trkt0Err",100,0.,10.);
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
   _hEoverP = tfs->make<TH1F>("_hEoverP","E/p for shortest distance cluster", 150,0.,1.5);
   _hEOverPZoom = tfs->make<TH1F>("_hEOverPZoom","E/p for shortest distance cluster", 50,0.8,1.3);
   _hDistToTrack = tfs->make<TH1F>("_hDistToTrack","Distance of Cluster to Track at first disk",100,0.,1000.);

   _hNumberOfHits = tfs->make<TH1F>("_hNumberOfHits","Total Number of Hits",100,0.,5000.);
   _hNumberOfHitsAfter650Nsec = tfs->make<TH1F>("_hNumberOfHitsAfter650Nsec","Total Number of Hits After 650 nsec",100,0.,5000.);
   _hTimeDiff = tfs->make<TH1F>("_hTimeDiff","time diff to t0",100, -200.,200.);

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
      double density = calPhys->density();
      std::cout << "density is " << density << std::endl;

      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      _firstDiskZ  = cal.toTrackerFrame(cal.section(0).frontFaceCenter()).z();
      _secondDiskZ = cal.toTrackerFrame(cal.section(1).frontFaceCenter()).z();

      std::cout << "checking disk locations " << _firstDiskZ << " " <<  _secondDiskZ << std::endl;

      std::cout << " and other end " << cal.toTrackerFrame(cal.section(0).frontFaceCenter()) << " " << cal.toTrackerFrame(cal.section(1).frontFaceCenter())<<std::endl;
      //Get handle to the tracker
      if( ! geom->hasElement<TTracker>() ) return;
      //      TTracker const & tracker = *(GeomHandle<TTracker>());

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
	   std::cout << "cluster number = " << _nCluster << " " << clusterIt->time() << std::endl;
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

	   std::cout << " crystal position in z " << _cluCogZ[_nCluster] << std::endl; 
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
           ( fitMomErr          < 0.25   );

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
    std::cout << "number of tracks = " << numberOfTracks << std::endl;
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
	    _trkz0[ithtrk] = _trkz0[ithtrk] + 10175.;

	    //
	    //known offset from pezzullo, see his producer file
	    _trkt0[ithtrk] -= 1.4;
	  }	
    }


    for (int ithtrk =  0; ithtrk < numberOfTracks;    ++ithtrk) {

      // first thing to do is extrapolate track to the disks.  Does it get through first disk and end up striking second disk?

	//
	// extrapolate to first disk.  just use the helix since the spread on the 
	//tracker is small compared to the location from the calorimeter.  On the 
	//other hand, knowing from the track where the shower hit might be helpful 
	//in making energy corrections for energy going down cracks.  
	//So eventually will have to look at both

      std::cout << "got a track" << std::endl;

      double betaZ = _trkDip[ithtrk]/sqrtOrThrow(1. + _trkDip[ithtrk]*_trkDip[ithtrk],0.00001);

     for (int ithclus = 0; ithclus < numberOfClusters; ++ithclus){

       std::cout << "got a cluster" << std::endl;
       //if this particle is being fitted under an e+/e- hypothesis, want to make a correction;
       //if under a muon hypothesis, correction is zero
       Float_t showerCorr;
       if (_tpart == TrkParticle::e_minus || _tpart == TrkParticle::e_plus)
	 { 

	   std::cout << "in here" << std::endl;
	   showerCorr = (  log(_cluEnergy[ithclus]/_eCrit) - 0.1 )*20.3; //adjusts for shower development, formula out of Rossi
	 } else {showerCorr = 0.;}

	    // showerCorr along axis; how far to project in z?
	    showerCorr *= betaZ;
	    double zClus = _cluCogZ[ithclus]; //need to add code for which disk this is, use right coordinates, etc.  placeholder

	    std::cout << "zclus is " << zClus << std::endl;
	    zClus += showerCorr;

	    // l parameter 
	    //
	    //switch z0 to Mu2e coordinates
	    double flightLength = (zClus - _trkz0[0])/_trkDip[0]; 
	    //
	    // some random code to avoid compile error
	    flightLength += 0.;

 

      }    
    }
   
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::KineticFracAnalysis);


