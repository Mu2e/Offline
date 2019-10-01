//
// an EDAnalyzer module that writes out calorimeter and tracker information so that it can be used to study combined track/calo information.  Largely based on Echenard's CaloExample_module.cc


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

#include "TrackerGeom/inc/Tracker.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CaloMC/inc/ClusterContentMC.hh"
#include "CaloMC/inc/CrystalContentMC.hh"

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
#include "RecoDataProducts/inc/TrackCaloAssnsCollection.hh"

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
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
    std::string _caloHitTruthModuleLabel;
    std::string _caloClusterTruthModuleLabel;
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
      ,*_hDistToTrackGoodCutC
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


    TH1F* _hZinter;
    TH1F* _hChisqDisk1GoodZ;
    TH1F* _hChisqDisk1BadZ;
    TH1F* _hDisk1GoodRadius;
    TH1F* _hDisk1BadRadius;
    TH1F* _hAngGood;
    TH1F* _hAngBad;
    TH1F* _hZForBadDisk1;

    TH1F* _hFracOnFrontFace;
    TH1F* _hFracOnEdgeFace;
    TH1F* _hFracOnEdgeFaceAndClusterInFront;
    TH1F* _hFracOnEdgeFaceAndClusterInBack;

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
    _caloHitTruthModuleLabel(pset.get<std::string>("caloHitTruthModuleLabel")),
    _caloClusterTruthModuleLabel(pset.get<std::string>("caloClusterTruthModuleLabel")),
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
    _htime = tfs->make<TH1F>("_htime","Time of Cluster - Time of Track at Calorimeter CoGZ",40, -20.,20.);
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
    _hDistToTrackGoodCutC = tfs->make<TH1F>("_hDistToTrackGoodCutC","Good Cut C, Distance of Cluster to Track at first disk",100,0.,1000.);

    _hFracOnFrontFace = tfs->make<TH1F>("_hFracOnFrontFace","E/p for front face hits disk 1", 150,0.,1.5);
    _hFracOnEdgeFace = tfs->make<TH1F>("_hFracOnEdgeFace","E/p for edge face hits disk 1", 150,0.,1.5);
    _hFracOnEdgeFaceAndClusterInFront = tfs->make<TH1F>("_hFracOnEdgeFaceAndClusterInFront","E/p for edge face hits disk 1, cluster front", 150,0.,1.5);
    _hFracOnEdgeFaceAndClusterInBack  = tfs->make<TH1F>("_hFracOnEdgeFaceAndClusterInBack", "E/p for edge face hits disk 1, cluster back", 150,0.,1.5);

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
    _hZinter = tfs->make<TH1F>("_hZinter","Intersection Z",350,1400.,2800.);
    _hChisqDisk1GoodZ = tfs->make<TH1F>("_hChisqDisk1GoodZ","Chisq for Hits on face of Front Disk",100,0.,100.);
    _hZForBadDisk1 = tfs->make<TH1F>("_hZForBadDisk1","Z of Best Cluster For Disk1 on inside",350,1400.,2800.);
    _hDisk1GoodRadius = tfs->make<TH1F>("_hDisk1GoodRadius","Radius for Hits on face of Front Disk",100,0.,1000.);
    //_hChisqDisk2GoodZ = tfs->make<TH1F>("_hChisqDisk2GoodZ",100,0.,100.);
    _hChisqDisk1BadZ = tfs->make<TH1F>("_hChisqDisk1BadZ","Chisq for hits on inside of Front Disk",100,0.,100.);
    _hDisk1BadRadius = tfs->make<TH1F>("_hDisk1BadRadius","Radius for Hits on inside of Front Disk",100,0.,1000.);

    _hAngGood = tfs->make<TH1F>("_hAngGood", "Angle wrt z at Front Face",20,0.,90.);
    _hAngBad  = tfs->make<TH1F>("_hAngBad",  "Angle wrt z at Inside Face",20,0.,90.);
    //_hChisqDisk2BadZ = tfs->make<TH1F>("_hChisqDisk2BadZ",100,0.,100.);
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

    if (_diagLevel > 0){std::cout << "******************new event in KineticFracAnalysis*******************" << std::endl;}


    //   ConditionsHandle<AcceleratorParams> accPar("ignored");
    //double _mbtime = accPar->deBuncherPeriod;
    //_toff.updateMap(event);

    ConditionsHandle<CalorimeterPhysicalConstants> calPhys("ignored");
    _density = calPhys->density();
    //    std::cout << "density is " << _density/(CLHEP::g/CLHEP::cm3) << std::endl;
    _eCritPos = calPhys->criticalEnergyPos();
    _eCritNeg = calPhys->criticalEnergyNeg();

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    _firstDiskZ  = cal.geomUtil().mu2eToTracker(cal.disk(0).geomInfo().frontFaceCenter()).z();
    _secondDiskZ = cal.geomUtil().mu2eToTracker(cal.disk(1).geomInfo().frontFaceCenter()).z();

    if (_diagLevel > 2){   std::cout << "checking disk locations " << _firstDiskZ << " " <<  _secondDiskZ << std::endl;}

    //Get handle to the tracker
    if( ! geom->hasElement<Tracker>() ) return;
	
    art::Handle<TrkCaloMatchCollection>  trkCaloMatchHandle;
    event.getByLabel(_trkCaloMatchModuleLabel, trkCaloMatchHandle);
    TrkCaloMatchCollection const& trkCaloMatches(*trkCaloMatchHandle);

    art::Handle<TrkCaloIntersectCollection>  trjIntersectHandle;
    event.getByLabel(_trkIntersectModuleLabel, trjIntersectHandle);
    TrkCaloIntersectCollection const& trkIntersect(*trjIntersectHandle);

    //Get generated particles
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);
    //GenParticleCollection const& genParticles(*gensHandle);


    //Calorimeter crystal hits (average from readouts)
    art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
    event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
    CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);

    //Calorimeter clusters
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    CaloClusterCollection const& caloClusters(*caloClustersHandle);


    //Calorimeter crystal truth assignment
    art::Handle<CaloClusterMCTruthAssns> caloClusterTruthHandle;
    event.getByLabel(_caloClusterTruthModuleLabel, caloClusterTruthHandle);
    const CaloClusterMCTruthAssns& caloClusterTruth(*caloClusterTruthHandle);

    //Calorimeter crystal truth assignment
    art::Handle<CaloHitMCTruthAssns> caloHitTruthHandle;
    event.getByLabel(_caloHitTruthModuleLabel, caloHitTruthHandle);
    const CaloHitMCTruthAssns& caloHitTruth(*caloHitTruthHandle);


    //Get virtual detector hits
    art::Handle<StepPointMCCollection> vdhits;
    event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);

    //Get tracks
    art::Handle<KalRepCollection> krepsHandle;
    event.getByLabel(_trkPatRecModuleLabel,_instanceName,krepsHandle);
    KalRepCollection const& kreps = *krepsHandle;


    //
    // handle to PDG
    //    GlobalConstantsHandle<ParticleDataTable> pdt;
    //ParticleDataTable const & pdt_ = *pdt;
    

    //      const double CrDensity = 4.9*(CLHEP::g/CLHEP::cm3);
    //const double CrMass    = CrDensity*cal.caloInfo().crystalVolume();

    double numberOfTracks = 0;
    double numberOfClusters = 0;



    //--------------------------  Do generated particles --------------------------------


    _evt = event.id().event();
    _run = event.run();

    if (_diagLevel > 2){std::cout << "processing event in Kinetic_Frac " << _nProcess << " run and event  = " << _run << " " << _evt << " with instance name = " << _instanceName << std::endl;}

    //--------------------------  Do calorimeter hits --------------------------------

    _nHits = _nSim = 0;
    _cryEtot = 0.0;

    for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
    {
        const CaloCrystalHit &hit     = caloCrystalHits.at(ic);
	int diskId                    = cal.crystal(hit.id()).diskId();
        CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position()); //in disk FF frame

        CrystalContentMC contentMC(cal, caloHitTruth, hit);

        _cryEtot             += hit.energyDep();
        _cryTime[_nHits]      = hit.time();
        _cryEdep[_nHits]      = hit.energyDep();
        _cryPosX[_nHits]      = crystalPos.x();
        _cryPosY[_nHits]      = crystalPos.y();
        _cryPosZ[_nHits]      = crystalPos.z();
        _cryId[_nHits]        = hit.id();
        _crySectionId[_nHits] = diskId;

        _crySimIdx[_nCluster] = _nCluSim;
        _crySimLen[_nCluster] = contentMC.simContentMap().size();

        for (const auto& contentMap : contentMC.simContentMap() )
	{	       
	    art::Ptr<SimParticle> sim = contentMap.first;
	    CaloContentSim       data = contentMap.second;

            _motId[_nSim]      = sim->id().asInt();
            _motPdgId[_nSim]   = sim->pdgId();
            _motmom[_nSim]     = data.mom();
            _motcrCode[_nSim]  = sim->creationCode();
       	    _motTime[_nSim]    = data.time();
            _motEdep[_nSim]    = data.edep();	       
            ++_nSim;
         }
        ++_nHits;
    }

    //--------------------------  Do clusters --------------------------------

    _nCluster = _nCluSim = 0;
    _cluList.clear();
    for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
    {

        ClusterContentMC contentMC(cal, caloClusterTruth, *clusterIt);

        std::vector<int> _list;
        for (int i=0;i<clusterIt->size();++i)
        {
            int idx = int(clusterIt->caloCrystalHitsPtrVector().at(i).get()- &caloCrystalHits.at(0));
            _list.push_back(idx);
        }

        _cluEnergy[_nCluster] = clusterIt->energyDep();
        _cluTime[_nCluster]   = clusterIt->time();
        _cluNcrys[_nCluster]  = clusterIt->size();
        _cluCogX[_nCluster]   = clusterIt->cog3Vector().x(); //in disk FF frame
        _cluCogY[_nCluster]   = clusterIt->cog3Vector().y();
        _cluCogZ[_nCluster]   = clusterIt->cog3Vector().z();
        _cluConv[_nCluster]   = (contentMC.hasConversion() ? 1 : 0);
        _cluList.push_back(_list);

        _cluSimIdx[_nCluster] = _nCluSim;
        _cluSimLen[_nCluster] = contentMC.simContentMap().size();

        for (const auto& contentMap : contentMC.simContentMap() )
	{	       
	    art::Ptr<SimParticle> sim = contentMap.first;
	    CaloContentSim       data = contentMap.second;

	    art::Ptr<SimParticle> smother(sim);
            while (smother->hasParent()) smother = smother->parent();
            int genIdx=-1;
            if (smother->genParticle()) genIdx = smother->genParticle()->generatorId().id();

            _clusimId[_nCluSim]     = sim->id().asInt();
            _clusimPdgId[_nCluSim]  = sim->pdgId();
            _clusimGenIdx[_nCluSim] = genIdx;
            _clusimTime[_nCluSim]   = data.time();
            _clusimEdep[_nCluSim]   = data.edep();
            _clusimMom[_nCluSim]    = data.mom();

            ++_nCluSim;
         }

        ++_nCluster;
    }


    numberOfClusters = _nCluster;

    //
    // for matching tracks to clusters
    int _numberOfTracks = kreps.size();
    std::vector<double>    bestChi2(_numberOfTracks,-1.);
    std::vector<double>    bestChi2Pos(_numberOfTracks,-1.);
    std::vector<double>    bestChi2Time(_numberOfTracks,-1.);
    std::vector<int>       bestCluster(_numberOfTracks,-1);
    std::vector<int>       bestMatch(_numberOfTracks,-1);
 
    for (int ithTrack = 0; ithTrack < _numberOfTracks; ++ithTrack)
      {
	double minChi2 = -1.; int minMatch = -1; int minCluster = -1; double minChi2Pos = -1.; double minChi2Time = -1.;
	int _nMatch(0);
	for (auto const& trkCaloMatch: trkCaloMatches)
	  {
	    //
	    // match this track only
	    if (trkCaloMatch.trkId() == ithTrack)
	      {
		_mTrkId[_nMatch]     = trkCaloMatch.trkId();
		_mCluId[_nMatch]     = trkCaloMatch.cluId();
		_mChi2[_nMatch]      = trkCaloMatch.chi2();
		_mChi2Pos[_nMatch]   = trkCaloMatch.chi2Pos();
		_mChi2Time[_nMatch]  = trkCaloMatch.chi2Time();
		if (trkCaloMatch.chi2() < minChi2 || minChi2 == -1)
		  {
		    minCluster = trkCaloMatch.cluId();
		    minChi2 = trkCaloMatch.chi2();
		    minChi2Pos = trkCaloMatch.chi2Pos();
		    minChi2Time = trkCaloMatch.chi2Time();
		    minMatch = _nMatch;
		  }
		if (_diagLevel > 3)
		  {std::cout << "nmatch, track ID, cluster ID, chi2, cluster time " << _nMatch << " " << _mTrkId[_nMatch] << " " 
			     << _mCluId[_nMatch] << " "  << _mChi2[_nMatch] << " " << _cluTime[_mCluId[_nMatch]]  << std::endl;
		  }
		++_nMatch;
	      }
	  }
	bestChi2[ithTrack] = minChi2;
	bestChi2Pos[ithTrack] = minChi2Pos;
	bestChi2Time[ithTrack] = minChi2Time;
	bestCluster[ithTrack] = minCluster;
	bestMatch[ithTrack] = minMatch;
	if (_diagLevel > 2)
	  {std::cout << "for track " << ithTrack << " nmatch, track ID, cluster ID, chi2, cluster time " << minMatch << " " << _mTrkId[minMatch] << " " 
		     << _mCluId[minMatch] << " "  << _mChi2[minMatch] << " " << _cluTime[_mCluId[minMatch]]  << std::endl;
	  }
	if (_diagLevel > 2){std::cout << "number of clusters = " << numberOfClusters << " number of matches = " << _nMatch << std::endl;}

	//
	// there is (2/23/2016) a code infelicity in TrackCaloMatchingBis.  If the track hits the first disk, the intersection module
	// still propagates the track to the second disk.  In that case the matching module matches
	// clusters to each intersection, and you get the same matches twice: _nCluster for the first disk, then the cluster list repeated
	// for the second disk.  Vice-versa for upstream tracks. In either case, take the first set and just use that.  The algorithm will be upgraded.
	// perhaps then there can be more than two intersections (one per disk) per track but the throw will warn we need a code upgrade.
	// 
	// and then sometimes you find more than one track (as of this writing the tracks are found in order of time, but that can't be assumed)
	//
	// hence I put in:
	if (_diagLevel > 2)
	  {
	    std::cout << " _nMatch = " << _nMatch << std::endl;
	    std::cout << " _nCluster = " << _nCluster << std::endl;
	    std::cout << " kreps.size() = " << kreps.size() << std::endl;
	  }
	if (_nMatch > 2*_nCluster*static_cast<int>(kreps.size()) )
	  { throw cet::exception("RANGE") << " _nMatch and _nCluster mismatch, from KineticFracAnalysis. Call Bob ...nmatch and ncluster are " <<  _nMatch << " " << _nCluster << std::endl;}

      }


    _nTrk=0;


	
    if (_diagLevel > 3)
      {
	std::cout << "number of tracks = "<< kreps.size() << std::endl;
	std::cout << "size of intersection collection = " << trkIntersect.size() << std::endl;
      }
    if (trkIntersect.size() == 0 ) return; // code purists look the other way
    if (trkIntersect.size()  > 2*kreps.size()) // can't have more than two intersections per track
      { 
	throw cet::exception("RANGE") << " more than two intersections for track and calorimeter, from KineticFracAnalysis. Call Bob or Bertrand..." << std::endl;
      }
 
    //
    // we now have matches to each track.  Have to line these up with TrackIntersectionCollection.  Trivial if one track
    // but some events (small to be sure) have multiple tracks.  Therefore loop over TrackIntersectionCollections and see which track
    // a given intersection goes with, then analyze that.

    int ithInt(0);
    for (auto const& trackCaloIntersect: trkIntersect)
      {
	if (ithInt > 0) {break;}
	++ithInt;
	int trackNumber = trackCaloIntersect.trkId();
	KalRep const &trk = *(trackCaloIntersect.trk());
	double pathLength = trackCaloIntersect.pathLengthEntrance();
	HepPoint point = trk.position(pathLength);
	std::cout << "path length = " << pathLength << std::endl;
	double trkArrivalTime = trk.arrivalTime(pathLength);
	if (_cluCogZ[bestCluster[trackNumber]] > 1800.) {break;}
	double deltaTime = _cluTime[bestCluster[trackNumber]] - trkArrivalTime;



	//
	// this is for comparing to clusters

	HelixTraj trkHel(trk.helix(pathLength).params(),trk.helix(pathLength).covariance());
	double trkTime                = trk.arrivalTime(pathLength);
	double tPhi     = trkHel.phi0();
	double td0     = trkHel.d0();
	double tz0     = trkHel.z0();
	double tOmega   = trkHel.omega();
	double tCosDip  = trkHel.cosDip();
	double tProb                  = trk.chisqConsistency().significanceLevel();
	int    tNhits                 = trk.nActive();
	int    tStatus                = trk.fitStatus().success();
	CLHEP::Hep3Vector trkMomentum = trk.momentum(pathLength);
	HepVector momvec(3);//establshing dimension for next line
	momvec = trkMomentum.unit();
	BbrVectorErr momCov = trk.momentumErr(0);
	double fitMomErr    = sqrt(momCov.covMatrix().similarity(momvec));


	double tanDip = sqrt( 1. - tCosDip*tCosDip)/tCosDip;

	      
	double angleWrtZ = acos(trkMomentum.z()/trkMomentum.mag())*(180/3.14159);

	std::cout << "delta time = " << _cluTime[bestCluster[trackNumber]] << " " << trkArrivalTime << " " << deltaTime << std::endl;
	//	if (deltaTime > -4){break;}
	if (_diagLevel > 2)
	  {
	    std::cout << "*****ith Int = " << ithInt << "**********" << std::endl;
	    std::cout << "trackID, pathLengthEntrance, point: " << trackNumber << " " << pathLength << " " << point << std::endl;
	  }

	//
	// only plot these for the disk that goes with the best match cluster

	double bestZ = _cluCogZ[bestCluster[trackNumber]];

	if (_diagLevel > 2){std::cout << " best Z is " << bestZ << " and point z = " << point.z() << std::endl;}

	_hZinter->Fill(point.z());
	double rad1 = sqrt(point.x()*point.x() + point.y()*point.y());

	// 
	//just look at front disk
	if (point.z() < 1800.)
	  {
	//    double rad2 = sqrt(pointFront.x()*pointBack.x() + pointFront.y()*pointBack.y());
	    if (abs(point.z() - bestZ) < 1) {
	      if (_diagLevel > 2)
		{
		  std::cout << "in good rad, chi2 = " << rad1 << " " << bestChi2[trackNumber]  
			    << std::endl;
		}
	      _hChisqDisk1GoodZ->Fill(bestChi2Time[trackNumber]); _hDisk1GoodRadius->Fill(rad1);
	      _hFracOnFrontFace->Fill(_cluEnergy[bestCluster[trackNumber]]/trk.momentum(pathLength).mag());
	      _hAngGood->Fill(angleWrtZ);
	    }
	    if (abs(point.z() - bestZ) > 30) 
	      {
		if (_diagLevel > 2)
		  {
		    std::cout << "in bad rad, chi2 = "  << rad1 << " " << bestChi2[trackNumber] << std::endl;
	    std::cout << "path lengths bad rad = " << trackCaloIntersect.pathLengthEntrance() << std::endl;
	    std::cout << "point bad rad = " << point << std::endl;
	    std::cout << " energy of best cluster match = " << _cluEnergy[bestCluster[trackNumber]] << " and momentum = " << trk.momentum(pathLength).mag() << std::endl; 
		  }
		_hAngBad->Fill(angleWrtZ);
		_hChisqDisk1BadZ->Fill(bestChi2Time[trackNumber]);
		_hDisk1BadRadius->Fill(rad1);
		_hZForBadDisk1->Fill(_cluCogZ[bestCluster[trackNumber]]);
		_hFracOnEdgeFace->Fill(_cluEnergy[bestCluster[trackNumber]]/trk.momentum(pathLength).mag());
		if (bestZ < 2000.){_hFracOnEdgeFaceAndClusterInFront->Fill(_cluEnergy[bestCluster[trackNumber]]/trk.momentum(pathLength).mag());}
		if (bestZ > 2000.){_hFracOnEdgeFaceAndClusterInBack->Fill(_cluEnergy[bestCluster[trackNumber]]/trk.momentum(pathLength).mag());
		}
	      }
	  }

	if (_diagLevel > 2)
	  {	
	    std::cout << "path lengths = " << trackCaloIntersect.pathLengthEntrance() << std::endl;
	    std::cout << "point = " << point << std::endl;
	  }
	      
	//
	// these are for cut set c comparisons; numbers aren't identical to Dave Brown's but are very close. Not worth going back into kreps
	//
	//after extendvalidrange is incorporated this fudge will be unnecessary.  Right now Bertrand has cast away original krep constancy so who knows
	//what we have.
	HelixTraj trkHel0(trk.helix(0).params(),trk.helix(0).covariance());

		  
	int    tNHits0    = trk.nActive();
	double tPhi0     = trkHel0.phi0();
	double td00     = trkHel0.d0();
	double tz00     = trkHel0.z0();
	double tOmega0   = trkHel0.omega();
	double tCosDip0  = trkHel0.cosDip();
	int    tStatus0     = trk.fitStatus().success();
	double tProb0    = trk.chisqConsistency().significanceLevel(); //is this right?
	double trkTime0 = trk.arrivalTime(0);
	CLHEP::Hep3Vector trkMomentum0 = trk.momentum(0);
	HepVector momvec0(3);
	momvec0 = trkMomentum0.unit();
	BbrVectorErr momCov0 = trk.momentumErr(0);
	double fitMomErr0    = sqrt(momCov0.covMatrix().similarity(momvec0));
 	
	if (tCosDip0 == 0)
	  { throw cet::exception("RANGE") << " cosDip = 0 for track, which makes no sense, from KineticFracAnalysis...it's Bobs fault" << std::endl;}

	double tanDip0 = sqrt( 1. - tCosDip0*tCosDip0)/tCosDip0;
	CLHEP::Hep3Vector momVec0 = trk.momentum(0);


	// Does this fit pass cut set C?
	bool cutC = ( tStatus0 >0)
	  && ( tProb0             > 2.e-3  )
	  && ( tNHits0    >= 20    ) 
	  && ( fitMomErr0          < 0.25   ) 
	  //      && ( trkTime0 > 700 && trkTime0 < 1695)
	  ;
	//    if (! cutC) return;  
	//
	// look at clusters to see if any match

	//double timeDiff = 1696;
	_htime->Fill(_cluTime[bestCluster[trackNumber]] - trkTime);
	if (_diagLevel > 2) {std::cout << " cluster - track time = " << _cluTime[bestCluster[trackNumber]] << " " << trkTime<< " "  << _cluTime[bestCluster[trackNumber]] - trkTime << std::endl;}

	//
	// what's the distance for this cluster?

	double _dist = sqrt(pow( _cluCogX[bestCluster[trackNumber]] - point.x() ,2) + pow(_cluCogY[bestCluster[trackNumber]] - point.y(),2) );
	_hDistToTrack->Fill(_dist);
	if (! cutC){_hDistToTrackGoodCutC->Fill(_dist);}

    HepPoint point0 = trk.position(0);


	if (_diagLevel >2 )
	  {
	    std::cout << "and the cluster match: best cluster = " << bestCluster[trackNumber] << " for track " << trackNumber << "  \n"
	      "cluCogX = " << _cluCogX[bestCluster[trackNumber]] << " cluCogY = " << _cluCogY[bestCluster[trackNumber]] << " cluCogZ " << _cluCogZ[bestCluster[trackNumber]] << " \n"
	      " point is " << point << std::endl;
	  }

	if (_diagLevel > 2 && _dist > 350)
	  {
	    std::cout << " big distance comparisons: dist = " << _dist << "\n"
		      << " tStatus    = " << tStatus << " " << tStatus0  << " \n" 
		      << " tProb      = " << tProb << " " << tProb0  << " \n"
		      << " nhits      = " << tNhits << " " << tNHits0  << " \n"
		      << " track time = " << trkTime << " " << trkTime0  << " \n"
		      << " omega      = " << tOmega << " " << tOmega0  << " \n"
		      << " phi        = " << tPhi << " " << tPhi0 << " \n"
		      << " d0         = " << td0 << " " << td00  << " \n"
		      << " z0         = " << tz0 << " " << tz00  << " \n"
		      << " tandip     = " << tanDip << " " << tanDip0 << " \n"
		      << " momentum   = " << trkMomentum << " " << trkMomentum0 << " \n"
		      << " magnitude  = " << trkMomentum.mag() << " " << trkMomentum0.mag() << "\n"
		      << " mom err    = " << fitMomErr <<  " " << fitMomErr0 << "\n"
		      << " cut set C  = " << cutC 
		      << std::endl;

	  }  
      }
 

     

    

    //--------------------------  Do tracks  --------------------------------
    _nTrkOk = 0;
    _nTrk = 0;

    for ( KalRepCollection::const_iterator i=kreps.begin(); i != kreps.end(); ++i ){

      KalRep const& krep = **i;

      // Arc length from center of tracker to the most upstream point on the fitted track.

      double s0 = krep.startValidRange(); // worst possible place to use, should start at tracker.  Dave has upgraded code
      //to extend track valid range all the way through calorimeter.  If you want to get track at s0 at intersection, have
      // to set valid range parameter in .fcl.  DOn't know s0 though until we've gotten the intersection!! so
      // probably setting this to 0 at center of tracker and then extrapolating from there is good enough.  Use Bertrand's extrapolation
      //since he does it all...
      s0= 0.;
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
      bool cutC = ( krep.fitStatus().success() >0) 
	&& ( fitCon             > 2.e-3  ) 
	&& ( krep.nActive()    >= 20    ) 
	&& ( fitMomErr          < 0.25   ) 
	//	( krep.t0().t0() > 700 && krep.t0().t0() < 1695)
	;

      if (!cutC && _diagLevel>2) {std::cout << "failed cutset" << std::endl;}

      //      cutC = true;
      
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


      if (_diagLevel > 4)
	{
	  std:: cout << "KineticFracAnalysis _nTrk, t0 = " << _nTrk << " " << krep.t0().t0() << " and momentum " << p0.mag() << std::endl;
	  std:: cout << "KineticFracAnalysis fitstatus = " << krep.fitStatus().success() << std::endl;
	}

      if (cutC)  ++_nTrkOk;
    }
    numberOfTracks = _nTrk;
    if (_diagLevel>3){std::cout << "ntracks = " << _nTrk << std::endl;}
    //--------------------------  Do tracker hits  --------------------------------

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


    _goodKfrac = 0;

    if (numberOfTracks   == 0){_goodKfrac =1;}
    if (numberOfClusters == 0 && numberOfTracks > 0){_goodKfrac = 2;}
    if (numberOfClusters > 0 && numberOfTracks == 0){_goodKfrac = 3;}
       std::cout << "number of tracks = " << numberOfTracks << std::endl;
    for (int ithtrk =  0; ithtrk < numberOfTracks;    ++ithtrk) {
    }
    // _Ntup->Fill();
  }
}// end namespace mu2e
DEFINE_ART_MODULE(mu2e::KineticFracAnalysis);


