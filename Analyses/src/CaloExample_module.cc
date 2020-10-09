//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
//
// Original author Bertrand Echenard
//

#include "CLHEP/Units/SystemOfUnits.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "CaloMC/inc/ClusterContentMC.hh"
#include "CaloMC/inc/CrystalContentMC.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

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

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>



namespace mu2e {


  class CaloExample : public art::EDAnalyzer {

     public:

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;



       explicit CaloExample(fhicl::ParameterSet const& pset);
       virtual ~CaloExample() { }

       virtual void beginJob();
       virtual void endJob();

       // This is called for each event.
       virtual void analyze(const art::Event& e);





     private:

       int _diagLevel;
       bool _doGenerated;
       int _nProcess;

       std::string _g4ModuleLabel;
       std::string _generatorModuleLabel;
       art::InputTag _simParticleTag;


       std::string _caloCrystalModuleLabel;
       std::string _caloClusterModuleLabel;
       std::string _caloHitTruthModuleLabel;
       std::string _caloClusterTruthModuleLabel;
       std::string _caloClusterAlgorithm;
       std::string _caloClusterSeeding;
       const std::string _producerName;
       std::string _virtualDetectorLabel;
       std::string _stepPointMCLabel;
       std::string _trkPatRecModuleLabel;
       SimParticleTimeOffset _toff;  // time offset smearing


       TH1F *_hcryE,*_hcryT,*_hcryX,*_hcryY,*_hcryZ;
       TH1F *_hcluE,*_hcluT,*_hcluX,*_hcluY,*_hcluZ,*_hcluE1Et,*_hcluE1E9,*_hcluE1E25,*_hcluEF;
       
       TH2F *_hxy;

       TTree* _Ntup;


       int   _evt,_run;

       int   _nGen,_genPdgId[16384],_genCrCode[16384];
       float _genmomX[16384],_genmomY[16384],_genmomZ[16384],_genStartX[16384],_genStartY[16384],_genStartZ[16384],_genStartT[16384];

       int   _nHits,_cryId[163840],_crySectionId[163840],_crySimIdx[163840],_crySimLen[163840];
       float _cryEtot,_cryTime[163840],_cryEdep[163840],_cryDose[163840],_cryPosX[163840],_cryPosY[163840],_cryPosZ[163840],_cryLeak[163840];

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





  };


  CaloExample::CaloExample(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _doGenerated(pset.get<bool>("doGenerated",false)),
    _nProcess(0),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel")),
    _simParticleTag(pset.get<std::string>("simParticleTag")),
    _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _caloHitTruthModuleLabel(pset.get<std::string>("caloHitTruthModuleLabel")),
    _caloClusterTruthModuleLabel(pset.get<std::string>("caloClusterTruthModuleLabel")),
    _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName")),
    _stepPointMCLabel(pset.get<std::string>("stepPointMCLabel")),
    _trkPatRecModuleLabel(pset.get<std::string>("trkPatRecModuleLabel")),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _Ntup(0)

  {
  }

  void CaloExample::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;

       _Ntup  = tfs->make<TTree>("Calo", "Calo");



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

       _hcryE     = tfs->make<TH1F>("cryEdep",  "Energy deposited / crystal", 100,    0., 50.   );
       _hcryT     = tfs->make<TH1F>("cryTime",  "Time of crystal hit",        100,    0., 2000. );
       _hcryX     = tfs->make<TH1F>("cryX",     "X coord of crystal hit",     100,  300., 700.  );
       _hcryY     = tfs->make<TH1F>("cryY",     "Y coord of crystal hit",     100,  300., 700.  );
       _hcryZ     = tfs->make<TH1F>("cryZ",     "Z coord of crystal hit",     100,11000., 13000.);
       _hcluE     = tfs->make<TH1F>("cluEdep",  "Energy deposited / clustal", 150,    0., 150.  );
       _hcluEF    = tfs->make<TH1F>("cluEdepF", "Energy deposited / clustal", 150,    0., 150.  );
       _hcluT     = tfs->make<TH1F>("cluTime",  "Time of clustal hit",        100,    0., 2000. );
       _hcluX     = tfs->make<TH1F>("cluX",     "X coord of clustal hit",     100,  300., 700.  );
       _hcluY     = tfs->make<TH1F>("cluY",     "Y coord of clustal hit",     100,  300., 700.  );
       _hcluZ     = tfs->make<TH1F>("cluZ",     "Z coord of clustal hit",     100,11000., 13000.);
       _hcluE1Et  = tfs->make<TH1F>("cluE1Et",  "E1/Etot",                    100,    0., 1.1   );
       _hcluE1E9  = tfs->make<TH1F>("cluE1E9",  "E1/E9",                      100,    0., 1.1   );
       _hcluE1E25 = tfs->make<TH1F>("cluE1E25", "E1/E25",                     100,    0., 1.1   );
       _hxy       = tfs->make<TH2F>("cryxy", "cryxy",                     350,-700,700,350,-700,700  );

  }



  void CaloExample::endJob(){
  }




  void CaloExample::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloExample =  "<<_nProcess <<std::endl;

      ConditionsHandle<AcceleratorParams> accPar("ignored");
      double _mbtime = accPar->deBuncherPeriod;
      _toff.updateMap(event);


      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
      CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);

      //Virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);

      //Calorimeter crystal truth assignment
      art::Handle<CaloClusterMCTruthAssns> caloClusterTruthHandle;
      event.getByLabel(_caloClusterTruthModuleLabel, caloClusterTruthHandle);
      const CaloClusterMCTruthAssns& caloClusterTruth(*caloClusterTruthHandle);

      //Calorimeter crystal truth assignment
      art::Handle<CaloHitMCTruthAssns> caloHitTruthHandle;
      event.getByLabel(_caloHitTruthModuleLabel, caloHitTruthHandle);
      const CaloHitMCTruthAssns& caloHitTruth(*caloHitTruthHandle);

//      art::Handle<KalRepPtrCollection> trksHandle;
//      event.getByLabel(_trkPatRecModuleLabel, trksHandle);
//      const KalRepPtrCollection& trks = *trksHandle;


      const double CrDensity = 4.9*(CLHEP::g/CLHEP::cm3);
      const double CrMass    = CrDensity*cal.caloInfo().crystalVolume();


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
       } else {_nGen=0;}


       //--------------------------  Do calorimeter hits --------------------------------

       _nHits = _nSim = 0;
       _cryEtot = 0.0;

       for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
       {
           const CaloCrystalHit &hit     = caloCrystalHits.at(ic);
	   int diskId                    = cal.crystal(hit.id()).diskId();
           CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  //in disk FF frame
 
           CrystalContentMC contentMC(cal, caloHitTruth, hit);

           _cryEtot             += hit.energyDep();
           _cryTime[_nHits]      = hit.time();
           _cryEdep[_nHits]      = hit.energyDep();
           _cryDose[_nHits]      = hit.energyDep() / CrMass / (CLHEP::joule/CLHEP::kg); //dose
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
               
	       auto parent(sim);
               while ( parent->hasParent()) parent = parent->parent();               
               
	       _motId[_nSim]      = sim->id().asInt();
               _motPdgId[_nSim]   = sim->pdgId();
               _motmom[_nSim]     = data.mom();
               _motcrCode[_nSim]  = sim->creationCode();
       	       _motTime[_nSim]    = data.time();
               _motEdep[_nSim]    = data.edep();	                      
	       
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
           _cluE1[_nCluster]     = clusterIt->e1();
           _cluE9[_nCluster]     = clusterIt->e9();
           _cluE25[_nCluster]    = clusterIt->e25();
           _cluSecMom[_nCluster] = clusterIt->secondMoment();
           _cluSplit[_nCluster]  = clusterIt->isSplit();
           _cluConv[_nCluster]   = (contentMC.hasConversion() ? 1 : 0);
           _cluList.push_back(_list);

           _cluSimIdx[_nCluster] = _nCluSim;
           _cluSimLen[_nCluster] = contentMC.simContentMap().size();

           for (const auto& contentMap : contentMC.simContentMap() )
	   {	       
	       art::Ptr<SimParticle> sim = contentMap.first;
	       CaloContentSim       data = contentMap.second;

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
		  CLHEP::Hep3Vector simPos = cal.geomUtil().mu2eToDiskFF(clusterIt->diskId(), vdMapEntry->second->position());		  
	       } 

               _clusimId[_nCluSim]     = sim->id().asInt();
               _clusimPdgId[_nCluSim]  = sim->pdgId();
               _clusimGenId[_nCluSim]  = genId;
	       _clusimGenPdg[_nCluSim] = genPdg;
               _clusimCrCode[_nCluSim] = sim->creationCode();
               _clusimTime[_nCluSim]   = data.time();
               _clusimEdep[_nCluSim]   = data.edep();
               _clusimMom[_nCluSim]    = data.mom();
               _clusimMom2[_nCluSim]   = simMom;
               _clusimPosX[_nCluSim]   = simPos.x(); // in disk FF frame
               _clusimPosY[_nCluSim]   = simPos.y();
               _clusimPosZ[_nCluSim]   = simPos.z();  
               _clusimStartX[_nCluSim] = sim->startPosition().x(); // in disk FF frame
               _clusimStartY[_nCluSim] = sim->startPosition().y();
               _clusimStartZ[_nCluSim] = sim->startPosition().z();  

               ++_nCluSim;
            }

           _hcluE->Fill(clusterIt->energyDep());
           _hcluT->Fill(clusterIt->time());
           _hcluX->Fill(clusterIt->cog3Vector().x());
           _hcluY->Fill(clusterIt->cog3Vector().y());
           _hcluZ->Fill(clusterIt->cog3Vector().z());
	   
	  
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

               double hitTimeUnfolded = _toff.timeWithOffsetsApplied(hit);
   	       double hitTime         = fmod(hitTimeUnfolded,_mbtime);

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


       //--------------------------  Do tracks  --------------------------------
       
       _nTrk = 0;
/*
       for ( size_t itrk=0; itrk< trks.size(); ++itrk )
       {
         KalRep const* krep = trks.at(itrk).get();

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
*/        


 
        _Ntup->Fill();





  }



}  

DEFINE_ART_MODULE(mu2e::CaloExample);


/*
	   if (!contentMC.hasConversion() && clusterIt->energyDep() > 70)
	   {
	         std::cout<<"Cluster content "<<std::endl;

        	 for (const auto& contentMap : contentMC.simContentMap() )
		 {	       
		   art::Ptr<SimParticle> sim = contentMap.first;
        	   CaloContentSim       data = contentMap.second;
		   
		   std::cout<<sim->id().asInt()<<" "<<sim->startPosition()<<" "<<sim->endPosition()<<"   "<<data.time()<<" "<<data.edep()<<std::endl;
		 }		 
		 
		 art::Handle<GenParticleCollection> gensHandle;
		 event.getByLabel("generate", gensHandle);
		 GenParticleCollection const& genParticles(*gensHandle);

		 //art::Handle<SimParticleCollection> simsHandle;
		 //event.getByLabel("g4run", simsHandle);
		 //SimParticleCollection const& simParticles2(*simsHandle);
		 
		 std::cout<<"Generated"<<std::endl;
		 for (const auto& gen : genParticles) std::cout<<gen.generatorId().name()<<" "<<gen.pdgId()<<" "<<gen.momentum().rho()<<"   "<<gen.position()<<std::endl;


		 //std::cout<<"Simulated"<<std::endl;
		 //for (const auto& sim2 : simParticles2)
		 //{
		 //  const SimParticle& simm(sim2.second);
		 //  std::cout<<simm.id().asInt()<<" "<<simm.pdgId()<<" "<<simm.creationCode()<<"   "<<simm.startMomentum().mag()<<" "<<simm.startPosition()<<" / "<<simm.endPosition()<<"   ";
                 //  if (simm.hasParent()) std::cout<<"p="<<simm.parent()->id().asInt()<<"   ";
		 //  if (simm.genParticle()) std::cout<<"g="<<simm.genParticle()->generatorId().id();
		   
		 //  std::cout<<std::endl;
		 //}    
	   
	   }

*/
