 //
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
// $Id: CaloExample_module.cc,v 1.4 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Original author Bertrand Echenard
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
#include "CaloCluster/inc/CaloContentMC.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/HelixTraj.hh"

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

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

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>




using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;



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
       
       typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
       typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;
       typedef art::Ptr<SimParticle> SimParticlePtr;


       int _diagLevel;
       int _nProcess;

       std::string _g4ModuleLabel;
       std::string _generatorModuleLabel;
       std::string _caloReadoutModuleLabel;
       std::string _caloCrystalModuleLabel;
       std::string _caloHitMCCrystalPtrLabel;
       std::string _caloClusterModuleLabel;
       std::string _caloClusterAlgorithm;
       std::string _caloClusterSeeding;
       const std::string _producerName;
       std::string _virtualDetectorLabel;
       std::string _trkPatRecModuleLabel;
       std::string _instanceName;

       TrkParticle _tpart;
       TrkFitDirection _fdir;

       
       
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
       float _cluEnergy[16384],_cluTime[16384],_cluCogX[16384],_cluCogY[16384],_cluCogZ[16384],_cluDist[16384];
       int   _cluDad[16384],_cluConv[16384],_cluSimIdx[16384],_cluSimLen[16384];
       std::vector<std::vector<int> > _cluList;	 

       int   _clusimId[16384],_clusimPdgId[16384],_clusimGenIdx[16384];
       float _clusimMom[16384],_clusimPosX[16384],_clusimPosY[16384],_clusimPosZ[16384],_clusimTime[16384],_clusimEdep[16384];

       int   _nVd,_vdId[16384],_vdPdgId[16384];
       float _vdTime[16384],_vdPosX[16384],_vdPosY[16384],_vdPosZ[16384],_vdMom[16384];

       int _nTrkOk,_nTrk,_trkOk[8192],_trkstat[8192],_trknHit[8192];
       float _trkDip[8192],_trkpt[8192],_trkcon[8192],_trkmomErr[8192];



       

  };


  CaloExample::CaloExample(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel","g4run")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel","generate")),
    _caloReadoutModuleLabel(pset.get<string>("caloReadoutModuleLabel","CaloReadoutHitsMaker")),
    _caloCrystalModuleLabel(pset.get<string>("caloCrystalModuleLabel","CaloCrystalHitsMaker")),
    _caloHitMCCrystalPtrLabel(pset.get<string>("calorimeterHitMCCrystalPtr","CaloHitMCCrystalPtr")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel","CaloClusterMakerNew")),
    _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
    _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
    _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
    _virtualDetectorLabel(pset.get<string>("virtualDetectorName","virtualdetector")),
    _trkPatRecModuleLabel(pset.get<string>("trkPatRecModuleLabel")),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _Ntup(0)

  {
    _instanceName = _fdir.name() + _tpart.name();
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
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomX[nGen]/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomX[nGen]/F");
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
    _Ntup->Branch("cluDist",      &_cluDist ,     "cluDist[nCluster]/F");	
    _Ntup->Branch("cluDad",       &_cluDad ,      "cluDad[nCluster]/I");	
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

    _Ntup->Branch("nTrkOk",       &_nTrkOk ,      "nTrkOk/I");
    _Ntup->Branch("nTrk",         &_nTrk ,        "nTrk/I");
    _Ntup->Branch("trkDip",       &_trkDip ,      "trkDip[nTrk]/F");
    _Ntup->Branch("trkOk",        &_trkOk ,       "trkOk[nTrk]/I");
    _Ntup->Branch("trkpt",        &_trkpt ,       "trkpt[nTrk]/F");
    _Ntup->Branch("trkstat",      &_trkstat ,     "trkstat[nTrk]/I");
    _Ntup->Branch("trkcon",       &_trkcon ,      "trkcon[nTrk]/F");
    _Ntup->Branch("trknHit",      &_trknHit ,     "trknHit[nTrk]/I");
    _Ntup->Branch("trkmomErr",    &_trkmomErr ,   "trkmomErr[nTrk]/F");

   

  }

	

  void CaloExample::endJob(){
  }




  void CaloExample::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%10==0) std::cout<<"Processing event "<<_nProcess<<std::endl;
      
      
      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  

      //Get generated particles
      art::Handle<GenParticleCollection> gensHandle;
      event.getByLabel(_generatorModuleLabel, gensHandle);
      GenParticleCollection const& genParticles(*gensHandle);

      //Get simulated particles by Geant4
      art::Handle<SimParticleCollection> simsHandle;
      event.getByLabel(_g4ModuleLabel, simsHandle);            
      SimParticleCollection const& simParticles(*simsHandle);

      
      
      //Get calorimeter readout hits (2 readout / crystal as of today)
      art::Handle<CaloHitCollection> caloHitsHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitsHandle);
      CaloHitCollection const& caloHits(*caloHitsHandle);
      
      //Get calorimeter readout hits MC level - energy/time/type
      art::Handle<CaloHitMCTruthCollection> caloHitMCTruthHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitMCTruthHandle);
      CaloHitMCTruthCollection const& caloHitsMCTruth(*caloHitMCTruthHandle);
      
      //Get stepPointMC for crystal readout hits
      art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
      event.getByLabel(_caloReadoutModuleLabel,_caloHitMCCrystalPtrLabel,mcptrHandle);
      PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
      
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
      event.getByLabel(_caloClusterModuleLabel,_producerName, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);

      //Get virtual detector hits
      art::Handle<StepPointMCCollection> vdhits;
      event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);

      //Get tracks
      art::Handle<KalRepCollection> krepsHandle;
      event.getByLabel(_trkPatRecModuleLabel,_instanceName,krepsHandle);
      KalRepCollection const& kreps = *krepsHandle;

      
      //Utility to match  cloHits with MCtruth, simParticles and StepPoints
      CaloHitMCNavigator caloHitNavigator(caloHits, caloHitsMCTruth, *hits_mcptr, caloHitSimPartMC);


      const double CrDensity = 4.9*(CLHEP::g/CLHEP::cm3);
      const double CrMass    = CrDensity*cal.caloGeomInfo().crystalVolume();
	   







       //--------------------------  Do generated particles --------------------------------

  
       _evt = event.id().event();
       _run = event.run();
              
       _nGen = genParticles.size();
       for (unsigned int i=0; i < genParticles.size(); ++i)
       {
           GenParticle const* gen = &genParticles[i];
	   _genPdgId[i]   = gen->pdgId();
	   _genCrCode[i]  = gen->generatorId().id();
	   _genmomX[i]    = gen->momentum().vect().x() + 3904;
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
	   CLHEP::Hep3Vector crystalPos = cal.crystalOrigin(hit.id());
           
           CaloHit const& caloHit         = *(hit.readouts().at(0));
	   CaloHitSimPartMC const& hitSim = caloHitNavigator.sim(caloHit);
           int nPartInside                = hitSim.simParticles().size();

	   _cryEtot             += hit.energyDep();
	   _cryTime[_nHits]      = hit.time();
	   _cryEdep[_nHits]      = hit.energyDep();
	   _cryDose[_nHits]      = hit.energyDep() / CrMass / (CLHEP::joule/CLHEP::kg); //dose
	   _cryPosX[_nHits]      = crystalPos.x()+ 3904;
	   _cryPosY[_nHits]      = crystalPos.y();
	   _cryPosZ[_nHits]      = crystalPos.z()-10200;
	   _cryId[_nHits]        = hit.id();
	   _crySectionId[_nHits] = cal.caloSectionId(hit.id());
           _crySimIdx[_nHits]    = _nSim;   	              
           _crySimLen[_nHits]    = nPartInside;


	   for (int ip=0; ip<nPartInside;++ip) 
	   { 
	   
	     art::Ptr<SimParticle> const& mother = hitSim.simParticles().at(ip);	   		
             art::Ptr<StepPointMC> const& mchit  = hitSim.stepPoints().at(ip);
 	     
	     art::Ptr<SimParticle> grandMother = mother;
             while (grandMother->hasParent()) grandMother = grandMother->parent();
	     GenParticle const* generated = grandMother->genParticle() ? grandMother->genParticle().get() : 0;
             
	     _motId[_nSim]      = mother->id().asInt();
	     _motPdgId[_nSim]   = mother->pdgId();
	     _motmom[_nSim]     = mchit->momentum().mag();
	     _motcrCode[_nSim]  = mother->creationCode();
	     _motStartX[_nSim]  = mother->startPosition().x()+ 3904.;
	     _motStartY[_nSim]  = mother->startPosition().y();
	     _motStartZ[_nSim]  = mother->startPosition().z() - 10200;
	     _motStartT[_nSim]  = mother->startGlobalTime();
	     _motPosX[_nSim]    = mchit->position().x() + 3904.;  //value used to shift in tracker coordinate system
	     _motPosY[_nSim]    = mchit->position().y();
	     _motPosZ[_nSim]    = mchit->position().z() - 10200;  //value used to shift in tracker coordinate system
	     _motTime[_nSim]    = mchit->time();
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
	  std::vector<art::Ptr<StepPointMC> > const& step = clutil.stepPoints();
	  std::vector<double>                 const& edep = clutil.edepTot();

	  std::vector<int> _list;
	  for (unsigned int i=0;i<clusterIt->size();++i) 
	  {
	    int idx = int(clusterIt->caloCrystalHitsPtrVector().at(i).get()- &caloCrystalHits.at(0));
	    _list.push_back(idx);	    
          }
	  	  
	  _cluEnergy[_nCluster] = clusterIt->energyDep();
	  _cluTime[_nCluster]   = clusterIt->time();
	  _cluNcrys[_nCluster]  = clusterIt->size();
          _cluCogX[_nCluster]   = clusterIt->cog3Vector().x()+ 3904.;
          _cluCogY[_nCluster]   = clusterIt->cog3Vector().y();
          _cluCogZ[_nCluster]   = clusterIt->cog3Vector().z();
	  _cluDist[_nCluster]   = clusterIt->distance();
	  _cluDad[_nCluster]    = clusterIt->daddy();	  
	  _cluConv[_nCluster]   = clutil.hasConversion();	  
          _cluSimIdx[_nCluster] = _nCluSim;   	              
	  _cluSimLen[_nCluster] = sim1.size();
	  _cluList.push_back(_list);
	  	  
	  //note: a signal cluster is defined as _cluConv==1 and _cluSimLen==1 (if it has a signal and there is only one inside)
	  
	  for (unsigned int ip=0; ip<sim1.size();++ip) 
	  {
              art::Ptr<SimParticle> smother = sim1[ip];
              while (smother->hasParent()) smother = smother->parent();
	      int genIdx=-1;
	      if (smother->genParticle()) genIdx = smother->genParticle()->generatorId().id();
            
	     _clusimId[_nCluSim]     = sim1[ip]->id().asInt();
             _clusimPdgId[_nCluSim]  = sim1[ip]->pdgId();
             _clusimGenIdx[_nCluSim] = genIdx;
             _clusimMom[_nCluSim]    = step[ip]->momentum().mag();
 	     _clusimPosX[_nCluSim]   = step[ip]->position().x() + 3904.;  //value used to shift in tracker coordinate system
	     _clusimPosY[_nCluSim]   = step[ip]->position().y();
	     _clusimPosZ[_nCluSim]   = step[ip]->position().z() - 10200;  //value used to shift in tracker coordinate system
	     _clusimTime[_nCluSim]   = step[ip]->time();
	     _clusimEdep[_nCluSim]   = edep[ip];

   	     ++_nCluSim;
	  }

	  ++_nCluster;
       }


       
       //--------------------------  Do virtual detectors --------------------------------
       //73/74/77/78 front back inner outer edges disk 0
       //75/76/79/80 front back inner outer edges disk 1
       
       for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter)
       {	   
	  const StepPointMC& hit = *iter;

	  if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;

	  _vdId[_nVd]    = hit.volumeId();
	  _vdPdgId[_nVd] = hit.simParticle()->pdgId();
	  _vdTime[_nVd]  = hit.time();
	  _vdPosX[_nVd]  = hit.position().x()+ 3904;
	  _vdPosY[_nVd]  = hit.position().y();
	  _vdPosZ[_nVd]  = hit.position().z()-10200;
	  _vdMom[_nVd]   = hit.momentum().mag();
          ++_nVd;	   
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
         double fitCon       = krep.chisqConsistency().significanceLevel();

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
	  ++_nTrk;

	   if (cutC)  ++_nTrkOk;
        }


  	_Ntup->Fill();
  




  }



}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloExample);


