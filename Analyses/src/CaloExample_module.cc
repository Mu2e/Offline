//
// An EDAnalyzer module that reads back the hits created by G4 and makes histograms.
//
// $Id: CaloExample_module.cc,v 1.1 2013/03/08 01:22:31 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:31 $
//
// Original author Bertrand Echenard

#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"

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
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"

#include "CaloCluster/inc/CaloClusterer.hh"
#include "CaloCluster/inc/CaloContentMC.hh"

// Mu2e includes.
#include "KalmanTests/inc/KalRepCollection.hh"

// BaBar Kalman filter includes
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/HelixTraj.hh"


#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>



using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;




  // Anonymous namespace to hold some helper classes.
  namespace {

       class CaloMCCrystalInfo {

           public:

             int     _crystalId;
	     mu2e::SimParticle const* _sim;
             mu2e::StepPointMC const* _step;                

             CaloMCCrystalInfo(int id, mu2e::SimParticle const* sim, mu2e::StepPointMC const* step): 
	       _crystalId(id), _sim(sim), _step(step) {}
       };
  }



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
       std::string _g4ModuleLabel;
       std::string _generatorModuleLabel;
       std::string _caloReadoutModuleLabel;
       std::string _caloCrystalModuleLabel;
       std::string _stepPoints;

       std::string _caloClusterModuleLabel;
       string _caloClusterAlgorithm;
       string _caloClusterSeeding;
       const string _producerName;
       string _trkPatRecModuleLabel;
       string _instanceName;

       TTree* _Ntup;
       int _nProcess;

       int _evt,_run,_nHits,_nSim,_nCluster,_nCluSim,_nGen;
       int _genPdgId[16384],_genCrCode[16384];
       int _cryId[16384],_crySectionId[16384],_crySimIdx[16384],_crySimLen[16384];
       int _motId[16384],_motPdgId[16384],_motcrCode[16384],_motGenIdx[16384];
       int _cluNcrys[128];

       float _genmomX[16384],_genmomY[16384],_genmomZ[16384],_genStartX[16384],_genStartY[16384],_genStartZ[16384],_genStartT[16384];

       float _cryEtot,_cryTime[16384],_cryEdep[16384],_cryDose[16384],_cryPosX[16384],_cryPosY[16384],_cryPosZ[16384];
       
       float _motmom[16384],_motStartX[16384],_motStartY[16384],_motStartZ[16384],_motStartT[16384];
       float _motTime[16384],_motEdep[16348],_motPosX[16384],_motPosY[16384],_motPosZ[16384];

       float _cluEnergy[2048],_cluTime[2048],_cluCogX[2048],_cluCogY[2048],_cluCogZ[2048],_cluDist[2048];
       int   _cluDad[2048],_cluConv[2048],_cluSimIdx[2048],_cluSimLen[2048];
       std::vector<std::vector<int> > _cluList;	 

       int   _clusimId[2048],_clusimPdgId[2048],_clusimGenIdx[2048];
       float _clusimMom[2048],_clusimPosX[2048],_clusimPosY[2048],_clusimPosZ[2048],_clusimTime[2048],_clusimEdep[2048];

       int _nTrkOk,_nTrk,_trkOk[8192];
       float _trkDip[8192],_trkpt[8192]; 
  };


  CaloExample::CaloExample(fhicl::ParameterSet const& pset) :
    _diagLevel(pset.get<int>("diagLevel",0)),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _caloReadoutModuleLabel(pset.get<string>("caloReadoutModuleLabel","CaloReadoutHitsMaker")),
    _caloCrystalModuleLabel(pset.get<string>("caloCrystalModuleLabel","CaloCrystalHitsMaker")),
    _stepPoints(pset.get<string>("calorimeterStepPoints","calorimeter")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel","CaloClusterMakerNew")),
    _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
    _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
    _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
    _trkPatRecModuleLabel(pset.get<string>("trkPatRecModuleLabel")),
    _instanceName(pset.get<string>("instanceName")),
    _Ntup(0),_nProcess(0)

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
    _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosZ[nCry]/F");
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



    _Ntup->Branch("nTrkOk",       &_nTrkOk ,      "nTrkOk/I");
    _Ntup->Branch("nTrk",         &_nTrk ,        "nTrk/I");
    _Ntup->Branch("trkDip",       &_trkDip ,      "trkDip[nTrk]/F");
    _Ntup->Branch("trkOk",        &_trkOk ,       "trkOk[nTrk]/I");
    _Ntup->Branch("trkpt",        &_trkpt ,       "trkpt[nTrk]/F");

  }



  void CaloExample::endJob(){
  }




  void CaloExample::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%1==0) std::cout<<"Processing event "<<_nProcess<<std::endl;
      //Get handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      const double CrDensity = 7.4*(CLHEP::g/CLHEP::cm3);
      const double CrMass    = CrDensity*cal.crystalVolume();



      art::Handle<GenParticleCollection> gensHandle;
      event.getByLabel(_generatorModuleLabel, gensHandle);
      GenParticleCollection const& genParticles(*gensHandle);

      art::Handle<SimParticleCollection> simsHandle;
      event.getByLabel(_g4ModuleLabel, simsHandle);            //"neutronMixer"
      SimParticleCollection const& simParticles(*simsHandle);

      art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
      event.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);
      PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();

      //data about hits in the calorimeter crystals
      art::Handle<CaloHitCollection> caloHitsHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitsHandle);
      CaloHitCollection const& caloHits(*caloHitsHandle);

      art::Handle<CaloHitMCTruthCollection> caloHitMCTruthHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitMCTruthHandle);
      CaloHitMCTruthCollection const& caloHitsMCTruth(*caloHitMCTruthHandle);

      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
      CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);


      //data about clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(_caloClusterModuleLabel,_producerName, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);

      art::Handle<CaloHitSimPartMCCollection> caloHitSimMCHandle;
      event.getByLabel(_caloReadoutModuleLabel, caloHitSimMCHandle);
      CaloHitSimPartMCCollection const& caloHitSimPartMC(*caloHitSimMCHandle);
      
      CaloHitMCNavigator caloHitNavigator(caloHits, caloHitsMCTruth, *hits_mcptr, caloHitSimPartMC);






/*
      for (GenParticleCollection::const_iterator it = genParticles.begin(); it != genParticles.end(); ++it) 
        std::cout<<it->pdgId()<<"  "<<it->generatorId()<<" "<<it->time()<<" "<<it->position()<<std::endl;    

      std::cout<<"Bertrand sim part "<<std::endl;
      for (SimParticleCollection::const_iterator it = simParticles.begin(); it != simParticles.end(); ++it) {
        SimParticle const& sim = it->second;
        //std::cout<<sim.id()<<" "<<sim.pdgId()<<"  "<<sim.startGlobalTime()<<"  "<<sim.startPosition()<<" "<<sim.endPosition()<<" "<<sim.startVolumeIndex()<<"    ";   
        std::cout<<sim.id()<<" "<<sim.pdgId()<<"  "<<sim.startVolumeIndex()<<" "<<sim.endVolumeIndex()<<" "<<sim.startPosition()<<" "<<sim.endPosition()<<" "<<
	sim.startGlobalTime()<<"  "<<sim.endGlobalTime()<<"  "<<sim.fromGenerator()<<"  "<<sim.startMomentum().vect().mag()<<"   ";   
        if (sim.hasParent()) std::cout<<sim.parentId();
        std::cout<<std::endl;  
      }
*/

 
 
 
 
 
 
  
       _evt = event.id().event();
       _run = event.run();
              
       _nGen = genParticles.size();
       for (unsigned int i=0; i < genParticles.size(); ++i){
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


       _nCluster = _nCluSim = 0;
       _cluList.clear();
       for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt) 
       {

	  CaloContentMC clutil(caloHitNavigator,*clusterIt);
	  std::vector<art::Ptr<SimParticle> > const& sim1 = clutil.simPart();
	  std::vector<art::Ptr<StepPointMC> > const& step = clutil.stepPoints();
	  std::vector<double >                const& edep = clutil.edepTot();

	  std::vector<int> _list;
	  for (unsigned int i=0;i<clusterIt->size();++i) 
	  {
	    size_t idx = clusterIt->caloCrystalHitsPtrVector().at(i).get()- &caloCrystalHits.at(0);
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

 //std::cout<<"ip="<<ip<<"   "<<sim1[ip]->id().asInt()<<" "<<sim1[ip]->pdgId()<<" "<<step[ip]->momentum().mag()
 //<<" "<<_clusimPosX[_nCluSim]<<" "<<_clusimPosY[_nCluSim]<<" "<<_clusimPosZ[_nCluSim]<<" "<<step[ip]->time()<<" "<<edep[ip]<<std::endl;  
  	     ++_nCluSim;
	  }



	  ++_nCluster;
       }



       art::Handle<KalRepCollection> krepsHandle;
       event.getByLabel(_trkPatRecModuleLabel,_instanceName,krepsHandle);
       KalRepCollection const& kreps = *krepsHandle;

       _nTrkOk = 0;
       _nTrk = 0;
       for ( KalRepCollection::const_iterator i=kreps.begin(); i != kreps.end(); ++i ){

         // Reference to one fitted track.
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
         bool cutC = ( krep.fitStatus().success() == 1) &&
           ( fitCon             > 1.e-3  ) &&
           ( krep.nActive()    >= 25    ) &&
           ( fitMomErr          < 0.18   );

	  _trkDip[_nTrk] = tanDip;
	  _trkOk[_nTrk]  = cutC? 1 : 0;
	  _trkpt[_nTrk]  = fitmompt;
	  ++_nTrk;

	   if (cutC)  ++_nTrkOk;
           //if (cutC && fitmompt > 54) ++_nTrkOk;

        }






  	_Ntup->Fill();
  




  }



}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloExample);




















/*  //in many function
    
    art::ProductInstanceNameSelector getCrystalSteps(_stepPoints);
    HandleVector crystalStepsHandles;
    event.getMany( getCrystalSteps, crystalStepsHandles);

    HitMap hitmapCrystal;        
    fillMapById( hitmapCrystal, crystalStepsHandles);

    for(HitMap::const_iterator crIter = hitmapCrystal.begin(); crIter != hitmapCrystal.end(); ++crIter ) {
      
	int crid = crIter->first;
	StepPtrs const& isteps = crIter->second;

        unsigned int minHit(1000000000);
	for ( StepPtrs::const_iterator i=isteps.begin(), e=isteps.end(); i != e; ++i ){

              StepPointMC const& h = **i;

              if ( h.eDep()<=0.0 ) continue;
              if (h.simParticle()->id().asInt()<minHit) minHit = h.simParticle()->id().asInt();

              CLHEP::Hep3Vector pos0 = h.position();
              SimParticle const* part(h.simParticle().get());
              while (part->hasParent()) {
                   _hdist->Fill((part->parent()->startPosition()-pos0).mag()); 
		   _htime->Fill(h.simParticle()->startGlobalTime()-part->startGlobalTime());
		   part = part->parent().get();
              }
   	      
              double dist(0);
              CLHEP::Hep3Vector posParticle = h.position();
              SimParticle const* mother = h.simParticle().get();
              while (mother->hasParent() && dist < 200) {
                   dist = (mother->parent()->startPosition()-posParticle).mag();
		   mother = mother->parent().get();
              }	      
	 }

     }



  void fillMapById(HitMap& map,HandleVector const& crystalStepsHandles);

  void CaloExample::fillMapById( HitMap& hitmap,HandleVector const& crystalStepsHandles)
  {
      for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end(); i != e; ++i ){

	art::Handle<StepPointMCCollection> const& handle(*i);
	StepPointMCCollection const& steps(*handle);

	StepPointMCCollection::const_iterator j0=steps.begin();
	for ( StepPointMCCollection::const_iterator j=j0, je=steps.end(); j != je; ++j ){
          StepPointMC const& step(*j);
          hitmap[step.volumeId()].push_back( StepPtr(handle,j-j0) );
	}
      }

  } 
*/



/*
       int nclus(0);
       for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt) {

	     std::map<unsigned int, double> clutserContent;

	     for (unsigned int i=0;i<clusterIt->size();++i) {

		  int cid = clusterIt->caloCrystalHitsPtrVector().at(i)->id();

		  std::map<unsigned int,unsigned int>::const_iterator toGenerated = caloHitToGenMother.find(cid);
		  if (toGenerated==caloHitToGenMother.end()) continue;
		  unsigned int generated = toGenerated->second;
		  (clutserContent[generated]) += clusterIt->caloCrystalHitsPtrVector().at(i)->energyDep();	      
	     }   

	     //std::cout<<"For cluster "<<nclus<<" we have as parent / true cog position"<<std::endl;

	     for (std::map<unsigned int, double>::const_iterator it=clutserContent.begin();it!=clutserContent.end();++it) {

	       double minTime(1e6);
	       Hep3Vector cogPosition(0,0,0); 
	       for (unsigned int i=0; i<caloCrystalHits.size();++i) {

        	  PtrStepPointMCVector const& mcptr(hits_mcptr->at(i));      
		  for ( PtrStepPointMCVector::const_iterator j=mcptr.begin(), je=mcptr.end(); j != je; ++j ){
        	     StepPointMC const& step(**j);
 		     if (!step.simParticle()->hasParent() || step.simParticle()->parent()->id().asInt() != it->first) continue;
		     if (step.time() > minTime) continue;
		     minTime = step.time();
		     cogPosition = step.position();
		  }

		}

		//std::cout<<"True cog "<<cogPosition<<std::endl;	      
	     }   

	     ++nclus;	   
       }
*/

/*
 	  // Remove SimParticles having their mother in the same set (typical secondary interactions)
 	  if (simParents.size()>1) {
	     std::set<SimParticle const* >::const_iterator it = simParents.begin(),tmp;
	     while ( it != simParents.end() ) {
	         tmp = it++;
	         if ( (*tmp)->hasParent() && simParents.find((*tmp)->parent().get()) != simParents.end()) simParents.erase(tmp);	     
	     }
	  }
*/




/*
       // Find the SimParticle(s) associated to each hits and the momentum of the incident particle at the calorimeter boundary
       //  - for each StepPointMC, find mother by going up the tree of SimParticles until one starts outside the calo (jums between calo sections are ok) 
       //  - record the momentum for the mother particles since momentum is not recorded in SimParticles at the calorimeter boundary
       //    so must use StePointMC to estimate the momentum at the entrance of the calo
 
      
       std::cout<<"Number of crystals = "<<caloCrystalHits.size()<<std::endl;
       
       std::vector<CaloMCCrystalInfo> hitMCInfo;
       
       
       for (unsigned int ic=0; ic<caloCrystalHits.size();++ic) {

	   CaloCrystalHit const& hit = caloCrystalHits.at(ic);
	   PtrStepPointMCVector const& mcptr( hits_mcptr->at(cal.getROBaseByCrystal(ic)) );


	   std::map<int, SimParticle const* >    simParents;
	   std::map<SimParticle const* ,StepPointMC const*> momentumParents;

	   for (unsigned int i=0; i<mcptr.size(); ++i) 
	   {
 //if (mcptr[i]->simParticle()->parent()) std::cout<<mcptr[i]->simParticle()->id().asInt()<<" "<<mcptr[i]->simParticle()->parent()->id().asInt()<<" "<<mcptr[i]->simParticle()->startMomentum().rho()<<" "<<mcptr[i]->time()<<" "<<mcptr[i]->momentum().mag()<<std::endl;
 //else std::cout<<mcptr[i]->simParticle()->id().asInt()<<" "<<0<<" "<<mcptr[i]->simParticle()->startMomentum().rho()<<" "<<mcptr[i]->time()<<" "<<mcptr[i]->momentum().mag()<<std::endl;	     
 //std::cout<<"-"<<std::endl;

	       StepPointMC const* mchit = mcptr[i].get();
	       int thisSimpartId        = mchit->simParticle()->id().asInt();
	       int parentSimpartId      = (mchit->simParticle()->parent()) ? mchit->simParticle()->parent()->id().asInt() : -1;


	       //already in the list
	       std::map<int, SimParticle const* >::iterator ifind = simParents.find(thisSimpartId); 
	       if ( ifind != simParents.end() )
	       {
		   SimParticle const* mommy = ifind->second;
		   std::map<SimParticle const* , StepPointMC const*>::iterator mfind = momentumParents.find(mommy);		 
		   if (mfind != momentumParents.end() && mchit->momentum().mag() > mfind->second->momentum().mag()) mfind->second = mchit; 
		   continue;		 
	       }	    

               //mother already in the list
	       std::map<int, SimParticle const* >::iterator ifindParent = simParents.find(parentSimpartId);              
	       if ( ifindParent != simParents.end() )
	       {
		   SimParticle const* mommy = ifindParent->second;
		   simParents.insert( std::pair<int,SimParticle const*>(thisSimpartId, mommy) );

		   std::map<SimParticle const* , StepPointMC const*>::iterator mfind = momentumParents.find(mommy);		 
		   if (mfind!=momentumParents.end() && mchit->momentum().mag() > mfind->second->momentum().mag()) mfind->second = mchit; 
		   continue;

	       }	    

	       //must be added to the list
	       SimParticle const* mother = mchit->simParticle().get();
               while (mother->hasParent() && !cal.hasCrossedCalo(mother->startPosition(),mother->endPosition()) ) mother = mother->parent().get();

	       simParents.insert( std::pair<int, SimParticle const*>(thisSimpartId, mother) );

   	       std::map<SimParticle const* , StepPointMC const*>::iterator mfind = momentumParents.find(mother);		 
	       if  (mfind!=momentumParents.end() && mchit->momentum().mag() > mfind->second->momentum().mag()) mfind->second = mchit; 
	       else momentumParents.insert(std::pair<SimParticle const* , StepPointMC const*>(mother, mchit));

	    }   

            //fill information for the cluster, keep results of each crystal
            for (std::map<SimParticle const* , StepPointMC const* >::const_iterator it=momentumParents.begin();it!=momentumParents.end();++it)
               hitMCInfo.push_back(CaloMCCrystalInfo(hit.id(),it->first, it->second));


//if (momentumParents.size()>1) {
// std::cout<<"Refined mother list ";
// for (std::map<int, SimParticle const*  >::const_iterator it=simParents.begin();it!=simParents.end();++it) std::cout<<it->first<<" ";
// std::cout<<std::endl; 
// std::cout<<"Simparticle momentum ";
// for (std::map<SimParticle const* , double >::const_iterator it=momentumParents.begin();it!=momentumParents.end();++it) std::cout<<it->first->id().asInt()<<" "<<it->second<<"     ";
// std::cout<<std::endl; 
//}


       }


       
       std::map<SimParticle const*, StepPointMC const* > clusterMomentumContent;
       std::set<SimParticle const* > clusterSimpartBase;
       std::set<GenParticle const* > clusterGenpartBase;


       for (std::vector<CaloMCCrystalInfo>::const_iterator it=hitMCInfo.begin(); it!=hitMCInfo.end();++it)
       {         
	  std::map<SimParticle const*, StepPointMC const* >::iterator ifind = clusterMomentumContent.find(it->_sim);		 
          if (ifind!=clusterMomentumContent.end() && it->_step->momentum().mag() > ifind->second->momentum().mag()) ifind->second = it->_step;
	  else clusterMomentumContent.insert(std::pair<SimParticle const* , StepPointMC const*>( it->_sim, it->_step));       
       
 	  SimParticle const* mother = it->_sim;
          while (mother->hasParent()) mother = mother->parent().get();
         
	  clusterSimpartBase.insert(mother);
	  if (mother->genParticle())  clusterGenpartBase.insert(mother->genParticle().get());
       
       }









bool mustReprint(false);
for (std::vector<CaloMCCrystalInfo>::const_iterator it=hitMCInfo.begin(); it!=hitMCInfo.end();++it)
{       
  std::cout<<it->_crystalId<<" "<<it->_step->momentum().mag()<<" "<<it->_sim->id().asInt()<<" "<<it->_sim->pdgId()<<" "<<it->_sim->hasParent()<<"   "<<it->_sim->creationCode()<<"  "<<it->_sim->fromGenerator();
  if (it->_sim->genParticle())  std::cout<<it->_sim->genParticle()->generatorId();  
  std::cout<<std::endl;  
//  if (it->_sim->id().asInt()>2 )mustReprint=true;
}

  std::cout<<"Summary ";
  std::cout<<clusterSimpartBase.size()<<" "<<clusterGenpartBase.size()<<endl;

  for (std::set<SimParticle const* >::const_iterator it=clusterSimpartBase.begin(); it!=clusterSimpartBase.end();++it)
   std::cout<<(*it)->id().asInt()<<" "<<(*it)->pdgId()<<" "<<(*it)->genParticle()<<std::endl;
  
  for (std::set<GenParticle const* >::const_iterator it=clusterGenpartBase.begin(); it!=clusterGenpartBase.end();++it)
   std::cout<<(*it)->generatorId()<<endl;

  std::cout<<"-----";









if (clusterSimpartBase.size() > 1){

      for (GenParticleCollection::const_iterator it = genParticles.begin(); it != genParticles.end(); ++it) 
        std::cout<<it->pdgId()<<"  "<<it->generatorId()<<" "<<it->time()<<" "<<it->position()<<std::endl;    


      std::cout<<"Bertrand sim part "<<std::endl;
      for (SimParticleCollection::const_iterator it = simParticles.begin(); it != simParticles.end(); ++it) {
        SimParticle const& sim = it->second;
        //std::cout<<sim.id()<<" "<<sim.pdgId()<<"  "<<sim.startGlobalTime()<<"  "<<sim.startPosition()<<" "<<sim.endPosition()<<" "<<sim.startVolumeIndex()<<"    ";   
        std::cout<<sim.id()<<" "<<sim.pdgId()<<"  "<<sim.startVolumeIndex()<<" "<<sim.endVolumeIndex()<<" "<<sim.startPosition()<<" "<<sim.endPosition()<<" "<<sim.startGlobalTime()<<"  "<<sim.endGlobalTime()<<"  "<<sim.fromGenerator()<<"   ";   
        if (sim.hasParent()) std::cout<<sim.parentId();
        std::cout<<std::endl;  

      }
      std::cout<<"--"<<std::endl;
}


*/


/*
double motherMom(-1);
CLHEP::Hep3Vector momMot1InSection(0,0,0);
CLHEP::Hep3Vector momMot2InSection(0,0,0);
CLHEP::Hep3Vector posMotInSection(0,0,0);
if (mother) motherMom=mother->startMomentum().vect().mag();

 
std::map<SimParticle const*, StepPointMC const*>::iterator ifind = clusterStepPointContent.find(mother);		 

if (it != clusterStepPointContent.end()) {
 motherMom = it->second->momentum().mag();
 momMot1InSection = cal.toSectionFrame(cal.getCaloSectionId(hit.id()),it->second->position());
 momMot2InSection = cal.toSectionFrame(cal.getCaloSectionId(hit.id()),it->second->position()+it->second->momentum());
 posMotInSection = cal.toSectionFrame(cal.getCaloSectionId(hit.id()),it->second->position());
}
CLHEP::Hep3Vector momMotInSection = momMot2InSection-momMot1InSection;



	   if (mother) {
	     _motSecPosX[_nHits]  = posMotInSection.x();
	     _motSecPosY[_nHits]  = posMotInSection.y();
	     _motSecPosZ[_nHits]  = posMotInSection.z();
	     _motSecMomX[_nHits]  = momMotInSection.x();
	     _motSecMomY[_nHits]  = momMotInSection.y();
	     _motSecMomZ[_nHits]  = momMotInSection.z();             
	   //cout<<"mot in mu2e "<<it->second->position()<<endl;
	   } else {
	     _motSecPosX[_nHits] = _motSecPosY[_nHits] = _motSecPosZ[_nHits] = 0;
	     _motSecMomX[_nHits] = _motSecMomY[_nHits] = _motSecMomZ[_nHits] = 0;	   
	   }

*/


/*
       //cluster level info
       std::map<SimParticle const*, StepPointMC const* > clusterStepPointContent;
       std::set<SimParticle const* > clusterSimpartBase;
       for (std::vector<CaloMCCrystalInfo>::const_iterator it=hitMCInfo.begin(); it!=hitMCInfo.end();++it)
       {         
	  std::map<SimParticle const*, StepPointMC const*>::iterator ifind = clusterStepPointContent.find(it->_sim);		 
          if (ifind!=clusterStepPointContent.end() && it->_step->momentum().mag() > ifind->second->momentum().mag()) ifind->second = it->_step;
	  else clusterStepPointContent.insert(std::pair<SimParticle const*,StepPointMC const*>( it->_sim, it->_step));       
       
 	  SimParticle const* mother = it->_sim;
          while (mother->hasParent()) mother = mother->parent().get();
  	  clusterSimpartBase.insert(mother);                
       }

*/
/*

if (reprint){


       for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt) 
       {
          std::cout<<"CluEnergy "<<clusterIt->energyDep()<<std::endl;
	  std::cout<<"Cluster content "<<std::endl;
	  for (unsigned int i=0;i<clusterIt->size();++i) 
	    std::cout<<clusterIt->caloCrystalHitsPtrVector().at(i)->id()<<" "; std::cout<<std::endl;

       }



  std::cout<<"Bertrand gen part "<<std::endl;

      for (GenParticleCollection::const_iterator it = genParticles.begin(); it != genParticles.end(); ++it) 
        std::cout<<it->pdgId()<<"  "<<it->generatorId()<<" "<<it->time()<<" "<<it->position()<<std::endl;    

  std::cout<<"Bertrand sim part "<<std::endl;
  for (SimParticleCollection::const_iterator it = simParticles.begin(); it != simParticles.end(); ++it) {
    SimParticle const& sim = it->second;
    //std::cout<<sim.id()<<" "<<sim.pdgId()<<"  "<<sim.startGlobalTime()<<"  "<<sim.startPosition()<<" "<<sim.endPosition()<<" "<<sim.startVolumeIndex()<<"    ";   
    std::cout<<sim.id()<<" "<<sim.pdgId()<<"  "<<sim.startVolumeIndex()<<" "<<sim.endVolumeIndex()<<" "<<sim.startPosition()<<" "<<sim.endPosition()<<" "<<
    sim.startGlobalTime()<<"  "<<sim.endGlobalTime()<<"  "<<sim.fromGenerator()<<"  "<<sim.startMomentum().vect().mag()<<"   ";   
    if (sim.hasParent()) std::cout<<sim.parentId();
    std::cout<<std::endl;  
  }



std::cout<<"Crystal hit inside"<<std::endl;
       for (unsigned int ic=0; ic<caloCrystalHits.size();++ic) 
       {
	   
	   CaloCrystalHit const& hit = caloCrystalHits.at(ic);
           for (unsigned int i=0; i < hit.readouts().size(); ++i) 
           {
              CaloHit const& caloHit = *(hit.readouts().at(i));
              PtrStepPointMCVector const& mcptr = caloHitNavigator.ptrs(caloHit);
               for (unsigned int ii=0; ii<mcptr.size(); ++ii) 
               {
               StepPointMC const* mchit = mcptr[i].get();
               cout<<hit.id()<<" "<<ii<<" "<<mchit->simParticle()->id().asInt()<<" "<<hit.energyDep()<<std::endl;;
	      
              }      
	   }   
       } 




}


*/
