//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
// $Id: CaloTrig_module.cc,v 1.4 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
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
#include "MCDataProducts/inc/StrawDigiMC.hh"

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


  class CaloTrig : public art::EDAnalyzer {

     public:

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;



       explicit CaloTrig(fhicl::ParameterSet const& pset);
       virtual ~CaloTrig() { }

       virtual void beginJob();
       virtual void endJob();

       // This is called for each event.
       virtual void analyze(const art::Event& e);





     private:

       int _diagLevel;
       int _nProcess;

       std::string   _g4ModuleLabel;
       art::InputTag _simParticleTag;
       std::string   _caloCrystalModuleLabel;
       std::string   _caloClusterModuleLabel;
       std::string   _caloClusterTruthModuleLabel;
       std::string   _stepPointMCLabel;
       const std::string _producerName;


       TTree* _Ntup;


       int   _evt,_run;
       float _cluTot,_cryTot;
       int   _nHits,_cryId[163840],_crySectionId[163840];
       float _cryTime[163840],_cryEdep[163840],_cryEdepErr[163840],_cryPosX[163840],_cryPosY[163840],_cryPosZ[163840];
       int   _nCluster,_nCluSim,_cluNcrys[16384];
       float _cluEnergy[16384],_cluTime[16384],_cluCogX[16384],_cluCogY[16384],_cluCogZ[16384],_cluE1[16384],_cluE9[16384],_cluE25[16384],_cluSecMom[16384];
       int   _cluSplit[16384],_cluConv[16384],_cluSimIdx[16384],_cluSimLen[16384];


  };


  CaloTrig::CaloTrig(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _caloClusterTruthModuleLabel(pset.get<std::string>("caloClusterTruthModuleLabel")),
    _stepPointMCLabel(pset.get<std::string>("stepPointMCLabel")),
    _Ntup(0)

  {
  }

  void CaloTrig::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;

       _Ntup  = tfs->make<TTree>("Calo", "Calo");



       _Ntup->Branch("evt",          &_evt ,        "evt/I");
       _Ntup->Branch("run",          &_run ,        "run/I");
       _Ntup->Branch("cryTot",       &_cryTot ,     "cryTot/F");
       _Ntup->Branch("cluTot",       &_cluTot ,     "cluTot/F");


       _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
       _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
       _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
       _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
       _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
       _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
       _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
       _Ntup->Branch("cryEdepErr",   &_cryEdepErr ,  "cryEdepErr[nCry]/F");
       _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");


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

  }



  void CaloTrig::endJob(){
  }




  void CaloTrig::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloTrig =  "<<_nProcess <<std::endl;


      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if( ! geom->hasElement<Calorimeter>() ) return;
      Calorimeter const & cal = *(GeomHandle<Calorimeter>());

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
      CaloClusterCollection const& caloClusters(*caloClustersHandle);

      //Calorimeter crystal truth assignment
      //art::Handle<CaloClusterMCTruthAssns> caloClusterTruthHandle;
      //event.getByLabel(_caloClusterTruthModuleLabel, caloClusterTruthHandle);
      //const CaloClusterMCTruthAssns& caloClusterTruth(*caloClusterTruthHandle);

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
      event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
      CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);


       //--------------------------  Do generated particles --------------------------------


       _evt = event.id().event();
       _run = event.run();

       if (_diagLevel == 3){std::cout << "processing event in calo_example " << _nProcess << " run and event  = " << _run << " " << _evt << std::endl;}

       
       _nHits  = 0;
       _cryTot = 0.0;
       for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
       {
           const CaloCrystalHit &hit     = caloCrystalHits.at(ic);
	   int diskId                    = cal.crystal(hit.id()).diskId();
           CLHEP::Hep3Vector crystalPos  = cal.crystal(hit.id()).localPosition();  //in disk FF frame
          
           _cryTot              += hit.energyDep();
           _cryTime[_nHits]      = hit.time();
           _cryEdep[_nHits]      = hit.energyDep();
           _cryEdepErr[_nHits]   = hit.energyDepErr();
           _cryPosX[_nHits]      = crystalPos.x();
           _cryPosY[_nHits]      = crystalPos.y();
           _cryPosZ[_nHits]      = crystalPos.z();
           _cryId[_nHits]        = hit.id();
           _crySectionId[_nHits] = diskId;

           ++_nHits;
       }


       //--------------------------  Do clusters --------------------------------
       _nCluster = _nCluSim = 0;
       _cluTot = 0;
       for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
       {       
           //ClusterContentMC contentMC(cal, caloClusterTruth, *clusterIt);
           _cluTot += clusterIt->energyDep();
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
           //_cluConv[_nCluster]   = (contentMC.hasConversion() ? 1 : 0);
	  
           ++_nCluster;
       }





        


 
        _Ntup->Fill();





  }



}  

DEFINE_ART_MODULE(mu2e::CaloTrig);


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
