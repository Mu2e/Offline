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


  class CaloMini : public art::EDAnalyzer {

     public:

       typedef art::Ptr<StepPointMC> StepPtr;
       typedef std::vector<StepPtr>  StepPtrs;
       typedef std::map<int,StepPtrs > HitMap;



       explicit CaloMini(fhicl::ParameterSet const& pset);
       virtual ~CaloMini() { }

       virtual void beginJob();
       virtual void endJob();

       // This is called for each event.
       virtual void analyze(const art::Event& e);





     private:

       int _diagLevel;
       int _nProcess;

       std::string _caloCrystalModuleLabel;
       std::string _caloClusterModuleLabel;



       TTree* _Ntup;
       int   _evt,_run;

       int   _nHits,_cryId[163840],_crySectionId[163840];
       float _cryEtot,_cryTime[163840],_cryEdep[163840],_cryDose[163840],_cryPosX[163840],_cryPosY[163840],_cryPosZ[163840],_cryLeak[163840];


       int   _nCluster,_nCluSim,_cluNcrys[16384];
       float _cluEnergy[16384],_cluTime[16384],_cluCogX[16384],_cluCogY[16384],_cluCogZ[16384],_cluE1[16384],_cluE9[16384],_cluE25[16384],_cluSecMom[16384];
       int   _cluSplit[16384],_cluConv[16384];
       std::vector<std::vector<int> > _cluList;


  };


  CaloMini::CaloMini(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _Ntup(0)

  {
  }

  void CaloMini::beginJob(){

       art::ServiceHandle<art::TFileService> tfs;

       _Ntup  = tfs->make<TTree>("Calo", "Calo");


       _Ntup->Branch("evt",          &_evt ,        "evt/I");
       _Ntup->Branch("run",          &_run ,        "run/I");
       _Ntup->Branch("cryEtot",      &_cryEtot ,    "cryEtot/F");

       _Ntup->Branch("nCry",         &_nHits ,       "nCry/I");
       _Ntup->Branch("cryId",        &_cryId ,       "cryId[nCry]/I");
       _Ntup->Branch("crySectionId", &_crySectionId, "crySectionId[nCry]/I");
       _Ntup->Branch("cryPosX",      &_cryPosX ,     "cryPosX[nCry]/F");
       _Ntup->Branch("cryPosY",      &_cryPosY ,     "cryPosY[nCry]/F");
       _Ntup->Branch("cryPosZ",      &_cryPosZ ,     "cryPosZ[nCry]/F");
       _Ntup->Branch("cryEdep",      &_cryEdep ,     "cryEdep[nCry]/F");
       _Ntup->Branch("cryTime",      &_cryTime ,     "cryTime[nCry]/F");
       _Ntup->Branch("cryDose",      &_cryDose ,     "cryDose[nCry]/F");

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
       _Ntup->Branch("cluList",      &_cluList);


  }



  void CaloMini::endJob(){
  }




  void CaloMini::analyze(const art::Event& event) {

      ++_nProcess;
      if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloMini =  "<<_nProcess <<std::endl;


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


      const double CrDensity = 4.9*(CLHEP::g/CLHEP::cm3);
      const double CrMass    = CrDensity*cal.caloInfo().crystalVolume();






       //--------------------------  Do generated particles --------------------------------


       _evt = event.id().event();
       _run = event.run();


       //--------------------------  Do calorimeter hits --------------------------------

       _nHits = 0;
       _cryEtot = 0.0;

       for (unsigned int ic=0; ic<caloCrystalHits.size();++ic)
       {
           const CaloCrystalHit &hit     = caloCrystalHits.at(ic);
	   int diskId                    = cal.crystal(hit.id()).diskId();
           CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.id()).position());  //in disk FF frame
 

           _cryEtot             += hit.energyDep();
           _cryTime[_nHits]      = hit.time();
           _cryEdep[_nHits]      = hit.energyDep();
           _cryDose[_nHits]      = hit.energyDep() / CrMass / (CLHEP::joule/CLHEP::kg); //dose
           _cryPosX[_nHits]      = crystalPos.x();
           _cryPosY[_nHits]      = crystalPos.y();
           _cryPosZ[_nHits]      = crystalPos.z();
           _cryId[_nHits]        = hit.id();
           _crySectionId[_nHits] = diskId;

           ++_nHits;
       }

       //--------------------------  Do clusters --------------------------------
       _nCluster = 0;
       _cluList.clear();
       for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt)
       {
       
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
           _cluList.push_back(_list);

           ++_nCluster;
       }

 
        _Ntup->Fill();





  }



}  

DEFINE_ART_MODULE(mu2e::CaloMini);
