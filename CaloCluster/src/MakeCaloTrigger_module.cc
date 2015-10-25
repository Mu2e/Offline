
// C++ includes.
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>


// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

//calorimeter packages
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/CaloClusterFinder.hh"
#include "CaloCluster/inc/CaloSeedManager.hh"
#include "CaloCluster/inc/CaloClusterAssociator.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

// Other includes.
#include "cetlib/exception.h"
#include "TH1D.h"







namespace mu2e {



class MakeCaloTrigger : public art::EDProducer {


     public:

             typedef std::vector<CaloCrystalHit const*>  CaloCrystalVec;
             typedef std::list<CaloCrystalHit const*>    CaloCrystalList;


             explicit MakeCaloTrigger(fhicl::ParameterSet const& pset) :
             _diagLevel(pset.get<int>("diagLevel",0)),
             _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
             _producerNameMain(pset.get<std::string>("mainClusterCollName")),
             _EminSeed(pset.get<double>("EminSeed")),
             _EnoiseCut(pset.get<double>("EnoiseCut")),
             _ExpandCut(pset.get<double>("ExpandCut")),
             _timeCut(pset.get<double>("timeCut")),
             _deltaTimePlus(pset.get<double>("deltaTimePlus")),
             _deltaTimeMinus(pset.get<double>("deltaTimeMinus")),
             _maxDistClu(pset.get<double>("maxDistClu")),
             _messageCategory("CLUSTER"),
	     _hE(0)
             {
                     produces<CaloProtoClusterCollection>(_producerNameMain);
             }


             virtual ~MakeCaloTrigger() { }
             virtual void beginJob();
             void produce( art::Event& e);



     private:

             int         _diagLevel;
             std::string _caloCrystalModuleLabel;
             std::string _producerNameMain;

             double      _EminSeed;
             double      _EnoiseCut;
             double      _ExpandCut;
             double      _timeCut;
             double      _deltaTimePlus;
             double      _deltaTimeMinus;
             double      _maxDistClu;

             const std::string _messageCategory;
             TH1F *_hE;


             void MakeCaloTriggers(CaloProtoClusterCollection& caloProtoClustersMain, 
                                   art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);



     };



     void MakeCaloTrigger::beginJob(){
      art::ServiceHandle<art::TFileService> tfs;
      _hE  = tfs->make<TH1F>("cluEner","Cluster energy",150,0,150.);

     }




     void MakeCaloTrigger::produce(art::Event& event) {

	 // Check that calorimeter geometry description exists
	 art::ServiceHandle<GeometryService> geom;
	 if( !(geom->hasElement<Calorimeter>()) ) return;

	 //Get handles to calorimeter crystal hits
	 art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
	 event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
	 if ( !caloCrystalHitsHandle.isValid()) return;

	 //Create a new CaloCluster collection and fill it
	 std::unique_ptr<CaloProtoClusterCollection> caloProtoClustersMain(new CaloProtoClusterCollection);
	 MakeCaloTriggers(*caloProtoClustersMain,caloCrystalHitsHandle);

	 event.put(std::move(caloProtoClustersMain), _producerNameMain);

     }





     void MakeCaloTrigger::MakeCaloTriggers(CaloProtoClusterCollection& caloProtoClustersMain, 
                                            art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle)
     {


            CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
            Calorimeter const & cal = *(GeomHandle<Calorimeter>());

            if (caloCrystalHits.empty()) return;

	    //declare and fill the hash map crystal_id -> list of CaloHits
            std::vector<CaloCrystalList>    mainClusterList;           
	    std::vector<CaloCrystalVec>     caloIdHitMap(cal.nCrystal());
	    CaloSeedManager                 seedList(cal.nCrystal());

	    for (auto const& hit : caloCrystalHits)
	    {
	        if (hit.energyDep() < _EnoiseCut) continue;
	        if (hit.time()      < _timeCut) continue;
	        caloIdHitMap[hit.id()].push_back(&hit);     
	        seedList.add(&hit);     
	    }   


	    //produce main clusters
	    while( CaloCrystalHit const* crystalSeed = seedList.seed() )
	    {
	       if (crystalSeed->energyDep() < _EminSeed) break;

	       CaloClusterFinder finder(cal,crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	       finder.formCluster(caloIdHitMap);	 

	       CaloCrystalList crystalsInCluster = finder.clusterList();	       
	       mainClusterList.push_back(crystalsInCluster);
	       seedList.checkSeedbyList(finder.inspected(),caloIdHitMap);
	    }  

            std::sort(mainClusterList.begin(),mainClusterList.end(),[](CaloCrystalList const& a, CaloCrystalList const& b){return (*a.begin())->time() < (*b.begin())->time();});

	    
	    
	    std::vector<double> clusterEnergy(mainClusterList.size(),0.0);
	    std::vector<int>    nCry(mainClusterList.size(),0);
	    std::vector<int>    associated(mainClusterList.size(),-1);

	    double enerMax(0);
	    for (unsigned int i=0;i<mainClusterList.size(); ++i)
	    {	   	     		
		if (associated[i]>-1) continue;
		
		for (auto const& il : mainClusterList[i]){clusterEnergy[i] += il->energyDep(); nCry[i]++;}	     

		CaloCrystalHit const* hitFirst = *(mainClusterList[i].begin());
		for (unsigned int j=i+1;j<mainClusterList.size();++j)
		{
		    if (associated[j]>-1) continue;
		    
		    CaloCrystalHit const* hitSecond = *(mainClusterList[j].begin());	     
	            if (hitSecond->time()  -  hitFirst->time() > _deltaTimePlus)  break;

		    CLHEP::Hep3Vector crystalPos1 = cal.crystalOrigin(hitFirst->id());
		    CLHEP::Hep3Vector crystalPos2 = cal.crystalOrigin(hitSecond->id());
        	    double dist = (crystalPos1-crystalPos2).mag();

		    if (dist > _maxDistClu) continue;

		    associated[j]=i;		   
		    for (auto const& il : mainClusterList[j]){clusterEnergy[i] += il->energyDep(); nCry[i]++;}
		}	  
	        if (clusterEnergy[i] > enerMax) enerMax = clusterEnergy[i];
	    }
	    
	    
	    _hE->Fill(enerMax);

 
     }




}



using mu2e::MakeCaloTrigger;
DEFINE_ART_MODULE(MakeCaloTrigger);
