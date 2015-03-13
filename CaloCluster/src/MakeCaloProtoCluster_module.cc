
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
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

// Other includes.
#include "cetlib/exception.h"

using std::string;






namespace mu2e {



class MakeCaloProtoCluster : public art::EDProducer {


     public:

             typedef std::vector<CaloCrystalHit const*>  CaloCrystalVec;
             typedef std::list<CaloCrystalHit const*>    CaloCrystalList;


             explicit MakeCaloProtoCluster(fhicl::ParameterSet const& pset) :
             _diagLevel(pset.get<int>("diagLevel",0)),
             _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
             _producerNameMain(pset.get<std::string>("mainClusterCollName")),
             _producerNameSplit(pset.get<std::string>("splitClusterCollName")),
             _EminSeed(pset.get<double>("EminSeed")),
             _EnoiseCut(pset.get<double>("EnoiseCut")),
             _ExpandCut(pset.get<double>("ExpandCut")),
             _deltaTimePlus(pset.get<double>("deltaTimePlus")),
             _deltaTimeMinus(pset.get<double>("deltaTimeMinus")),
             _messageCategory("CLUSTER")
             {
                     produces<CaloProtoClusterCollection>(_producerNameMain);
                     produces<CaloProtoClusterCollection>(_producerNameSplit);
             }


             virtual ~MakeCaloProtoCluster() { }
             virtual void beginJob();
             void produce( art::Event& e);



     private:

             int         _diagLevel;
             std::string _caloCrystalModuleLabel;
             std::string _producerNameMain;
             std::string _producerNameSplit;

             double      _EminSeed;
             double      _EnoiseCut;
             double      _ExpandCut;
             double      _deltaTimePlus;
             double      _deltaTimeMinus;

             const std::string _messageCategory;



             void MakeCaloProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain, 
	                                CaloProtoClusterCollection& caloProtoClustersSplit, 
                                	art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);

             void filterByTime(CaloCrystalVec& vec, std::vector<double>& clusterTime);
             void dump(std::vector<CaloCrystalVec >& caloIdHitMap);



     };



     void MakeCaloProtoCluster::beginJob(){

     }




     void MakeCaloProtoCluster::produce(art::Event& event) {


	 // Check that calorimeter geometry description exists
	 art::ServiceHandle<GeometryService> geom;
	 if( !(geom->hasElement<Calorimeter>()) ) return;

	 //Get handles to calorimeter crystal hits
	 art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
	 event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
	 if ( !caloCrystalHitsHandle.isValid()) return;

	 //Create a new CaloCluster collection and fill it
	 std::unique_ptr<CaloProtoClusterCollection> caloProtoClustersMain(new CaloProtoClusterCollection);
	 std::unique_ptr<CaloProtoClusterCollection> caloProtoClustersSplit(new CaloProtoClusterCollection);
	 MakeCaloProtoClusters(*caloProtoClustersMain,*caloProtoClustersSplit,caloCrystalHitsHandle);

	 event.put(std::move(caloProtoClustersMain),  _producerNameMain);
	 event.put(std::move(caloProtoClustersSplit), _producerNameSplit);

     }





     void MakeCaloProtoCluster::MakeCaloProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain, 
                                                      CaloProtoClusterCollection& caloProtoClustersSplit,
                                	              art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle)
     {



            CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
            Calorimeter const & cal = *(GeomHandle<Calorimeter>());

            if (caloCrystalHits.empty()) return;
            CaloCrystalHit const* caloCrystalHitBase = &caloCrystalHits.front();


	    //declare and fill the hash map crystal_id -> list of CaloHits
            std::vector<CaloCrystalList>    mainClusterList, splitClusterList;
	    std::vector<double>             clusterTime;
	    std::vector<CaloCrystalVec>     caloIdHitMap(cal.nCrystal());
	    CaloSeedManager                 seedList(cal.nCrystal());


	    for (auto const& hit : caloCrystalHits)
	    {
	       if (hit.energyDep() < _EnoiseCut) continue;
	       caloIdHitMap[hit.id()].push_back(&hit);     
	       seedList.add(hit);     
	    }   



	    //produce main clusters
	    while( CaloCrystalHit const* crystalSeed = seedList.seed() )
	    {
	       if (crystalSeed->energyDep() < _EminSeed) break;

	       CaloClusterFinder finder(cal,*crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	       finder.formCluster(caloIdHitMap);	 

	       CaloCrystalList crystalsInCluster = finder.clusterList();
	       mainClusterList.push_back(crystalsInCluster);
               clusterTime.push_back(crystalSeed->time());

	       seedList.checkSeedbyList(finder.inspected(),caloIdHitMap);

	    }  

            
	    //filter unneeded hits
	    for (unsigned i=0; i < caloIdHitMap.size(); ++i )
	    {
	       filterByTime(caloIdHitMap[i], clusterTime);
	       seedList.checkSeedbyId(i,caloIdHitMap[i]);
	     }



	    //produce split-offs clusters
	    while(CaloCrystalHit const* crystalSeed = seedList.seed())
	    {
	       CaloClusterFinder finder(cal,*crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	       
	       finder.formCluster(caloIdHitMap);	 	    
	       splitClusterList.push_back(finder.clusterList());
	       
	       seedList.checkSeedbyList(finder.inspected(),caloIdHitMap);
	    }





	    //save them
	    for (auto main : mainClusterList)
	    {
		if (_diagLevel) std::cout<<"This cluster contains "<<main.size()<<" crystals, id= ";

		bool isSplit(false);
		CaloCrystalHit const* seed  = (*main.begin());
		double seed_time            = seed->time();

		double totalEnergy = 0;
		std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitPtrVector;

		for (auto il = main.begin(); il !=main.end(); ++il)
		{
		  totalEnergy += (*il)->energyDep();
		  size_t idx = (*il - caloCrystalHitBase);
		  caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
		  if (_diagLevel ) std::cout<<(*il)->id()<<" "; 
		}

		CaloProtoCluster caloprotoCluster(seed_time,totalEnergy,caloCrystalHitPtrVector,isSplit);	      
		caloProtoClustersMain.push_back(caloprotoCluster);

		if (_diagLevel) std::cout<<" with energy="<<totalEnergy<<" and time="<<seed_time<<std::endl;;
	    }


	    for (unsigned int i=0;i<splitClusterList.size();++i)
	    {

		if (_diagLevel) std::cout<<"This split-off cluster contains "<<splitClusterList[i].size()<<" crystals, id= ";

		CaloCrystalList& thisList   = splitClusterList[i];

		bool isSplit(true);
		CaloCrystalHit const* seed  = (*thisList.begin());
		double seed_time            = seed->time();

		double totalEnergy = 0;
		std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitPtrVector;

		for (auto il = thisList.begin(); il !=thisList.end(); ++il)
		{
		  totalEnergy += (*il)->energyDep();
		  size_t idx = (*il - caloCrystalHitBase);
		  caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
	          if (_diagLevel ) std::cout<<(*il)->id()<<" "; 
		}

	        CaloProtoCluster caloprotoCluster(seed_time,totalEnergy,caloCrystalHitPtrVector,isSplit);	      
		caloProtoClustersSplit.push_back(caloprotoCluster);

		if (_diagLevel) std::cout<<" with energy="<<totalEnergy<<" and time="<<seed_time<<std::endl;;
	    }

     }





     void MakeCaloProtoCluster::filterByTime(CaloCrystalVec& vec, std::vector<double>& clusterTime)
     {

            for (unsigned int ivec=0; ivec < vec.size(); ++ivec)
	    {	        
		if (vec[ivec]==0) continue;
		CaloCrystalHit const* hit = vec[ivec];
		double timePlus           = _deltaTimeMinus + hit->time();
		double timeMinus          = hit->time() - _deltaTimePlus;

		auto itTime = clusterTime.begin();
		while (itTime != clusterTime.end())
		{
		  if ( timeMinus < *itTime && *itTime < timePlus) break;
		  ++itTime;
		}  

	       if (itTime == clusterTime.end() ) vec[ivec]=0;
	    }
     }

    
     void MakeCaloProtoCluster::dump(std::vector<CaloCrystalVec>& caloIdHitMap)
     {

       for (unsigned int i=0; i<caloIdHitMap.size(); ++i)
	 for (auto j :  caloIdHitMap[i]) std::cout<<i<<" "<<j->id()<<" "<<j->energyDep()<<" "<<j->time()<<std::endl;
     }


}



using mu2e::MakeCaloProtoCluster;
DEFINE_ART_MODULE(MakeCaloProtoCluster);
