//
// The clustering is performed in several stages, using the folowing data structures:
//  - bidimensional array that have a list of energy deposits for each crystal id (prefilter low energy deposits to speed up things)
//  - a map containing the potential seeds (the map is ordered by most energetic to lowest energetic hits)

// The clustering proceeds in three steps: form energetic proto-clusters, look at split-offs and form final clusters. The first two steps are done by MakeCaloProtoCluster
// and produce proto-clusters. The last step is done by MakeCaloCluster, and procudes a cluster

// 1. Energetic proto-clusters (clusters above some threshold energy)
//    - start from the most energetic seed (over some threshold)
//    - form a proto-cluster by adding all simply connected cluster to the seed (simply connected = any two hits in a cluster can be joined 
//      by a continuous path of clusters in the crystal)
//    - mark the correpsonding hits as used, and update the seed list
//
// 2. Split-offs: some clusters might have low-energy split-offs, and we need to find them
//    - filter the remaining unassigned hits to retain only those compatible with the time of the main clusters (there are a lot of background low energy deposits
//      incompatible with any clsuter, we don't want them)
//    - update the seed lists and form all remaining clusters as before
//
// 3. Cluster formation
//    - merge proto-clusters and split-offs to form clusters ifg they are "close" enough
//
// Note: 1. One can try to filter the hits first, before producing energetic proto-clusters (use two iterators on the time ordered crystal list, see commented code at the end). 
//          The performance difference is very small compared to this implementation, but more obscure, so for clarity, I kept this one.
//       2. I tried a bunch of other optimizations but the performance increase is small and the code harder to read, so I left it as is
//
// Original author: B. Echenard
//


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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

//calorimeter packages
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/CaloClusterFinder.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

// Other includes.
#include "cetlib/exception.h"


namespace
{
   struct caloSeedCompare {
     bool operator() (mu2e::CaloCrystalHit const* a, mu2e::CaloCrystalHit const* b) const 
     {
        if (std::abs(a->energyDep() - b->energyDep()) > 1e-6) return a->energyDep() > b->energyDep();
        if (a->id() != b->id()) return a->id() > b->id();
        return a->time() > b->time();
     }
   };
}




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
              _timeCut(pset.get<double>("timeCut")),
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

              int               _diagLevel;
              std::string       _caloCrystalModuleLabel;
              std::string       _producerNameMain;
              std::string       _producerNameSplit;
              double            _EminSeed;
              double            _EnoiseCut;
              double            _ExpandCut;
              double            _timeCut;
              double            _deltaTimePlus;
              double            _deltaTimeMinus;
              const std::string _messageCategory;



              void MakeCaloProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain, 
	                                 CaloProtoClusterCollection& caloProtoClustersSplit, 
                                	 art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);

              void filterByTime(CaloCrystalList& liste, 
	                        std::vector<double> const& clusterTime, 
	                        std::set<CaloCrystalHit const*, caloSeedCompare>& seedListe);

              void dump(std::vector<CaloCrystalList> const& caloIdHitMap);



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





      //----------------------------------------------------------------------------------------------------------
      void MakeCaloProtoCluster::MakeCaloProtoClusters(CaloProtoClusterCollection& caloProtoClustersMain, 
                                                       CaloProtoClusterCollection& caloProtoClustersSplit,
                                	               art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle)
      {



             CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
             Calorimeter const & cal = *(GeomHandle<Calorimeter>());

             if (caloCrystalHits.empty()) return;
             CaloCrystalHit const* caloCrystalHitBase = &caloCrystalHits.front();


	     //declare and fill the hash map crystal_id -> list of CaloHits
             std::vector<CaloCrystalList>                     mainClusterList, splitClusterList,caloIdHitMap(cal.nCrystal());
             std::set<CaloCrystalHit const*, caloSeedCompare> seedList;
	     std::vector<double>                              clusterTime;



	     //fill data structures
	     for (auto const& hit : caloCrystalHits)
	     {
	         if (hit.energyDep() < _EnoiseCut || hit.time() < _timeCut) continue;
	         caloIdHitMap[hit.id()].push_back(&hit);      
	         seedList.insert(&hit);
	     }   


	     //produce main clusters
	     while( !seedList.empty() )
	     {
	         CaloCrystalHit const* crystalSeed = *seedList.begin();
	         if (crystalSeed->energyDep() < _EminSeed) break;

	         CaloClusterFinder finder(cal,crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	         finder.formCluster(caloIdHitMap);	 

	         mainClusterList.push_back(finder.clusterList());
                 clusterTime.push_back(crystalSeed->time());

	         for (auto const& hit: finder.clusterList()) seedList.erase(seedList.find(hit));
	     }  


	     //filter unneeded hits
	     for (unsigned i=0; i < caloIdHitMap.size(); ++i ) filterByTime(caloIdHitMap[i], clusterTime, seedList);


	     //produce split-offs clusters
	     while(!seedList.empty())
	     {
                 CaloCrystalHit const* crystalSeed = *seedList.begin();
	         CaloClusterFinder finder(cal,crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);

	         finder.formCluster(caloIdHitMap);	 	    
	         splitClusterList.push_back(finder.clusterList());

	         for (auto const& hit: finder.clusterList()) seedList.erase(seedList.find(hit));
	     }





	     //save these clusters
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

		 caloProtoClustersMain.push_back(CaloProtoCluster(seed_time,totalEnergy,caloCrystalHitPtrVector,isSplit));

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

		 caloProtoClustersSplit.push_back(CaloProtoCluster(seed_time,totalEnergy,caloCrystalHitPtrVector,isSplit));

		 if (_diagLevel) std::cout<<" with energy="<<totalEnergy<<" and time="<<seed_time<<std::endl;;
	     }

	     std::sort(caloProtoClustersMain.begin(),  caloProtoClustersMain.end(), [](CaloProtoCluster const& a, CaloProtoCluster const& b) {return a.time() < b.time();});
	     std::sort(caloProtoClustersSplit.begin(), caloProtoClustersSplit.end(),[](CaloProtoCluster const& a, CaloProtoCluster const& b) {return a.time() < b.time();});

      }





      //----------------------------------------------------------------------------------------------------------
      void MakeCaloProtoCluster::filterByTime(CaloCrystalList& liste, std::vector<double> const& clusterTime, std::set<CaloCrystalHit const*, caloSeedCompare>& seedList)
      {

             for (auto it = liste.begin(); it != liste.end(); ++it)
	     {	        
		  CaloCrystalHit const* hit = *it;
		  double timePlus           =  _deltaTimeMinus + hit->time();
		  double timeMinus          =  hit->time() - _deltaTimePlus;

		  auto itTime = clusterTime.begin();
		  while (itTime != clusterTime.end())
		  {
		      if ( timeMinus < *itTime && *itTime < timePlus) break;
		      ++itTime;
		  }  

	          if (itTime == clusterTime.end() ) {seedList.erase(seedList.find(hit)); liste.erase(it); --it;} 		 
	     }


      }


      //----------------------------------------------------------------------------------------------------------
      void MakeCaloProtoCluster::dump(std::vector<CaloCrystalList> const& caloIdHitMap)
      {

	for (unsigned int i=0; i<caloIdHitMap.size(); ++i)
	{
	    if (caloIdHitMap[i].size()>0) std::cout<<"ProtoCluster crystal "<<i<<" has size "<<caloIdHitMap[i].size()<<std::endl;
	    for (auto j :  caloIdHitMap[i]) std::cout<<i<<" "<<j->id()<<" "<<j->energyDep()<<" "<<j->time()<<std::endl;
	}
      }


}



using mu2e::MakeCaloProtoCluster;
DEFINE_ART_MODULE(MakeCaloProtoCluster);







/*

// Just on case 
// This is a snippet of code to include only the hits compatible with the time of potential seeds in the map. 
// The gain in performance is small, and I find it obscures the code, so I left the old version

	    
            //fast forward iterator until first crystal in time
	    std::vector<CaloCrystalHit>::const_iterator allCrystal = caloCrystalHits.begin();
            while (allCrystal->time() < _timeCut && allCrystal != caloCrystalHits.end()) ++allCrystal;
	    if (allCrystal == caloCrystalHits.end()) return;

            
            //fast forward iterator until first high energy crystal
	    std::vector<CaloCrystalHit>::const_iterator highCrystal = allCrystal;
	    while (highCrystal->energyDep() < _EminSeed && highCrystal != caloCrystalHits.end()) ++highCrystal;
	    if (highCrystal == caloCrystalHits.end() ) return;


	    double maxTime = highCrystal->time() + _deltaTimePlus;
	    double minTime = highCrystal->time() - _deltaTimeMinus;
	    
	    while( allCrystal != caloCrystalHits.end() )
	    {
		if (allCrystal->time() > maxTime && highCrystal != caloCrystalHits.end())
		{
		     do ++highCrystal; while ( highCrystal != caloCrystalHits.end() && highCrystal->energyDep() < _EminSeed);
                     if (highCrystal == caloCrystalHits.end()) --highCrystal;
		     maxTime = highCrystal->time() + _deltaTimePlus;
		     minTime = highCrystal->time() - _deltaTimeMinus;		     
		     continue;
		}
	        
		if (allCrystal->energyDep() > _EnoiseCut && allCrystal->time() > minTime )
		{
		    CaloCrystalHit const* hit = &(*allCrystal);
		    caloIdHitMap[hit->id()].push_back(hit);     
		    seedList.add(hit);			 
		} 
		++allCrystal;
	    }
*/
