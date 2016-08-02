
// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CaloCluster/inc/ClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// Other includes.
#include "cetlib/exception.h"
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <memory>




namespace mu2e {


 class CaloClusterFast : public art::EDProducer {


      public:

              typedef std::vector<CaloCrystalHit const*>  CaloCrystalVec;
              typedef std::list<CaloCrystalHit const*>    CaloCrystalList;


              explicit CaloClusterFast(fhicl::ParameterSet const& pset) :
              _diagLevel(pset.get<int>("diagLevel",0)),
              _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
              _EminSeed(pset.get<double>("EminSeed")),
              _EnoiseCut(pset.get<double>("EnoiseCut")),
              _ExpandCut(pset.get<double>("ExpandCut")),
              _timeCut(pset.get<double>("timeCut")),
              _deltaTimePlus(pset.get<double>("deltaTimePlus")),
              _deltaTimeMinus(pset.get<double>("deltaTimeMinus")),
              _messageCategory("CLUSTER")
              {
                      produces<CaloClusterCollection>();
              }


              virtual ~CaloClusterFast() { }
              void produce( art::Event& e);



      private:

              int               _diagLevel;
              std::string       _caloCrystalModuleLabel;
              double            _EminSeed;
              double            _EnoiseCut;
              double            _ExpandCut;
              double            _timeCut;
              double            _deltaTimePlus;
              double            _deltaTimeMinus;
              const std::string _messageCategory;


              void makeClusters(CaloClusterCollection& caloClusters, 	                                 
                                art::Handle<CaloCrystalHitCollection> const& CaloCrystalHitsHandle);

      };





      void CaloClusterFast::produce(art::Event& event) {

	  // Check that calorimeter geometry description exists
	  art::ServiceHandle<GeometryService> geom;
	  if( !(geom->hasElement<Calorimeter>()) ) return;

	  //Get handles to calorimeter crystal hits
	  art::Handle<CaloCrystalHitCollection> CaloCrystalHitsHandle;
	  event.getByLabel(_caloCrystalModuleLabel, CaloCrystalHitsHandle);
	  if ( !CaloCrystalHitsHandle.isValid()) return;

	  //Create a new CaloCluster collection and fill it
	  std::unique_ptr<CaloClusterCollection> caloClusters(new CaloClusterCollection);
	  makeClusters(*caloClusters,CaloCrystalHitsHandle);

	  event.put(std::move(caloClusters));

      }





      //----------------------------------------------------------------------------------------------------------
      void CaloClusterFast::makeClusters(CaloClusterCollection& caloClusters,                                                        
                                	  art::Handle<CaloCrystalHitCollection> const& CaloCrystalHitsHandle)
      {

             CaloCrystalHitCollection const& CaloCrystalHits(*CaloCrystalHitsHandle);
             Calorimeter const & cal = *(GeomHandle<Calorimeter>());

             if (CaloCrystalHits.empty()) return;


	     //declare and fill the hash map crystal_id -> list of CaloHits
             std::vector<CaloCrystalList>       mainClusterList, splitClusterList,caloIdHitMap(cal.nCrystal());
             std::list<CaloCrystalHit const*>   seedList;
	     std::vector<double>                clusterTime;


	     //fill data structures
	     for (auto const& hit : CaloCrystalHits)
	     {
	         if (hit.energyDep() < _EnoiseCut || hit.time() < _timeCut) continue;
	         caloIdHitMap[hit.id()].push_back(&hit);      
	         seedList.push_back(&hit);
	     }   

	     seedList.sort([](CaloCrystalHit const* a, CaloCrystalHit const* b) {return a->energyDep() > b->energyDep();});


	     //produce main clusters
	     while( !seedList.empty() )
	     {
	         CaloCrystalHit const* crystalSeed = *seedList.begin();
	         if (crystalSeed->energyDep() < _EminSeed) break;

	         ClusterFinder finder(cal,crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	         finder.formCluster(caloIdHitMap);	 

	         mainClusterList.push_back(finder.clusterList());
                 clusterTime.push_back(crystalSeed->time());

	         for (auto const& hit: finder.clusterList()) seedList.remove(hit);
	     }  


	     //save these clusters
	     for (auto main : mainClusterList)
	     {
		 CaloCrystalHit const* seed  = (*main.begin());
		 double seed_time            = seed->time();
                 int seed_section   = cal.crystal(seed->id()).sectionId();
		 std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitsPtrVector;

		 double totalEnergy = 0;
		 for (auto il = main.begin(); il !=main.end(); ++il)totalEnergy += (*il)->energyDep();
                
		 caloClusters.push_back(CaloCluster(seed_section,seed_time,totalEnergy,caloCrystalHitsPtrVector,0));

	     }

	     std::sort(caloClusters.begin(),  caloClusters.end(), [](CaloCluster const& a, CaloCluster const& b) {return a.time() < b.time();});
      }









}



using mu2e::CaloClusterFast;
DEFINE_ART_MODULE(CaloClusterFast);



