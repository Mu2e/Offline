// To do list
// Must form a calo Cluster object
/*

- get the main cluster
  - loop over crystals, find the ones with a time compatible with the cluster time
  - use them to find the other clusters
repeat

- find all the main clusters
- filter by time
- form the clusters and connect them to main clusters

  - 
  
   


*/

// then must compute COG / direction / other crap
// then must check effciency
// then must redo everything to include splitting
// and finally get a beer for all these trouble
/*
 * MakeCaloClusterNew3_module.cc
 *
 *  Created on: Feb 10, 2012
 *      Author: echenard
 */

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <queue>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>


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
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

//calorimeter packages
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "CaloCluster/inc/CaloClusterCogCalculator.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "CaloCluster/inc/CaloClusterFinderNew.hh"
#include "CaloCluster/inc/CaloSeedManager.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"
#include "TMath.h"







#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"


#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"


using namespace std;

namespace mu2e {



class MakeCaloClusterNew3 : public art::EDProducer {


   public:
           typedef std::vector<CaloCrystalHit const*>                CaloCrystalVec;
           typedef std::list<CaloCrystalHit const*>                  CaloCrystalList;
           
           
           explicit MakeCaloClusterNew3(fhicl::ParameterSet const& pset) :

           // Parameters
           _diagLevel(pset.get<int>("diagLevel",0)),
           _maxFullPrint(pset.get<int>("maxFullPrint",5)),
           _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)),
           _deltaTimePlus(pset.get<double>("deltaTimePlus", 10.)),// ns
           _deltaTimeMinus(pset.get<double>("deltaTimeMinus", 10.)),// ns
           _nCryPerCluster(pset.get<int>("nCryPerCrystal", 0)),
           _EnoiseCut(pset.get<double>("EnoiseCut", 0.090)),//MeV 3 sigma noise
           _ExpandCut(pset.get<double>("ExpandCut", 0.090)),//MeV
           _EminSeed(pset.get<double>("EminSeed", 10)),//MeV
           _maxDist(pset.get<double>("maxDist", 10000)),//MeV
	   _doAssociate(pset.get<bool>("doAssociate", true)),
           _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
           _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
           _messageCategory("HITS"),
           _firstEvent(true)
           {
                   // Tell the framework what we make.
                   produces<CaloClusterCollection>();
           }
           
           
           virtual ~MakeCaloClusterNew3() { }
           virtual void beginJob();
           void produce( art::Event& e);



   private:
   
 
           // Diagnostics level.
           int _diagLevel;

           // Limit on number of events for which there will be full printout.
           int _maxFullPrint;

           // Name of the calorimeter StepPoint collection
           std::string _caloStepPoints;

           // Parameters
           double _minimumEnergy;  // minimum energy deposition of G4 step
           double _deltaTimePlus;
           double _deltaTimeMinus;
           double _nCryPerCluster;
           double _EnoiseCut;
           double _ExpandCut;
           double _EminSeed;
           double _maxDist;
	   bool   _doAssociate;
           string _g4ModuleLabel;  // Name of the module that made these hits.
           string _caloReadoutModuleLabel;
           string _caloCrystalModuleLabel;

           // A category for the error logger.
           const std::string _messageCategory;

           // Give some informationation messages only on the first event.
           bool _firstEvent;

           void makeCaloClusters(CaloClusterCollection& caloClusters, 
                                 art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);


           CLHEP::Hep3Vector calculateCog(std::list<CaloCrystalHit const*> cluster, Calorimeter const & cal, int mode=1); 
           double closestDistance(std::list<CaloCrystalHit const*> cluster, std::list<CaloCrystalHit const*> cluster2, Calorimeter const & cal) ;
           void cleanHits(std::vector<CaloCrystalVec >& CaloIdHitMap, std::vector<double>& clusterTime, std::set<int>& neighbors);
           void filterByTime(CaloCrystalVec& vec, std::vector<double>& clusterTime);

           void dump(std::vector<CaloCrystalVec >& caloIdHitMap);
           
	   double closestDistance(Calorimeter const & cal,CaloCrystalList& cluster, CaloCrystalList& cluster2);
           double distance(Calorimeter const & cal, CaloCrystalHit const* hit1, CaloCrystalHit const* hit2);
	   


};


   void MakeCaloClusterNew3::beginJob(){
           

   }




   void MakeCaloClusterNew3::produce(art::Event& event) {
          

       // Check that calorimeter geometry description exists
       art::ServiceHandle<GeometryService> geom;
       if( !(geom->hasElement<Calorimeter>()) ) return;

       //Get handles to calorimeter crystal hits
       art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
       event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
       if ( !caloCrystalHitsHandle.isValid()) return;

      //Create a new CaloCluster collection and fill it
       unique_ptr<CaloClusterCollection> caloClusters(new CaloClusterCollection);
       makeCaloClusters(*caloClusters,caloCrystalHitsHandle);

       event.put(std::move(caloClusters));

   }
   
   
                     
   
   
   void MakeCaloClusterNew3::makeCaloClusters(CaloClusterCollection& caloClusters,
                                      art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle) {


     
          CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
          Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  
          if (caloCrystalHits.empty()) return;
          CaloCrystalHit const* caloCrystalHitBase = &caloCrystalHits.front();



	  //declare and fill the hash map crystal_id -> list of CaloHits
          std::vector<CaloCrystalList>                mainClusterList, smallClusterList;
	  std::vector<double>                         clusterTime;
	  std::vector<CaloCrystalVec>                 caloIdHitMap(cal.nCrystal());
	  CaloSeedManager                             seedList(cal.nCrystal());
          std::set<int>                               clusterNeighborsList;	  


	  for (auto const& hit : caloCrystalHits)
	  {
	     //std::cout<<hit.id()<<" "<<hit.time()<<" "<<hit.energyDep()<<std::endl;
	     if (hit.energyDep() < _EnoiseCut) continue;
	     caloIdHitMap[hit.id()].push_back(&hit);     
	     seedList.add(hit);     
	  }   
	  	  
		  
		  	  	  
	  while(CaloCrystalHit const* crystalSeed = seedList.seed())
	  {
             if (crystalSeed->energyDep() < _EminSeed) break;

	     CaloClusterFinderNew finder(cal,*crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	     finder.formCluster(caloIdHitMap);	 
	     CaloCrystalList crystalsInCluster = finder.clusterList();
	     mainClusterList.push_back(crystalsInCluster);
             if (!crystalsInCluster.empty())  clusterTime.push_back((*crystalsInCluster.begin())->time());
	     seedList.checkSeedbyList(finder.inspected(),caloIdHitMap);

	     for (auto const& i :  crystalsInCluster)
	     {
	        for (auto const& j : cal.nextNeighbors(i->id())) clusterNeighborsList.insert(j);
	        for (auto const& j : cal.nextNextNeighbors(i->id())) clusterNeighborsList.insert(j);
	     }	
			    
	  }  


	  cleanHits(caloIdHitMap,clusterTime,clusterNeighborsList);
	  for (unsigned i=0; i < caloIdHitMap.size(); ++i ) seedList.checkSeedbyId(i,caloIdHitMap[i]);



	  while(CaloCrystalHit const* crystalSeed = seedList.seed())
	  {
	     CaloClusterFinderNew finder(cal,*crystalSeed,_deltaTimePlus,_deltaTimeMinus, _ExpandCut);
	     finder.formCluster(caloIdHitMap);	 	    
	     smallClusterList.push_back(finder.clusterList());
	     seedList.checkSeedbyList(finder.inspected(),caloIdHitMap);
	  }





	  //run associator
	  std::map<int,int> associatedId;
	  std::map<int,double> associatedDist;
          for (unsigned int i=0;i<smallClusterList.size(); ++i)
	  {	   
	     
	     CaloCrystalHit const* hitSmall = *(smallClusterList[i].begin());
	     
	     double minDist(1e6);
	     int jmin(-1);
	     for (unsigned int j=0;j<mainClusterList.size(); ++j)
	     {
  	        CaloCrystalHit const* hitMain = *(mainClusterList[j].begin());
                
	        if (hitSmall->time() -  hitMain->time()  > _deltaTimePlus)  continue;
        	if (hitMain->time()  -  hitSmall->time() > _deltaTimeMinus) continue;
		
		double dist =  closestDistance(cal,mainClusterList[j],smallClusterList[i]);
		
		if (dist < _maxDist && dist < minDist) {minDist=dist; jmin=j;} 
	     }	  
             associatedId[i]   = jmin;
             associatedDist[i] = minDist;
	  }








	  for (auto main : mainClusterList)
	  {

	    CaloCrystalHit const* seed  = (*main.begin());
	    double seed_time            = seed->time();
	    int seed_section            = cal.caloSectionId(seed->id());

	    double totalEnergy = 0;
	    std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitPtrVector;

	    for (auto il = main.begin(); il !=main.end(); ++il)
	    {
	      totalEnergy += (*il)->energyDep();
	      size_t idx = (*il - caloCrystalHitBase);
	      caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
	    }

	    CaloCluster caloCluster(seed_section,seed_time,totalEnergy,caloCrystalHitPtrVector);	      
	    caloCluster.SetDistance(0);	      
	    caloCluster.SetParentId(-1);	      
	    caloClusters.push_back(caloCluster);
	  }


	  for (unsigned int i=0;i<smallClusterList.size();++i)
	  {

	      CaloCrystalList& thisList   = smallClusterList[i];

	      CaloCrystalHit const* seed  = (*thisList.begin());
	      double seed_time            = seed->time();
	      int seed_section            = cal.caloSectionId(seed->id());

	      double totalEnergy = 0;
	      std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitPtrVector;

	      for (auto il = thisList.begin(); il !=thisList.end(); ++il)
	      {
		totalEnergy += (*il)->energyDep();
		size_t idx = (*il - caloCrystalHitBase);
		caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
	      }

	      CaloCluster caloCluster(seed_section,seed_time,totalEnergy,caloCrystalHitPtrVector);	      
	      caloCluster.SetDistance(associatedDist[i]);	      
	      caloCluster.SetParentId(associatedId[i]);	      
	      caloClusters.push_back(caloCluster);
	  }
  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   }

   
   

   
   
   
   void MakeCaloClusterNew3::cleanHits(std::vector<CaloCrystalVec >& caloIdHitMap, std::vector<double>& clusterTime, std::set<int>& neighbors)
   {
       
       for (unsigned i=0; i < caloIdHitMap.size(); ++i )
       {
           
	    if (_doAssociate && _maxDist<100 && neighbors.find(i) == neighbors.end()) {caloIdHitMap[i].clear(); continue;}

	    filterByTime(caloIdHitMap[i], clusterTime);
	                
       }          
   }   



   void MakeCaloClusterNew3::filterByTime(CaloCrystalVec& vec, std::vector<double>& clusterTime)
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


   double MakeCaloClusterNew3::closestDistance(Calorimeter const & cal, CaloCrystalList& cluster, CaloCrystalList& cluster2)
   {
      double minDistance(1e6);    
      for (auto const& hit : cluster)
      {        
	  CLHEP::Hep3Vector crystalPos = cal.crystalOrigin(hit->id());

	  for (auto const& hit2 : cluster2)
	  {
	     CLHEP::Hep3Vector crystalPos2 = cal.crystalOrigin(hit2->id());
	     double dist = (crystalPos-crystalPos2).mag();
	     if (dist<minDistance) minDistance = dist;
	  }	 
      }

     return minDistance;
   }


   double MakeCaloClusterNew3::distance(Calorimeter const & cal, CaloCrystalHit const* hit1, CaloCrystalHit const* hit2)
   {
 	  CLHEP::Hep3Vector crystalPos1 = cal.crystalOrigin(hit1->id());
 	  CLHEP::Hep3Vector crystalPos2 = cal.crystalOrigin(hit2->id());
	  return (crystalPos1-crystalPos2).mag();
   }






   void MakeCaloClusterNew3::dump(std::vector<CaloCrystalVec>& caloIdHitMap)
   {

     for (unsigned int i=0; i<caloIdHitMap.size(); ++i)
       for (auto j :  caloIdHitMap[i]) std::cout<<i<<" "<<j->id()<<" "<<j->energyDep()<<" "<<j->time()<<std::endl;
   }


}// end namespace mu2e



using mu2e::MakeCaloClusterNew3;
DEFINE_ART_MODULE(MakeCaloClusterNew3);
