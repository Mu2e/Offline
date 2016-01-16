
// C++ includes.
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <unordered_map>
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
#include "CaloCluster/inc/CaloClusterAssociator.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "CaloCluster/inc/CaloClusterMoments.hh"

// Other includes.
#include "cetlib/exception.h"
#include "TH1D.h"



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



	     

  class MakeCaloTrigger : public art::EDProducer {


     public:

             explicit MakeCaloTrigger(fhicl::ParameterSet const& pset) :
             _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel")),
             _EminSeed(pset.get<double>("EminSeed")),
             _EnoiseCut(pset.get<double>("EnoiseCut")),
             _ExpandCut(pset.get<double>("ExpandCut")),
             _timeCut(pset.get<double>("timeCut")),
             _deltaTimePlus(pset.get<double>("deltaTimePlus")),
             _deltaTimeMinus(pset.get<double>("deltaTimeMinus")),
             _maxDistClu(pset.get<double>("maxDistClu")),
             _associatorFancy(pset.get<bool>("associatorFancy")),
             _messageCategory("CaloTrigger"),
             _diagLevel(pset.get<int>("diagLevel",0)),
	     _hE(0)
             {
                     produces<CaloClusterCollection>();
             }


             virtual ~MakeCaloTrigger() { }
             virtual void beginJob();
             void produce( art::Event& e);



     private:
             
             typedef std::unordered_map<unsigned int,std::vector<unsigned int> >  AssoMap;
             typedef std::vector<CaloCrystalHit const*>                           CaloCrystalVec;
             typedef std::list<CaloCrystalHit const*>                             CaloCrystalList;
	     
             
             std::string _caloCrystalModuleLabel;
             double      _EminSeed;
             double      _EnoiseCut;
             double      _ExpandCut;
             double      _timeCut;
             double      _deltaTimePlus;
             double      _deltaTimeMinus;
             double      _maxDistClu;
	     bool        _associatorFancy;

             const std::string _messageCategory;
	     int         _diagLevel;
             TH1F *_hE;


             void MakeCaloTriggers(CaloClusterCollection& caloClusters, 
                                   art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);

             void associateFast(Calorimeter const& cal, AssoMap& associatedMainId, 
                                std::vector<CaloCrystalList> const& clusterColl, double deltaTime, double maxDist);  

             void associateFancy(Calorimeter const& cal, AssoMap& associatedMainId, 
                                 std::vector<CaloCrystalList> const& clusterColl, double deltaTime, double maxDist);  

             double closestDistance(Calorimeter const& cal, CaloCrystalList const& cluster, CaloCrystalList const& cluster2);


     };



     void MakeCaloTrigger::beginJob()
     {
        art::ServiceHandle<art::TFileService> tfs;
        _hE  = tfs->make<TH1F>("cluEner","Cluster energy",150,0,150.);
     }




     void MakeCaloTrigger::produce(art::Event& event)
     {

	 // Check that calorimeter geometry description exists
	 art::ServiceHandle<GeometryService> geom;
	 if( !(geom->hasElement<Calorimeter>()) ) return;

	 //Get handles to calorimeter crystal hits
	 art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
	 event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
	 if ( !caloCrystalHitsHandle.isValid()) return;

	 //Create a new CaloCluster collection and fill it
	 std::unique_ptr<CaloClusterCollection> caloClusters(new CaloClusterCollection);
	 MakeCaloTriggers(*caloClusters,caloCrystalHitsHandle);

	 event.put(std::move(caloClusters));
     }





     void MakeCaloTrigger::MakeCaloTriggers(CaloClusterCollection& caloClusters, 
                                             art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle)
     {


            CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
            Calorimeter const & cal = *(GeomHandle<Calorimeter>());

            if (caloCrystalHits.empty()) return;
            CaloCrystalHit const* caloCrystalHitBase = &caloCrystalHits.front();


	    //declare and fill the hash map crystal_id -> list of CaloHits
            std::vector<CaloCrystalList>    mainClusterList,caloIdHitMap(cal.nCrystal());           
            std::set<CaloCrystalHit const*, caloSeedCompare> seedList;


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
	         for (auto const& hit: finder.clusterList()) seedList.erase(seedList.find(hit));
	    }  

            std::sort(mainClusterList.begin(),mainClusterList.end(),[](CaloCrystalList const& a, CaloCrystalList const& b){return (*a.begin())->time() < (*b.begin())->time();});
	    
	    


	    AssoMap associatedMainId;
	    if (_associatorFancy ) associateFancy(cal,associatedMainId,mainClusterList ,_deltaTimePlus,_maxDistClu);
	    else                   associateFast(cal,associatedMainId,mainClusterList ,_deltaTimePlus,_maxDistClu);




	    std::vector<int> flagProto(mainClusterList.size(),0);
	    for (unsigned int iclu=0; iclu<mainClusterList.size(); ++iclu)
	    { 

		 if (flagProto[iclu]) continue;		  		  

		 bool isSplit(false);
		 double totalEnergy = 0;
		 std::vector<art::Ptr<CaloCrystalHit>> caloCrystalHitPtrVector;


		 CaloCrystalList list = mainClusterList[iclu];
		 for (auto const& il : list)
		 {
		     totalEnergy += il->energyDep();
		     size_t idx = (il - caloCrystalHitBase);
		     caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
		 }


		 for (int iassoc : associatedMainId[iclu])
		 {
		     flagProto[iassoc] = 1;
                     
		     isSplit = true;
		     CaloCrystalList listAssoc = mainClusterList[iassoc];
		     for (auto const& il : listAssoc)
		     {
			 totalEnergy += il->energyDep();
			 size_t idx = (il - caloCrystalHitBase);
			 caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
		     }
		 }

		 std::sort(caloCrystalHitPtrVector.begin(),caloCrystalHitPtrVector.end(),
			   [] (art::Ptr<CaloCrystalHit> const& lhs, art::Ptr<CaloCrystalHit> const& rhs) 
			       {return lhs->energyDep() > rhs->energyDep();} );

		 auto const& seed   = **caloCrystalHitPtrVector.begin();
		 double seed_time   = seed.time();
		 int seed_section   = cal.crystal(seed.id()).sectionId();

		 CaloCluster caloCluster(seed_section, seed_time, totalEnergy, caloCrystalHitPtrVector, isSplit);	      
  	         
		 CaloClusterMoments cogCalculator(cal,caloCluster, seed_section);
		 cogCalculator.calculate(CaloClusterMoments::Linear);
                 caloCluster.cog3Vector(cogCalculator.cog());
		 
		 caloClusters.push_back(caloCluster);


		 _hE->Fill(totalEnergy);
	    }
 
     }







     
     void MakeCaloTrigger::associateFast(Calorimeter const& cal, AssoMap& associatedMainId, std::vector<CaloCrystalList> const& clusterColl, 
                                          double deltaTime, double maxDist)  
     { 

	std::vector<int>  isAssociatedTo(clusterColl.size(),-1);

        for (unsigned int i=0;i<clusterColl.size(); ++i)
	{	   	     	       
	    associatedMainId[i].clear();
            CaloCrystalHit const* hitFirst = *(clusterColl[i].begin());

	    for (unsigned int j=i+1;j<clusterColl.size();++j)
	    {
 	         CaloCrystalHit const* hitSecond = *(clusterColl[j].begin());	     
		 if (isAssociatedTo[j] > -1) continue;
	      
	         if (hitSecond->time()  -  hitFirst->time() > deltaTime)  break;
		 CLHEP::Hep3Vector crystalPos1 = cal.crystal(hitFirst->id()).position();
		 CLHEP::Hep3Vector crystalPos2 = cal.crystal(hitSecond->id()).position();
        	 double dist = (crystalPos1-crystalPos2).mag();

		 if (dist > maxDist) continue;

		 isAssociatedTo[j] = i;
		 associatedMainId[i].push_back(j); 
	    }	  
	}
	
	return;
    } 




     void MakeCaloTrigger::associateFancy(Calorimeter const& cal, AssoMap& associatedMainId, std::vector<CaloCrystalList> const& clusterColl, 
                                           double deltaTime, double maxDist)  
     { 

	std::vector<int>                isAssociatedTo(clusterColl.size(),-1);
	std::vector<std::vector<int> >  associatedId(clusterColl.size());


        for (unsigned int i=0;i<clusterColl.size(); ++i)
	{	   	     	       
	    associatedMainId[i].clear();
            CaloCrystalHit const* hitFirst = *(clusterColl[i].begin());

	    for (unsigned int j=i+1;j<clusterColl.size();++j)
	    {
 	         CaloCrystalHit const* hitSecond = *(clusterColl[j].begin());	     
		 if (isAssociatedTo[j] > -1) continue;
	      
	         if (hitSecond->time()  -  hitFirst->time() > deltaTime)  break;

		 double dist = closestDistance(cal, clusterColl[i], clusterColl[j]);
		 if (dist > maxDist) continue;

		 isAssociatedTo[j] = i;
		 associatedId[i].push_back(j);
		 //associatedMainId[i].push_back(j); // and remove next snippet of code
	    }	  
	}


	for (unsigned int i=0;i<associatedId.size();++i)
	{	       
	     std::set<int> neighbors;
	     std::queue<int> list;
	     for (int id : associatedId[i]){neighbors.insert(id); list.push(id);}

	     while (!list.empty())
	     {
	          int nextId = list.front();
		  for (int id : associatedId[nextId]) {neighbors.insert(id); list.push(id);}
	          associatedId[nextId].clear();  
		  list.pop();	       
	     }
	     for (int id : neighbors) associatedMainId[i].push_back(id);
	}
	
	return;
    } 



    double MakeCaloTrigger::closestDistance(Calorimeter const& cal, CaloCrystalList const& cluster, CaloCrystalList const& cluster2)
    {
       double minDistance(1e6);    
       for (auto const& hit : cluster)
       {        
	   CLHEP::Hep3Vector crystalPos = cal.crystal(hit->id()).position();

	   for (auto const& hit2 : cluster2)
	   {
	      double dist = (crystalPos - cal.crystal(hit2->id()).position()).mag();
	      if (dist<minDistance) minDistance = dist;
	   }	 
       }

      return minDistance;
    }



}



using mu2e::MakeCaloTrigger;
DEFINE_ART_MODULE(MakeCaloTrigger);
