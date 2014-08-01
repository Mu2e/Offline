// To do list
// Must form a calo Cluster object

// then must compute COG / direction / other crap
// then must check effciency
// then must redo everything to include splitting
// and finally get a beer for all these trouble
//
// $Id: MakeCaloCluster_module.cc,v 1.10 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//

/*
 * MakeCaloCluster_module.cc
 *
 *  Created on: Feb 10, 2012
 *      Author: gianipez
 */

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <list>
#include <queue>
#include <vector>
#include <algorithm>

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
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "CaloCluster/inc/CaloClusterFinder.hh"
//#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

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


bool caloCrystalHitEnergyPredicate( CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();}


class MakeCaloCluster : public art::EDProducer {


   public:
           
	   typedef std::list<CaloCrystalHit const*>                  CaloCrystalList;
           typedef std::list<CaloCrystalHit const*>::iterator        CaloCrystalList_iter;
           typedef std::list<CaloCrystalHit const*>::const_iterator  CaloCrystalList_const_iter;
           
	   
	   explicit MakeCaloCluster(fhicl::ParameterSet const& pset) :

           // Parameters
           _diagLevel(pset.get<int>("diagLevel",0)),
           _maxFullPrint(pset.get<int>("maxFullPrint",5)),
           _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)),
           _deltaTime(pset.get<double>("deltaTime", 100.)),// ns
           _nCryPerCluster(pset.get<int>("nCryPerCrystal", 0)),
           _EnoiseCut(pset.get<double>("EnoiseCut", 0.500)),//MeV 3 sigma noise
           _ExpandCut(pset.get<double>("ExpandCut", 1.000)),//MeV
           _EminCluster(pset.get<double>("EminCluster", 10)),//MeV
           _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
           _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
           _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
           _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
           _producerName("Algo"+TOUpper(_caloClusterAlgorithm)+"SeededBy"+TOUpper(_caloClusterSeeding)),
           _messageCategory("HITS"),
           _firstEvent(true)
	   {
                   // Tell the framework what we make.
                   produces<CaloClusterCollection>(_producerName);
           }
           
	   
	   virtual ~MakeCaloCluster() { }
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
           double _deltaTime;
           double _nCryPerCluster;
           double _EnoiseCut;
           double _ExpandCut;
           double _EminCluster;
           string _g4ModuleLabel;  // Name of the module that made these hits.
           string _caloReadoutModuleLabel;
           string _caloCrystalModuleLabel;
           string _caloClusterAlgorithm;
           string _caloClusterSeeding;
           const string _producerName;

           // A category for the error logger.
           const std::string _messageCategory;

           // Give some informationation messages only on the first event.
           bool _firstEvent;

 	   void makeCaloClusters(CaloClusterCollection& caloClusters,
                	         art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);

	   std::list<CaloCrystalHit const*> formCluster(CaloCrystalList& hitList,Calorimeter const & cal);
	   std::list<CaloCrystalHit const*>::iterator findCrystalSeed(CaloCrystalList& hitList);
           CLHEP::Hep3Vector calculateCog(std::list<CaloCrystalHit const*> cluster, Calorimeter const & cal, int mode=1); 


   };


//must include a check of the parameters here
   void MakeCaloCluster::beginJob(){
           
	if ( (_caloClusterSeeding.compare("TIME") != 0 ) &&  (_caloClusterSeeding.compare("ENERGY") != 0 ) )	     
	   throw cet::exception("CaloCluster") << "caloClusterSeeding must be either \"energy\" or \"time\" \n";

	   cout<<"selected clustering algorithm--> "<<_caloClusterAlgorithm <<", seeded by "<< _caloClusterSeeding<<endl;
	   //        cout << "Diaglevel: "
           //                        << _diagLevel << " "
           //                        << _maxFullPrint
           //                        << endl;

   }




   void MakeCaloCluster::produce(art::Event& event) {
          

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

       event.put(std::move(caloClusters), _producerName);

   }
   
   
		     
   
   
   void MakeCaloCluster::makeCaloClusters(CaloClusterCollection& caloClusters,
                                      art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle) {


     
          CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
          Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  
  
          //Get a working copy of the CaloCrystalHits
	  std::list<CaloCrystalHit const*> caloCrystalHitsWork;
	  for ( CaloCrystalHitCollection::const_iterator i=caloCrystalHits.begin(); i!=caloCrystalHits.end(); ++i )
	    caloCrystalHitsWork.push_back( &(*i));


          //store the proto clusters here
	  std::vector< CaloCrystalList > protoClusterList;    

	  while( ! caloCrystalHitsWork.empty() ){       
	    CaloCrystalList cluster = formCluster(caloCrystalHitsWork,cal);
	    cluster.sort(mu2e::caloCrystalHitEnergyPredicate);
	    if ((*cluster.begin())->energyDep()<_EminCluster) break;

	    //std::cout<<"Cluster "<<protoClusterList.size()<<std::endl;
	    //for(CaloCrystalList_const_iter i = cluster.begin(); i!= cluster.end();++i)
	    // std::cout<<(*i)->id()<<" "<<(*i)->energyDep()<<" "<<(*i)->time()<<" "<<cal.getCrystalOrigin((*i)->id())<<std::endl;

	    protoClusterList.push_back(cluster); 
	  }


          CaloCrystalHit const* caloCrystalHitBase = &caloCrystalHits.front();


          //Very approximative method to produce the cluster, enough for now
          for (std::vector< CaloCrystalList >::const_iterator it = protoClusterList.begin(); it!=protoClusterList.end();++it){
	     
	      int isection(-1);
	      double totalEnergy(0),averageTime(0);
	      CaloCrystalList thisList = (*it);
              
	      std::vector< art::Ptr< CaloCrystalHit> > caloCrystalHitPtrVector;

	      for (CaloCrystalList_const_iter il = thisList.begin(); il !=thisList.end(); ++il){

		CaloCrystalHit const* hit = *il;
		totalEnergy += hit->energyDep();
		averageTime += hit->time();
		if (isection==-1) isection = cal.caloSectionId(hit->id());
                
		size_t idx = ( hit - caloCrystalHitBase );
		caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );

	      }
	      averageTime /= float(thisList.size());

              if (thisList.size()<2) continue;
	      
	      CLHEP::Hep3Vector cog = calculateCog(thisList,cal,2);

	      CaloCluster caloCluster(isection,averageTime,totalEnergy,caloCrystalHitPtrVector);	      
      	      caloCluster.SetCog3Vector(cog);
	      caloClusters.push_back(caloCluster);

	      //std::cout<<" Cluster cog "<<cog)<<std::endl;
	      
	  }


   }






   std::list<CaloCrystalHit const*> MakeCaloCluster::formCluster(CaloCrystalList& hitList, Calorimeter const & cal)
   {
       
       CaloCrystalList   crystalsInCluster;
       std::queue<CaloCrystalHit const*>  crystalToVisit;

       CaloCrystalList_iter iseed = findCrystalSeed(hitList);
       double seedTime = (*iseed)->time();

       crystalToVisit.push(*iseed);  
       crystalsInCluster.push_front(*iseed);
       hitList.erase(iseed);

       while (!crystalToVisit.empty()) {

	    std::vector<int> neighbors = cal.neighborsByLevel(crystalToVisit.front()->id(),1);
	    crystalToVisit.pop();

	    for (unsigned int iv=0;iv<neighbors.size();++iv)
	    {

		//check if this neighbor crystal is in the list
		for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
		{	       
		   if ( (*it)->id() != neighbors[iv]) continue;

		   //check if timing is ok
                   if (  std::abs((*it)->time() - seedTime)  > _deltaTime) continue;

		   //check crystal energy and expension properties
		   if ((*it)->energyDep() >= _ExpandCut) crystalToVisit.push(*it);
		   if ((*it)->energyDep() >= _EnoiseCut) crystalsInCluster.push_front(*it);
 
		   //remove hit from hitList
		   hitList.erase(it);
		   break;	       
		}
	    }          
       }  

       return crystalsInCluster;
   }  





   std::list<CaloCrystalHit const*>::iterator   MakeCaloCluster::findCrystalSeed(CaloCrystalList& hitList)
   {
       
	 CaloCrystalList_iter iseed = hitList.begin();
	 double Emax((*iseed)->energyDep());       
	 double Tmin((*iseed)->time());

	 if (_caloClusterSeeding.compare("TIME") == 0) 
	 {       	  
	    for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
               if ((*it)->time() > Tmin) {iseed=it, Tmin=(*it)->time();}
	    
	    return iseed;
	 } 

	 for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
	    if ((*it)->energyDep() > Emax) {iseed=it, Emax=(*it)->energyDep();}
	 
	 return iseed;
      
   }





   CLHEP::Hep3Vector MakeCaloCluster::calculateCog(std::list<CaloCrystalHit const*> cluster, Calorimeter const & cal, int mode) 
   {
       const double Offset(4.0);

       CLHEP::Hep3Vector aVector(0,0,0);
       double sumWeights(0);    

       for (std::list<CaloCrystalHit const*>::const_iterator it = cluster.begin(); it != cluster.end(); ++it)
       {        
	 double energy = (*it)->energyDep();
	 CLHEP::Hep3Vector crystalPos = cal.crystalOrigin((*it)->id());

	 if (energy > 1e-6) {
             double weight = 1;
	     if (mode==2) weight = Offset+log(energy);
	     aVector[0] += crystalPos.x()*weight;
	     aVector[1] += crystalPos.y()*weight;
	     aVector[2] += crystalPos.z()*weight;
	     sumWeights += weight;
	 }

      }

      aVector[0] /= sumWeights;
      aVector[1] /= sumWeights;
      aVector[2] /= sumWeights;

      return aVector;

   }






}// end namespace mu2e



using mu2e::MakeCaloCluster;
DEFINE_ART_MODULE(MakeCaloCluster);
