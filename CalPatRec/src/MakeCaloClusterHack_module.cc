// To do list
// Must form a calo Cluster object

// then must compute COG / direction / other crap
// then must check effciency
// then must redo everything to include splitting
// and finally get a beer for all these trouble
/*
 * MakeCaloClusterHack_module.cc (cloned from MakeCaloClusterew_module.cc)
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

#include "TROOT.h"
#include "TFolder.h"
#include "CalPatRec/inc/THackData.hh"

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
#include "CaloCluster/inc/CaloClusterer.hh"
#include "CaloCluster/inc/CaloClusterCogCalculator.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
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

bool caloCrystalHitEnergyPredicate( CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {
  return lhs->energyDep() > rhs->energyDep();
}

bool caloCrystalHitTimePredicate( CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {
  return lhs->time() < rhs->time();
}

//-----------------------------------------------------------------------------
// need to order clusters in the descending energy order, thus reverse 'less'
//-----------------------------------------------------------------------------
bool caloClusterEnergyPredicate(CaloCluster lhs, CaloCluster rhs) {
  return lhs.energyDep() > rhs.energyDep();
}


class MakeCaloClusterHack : public art::EDProducer {


   public:
           
	   typedef std::list<CaloCrystalHit const*>                  CaloCrystalList;
           typedef std::list<CaloCrystalHit const*>::iterator        CaloCrystalList_iter;
           typedef std::list<CaloCrystalHit const*>::const_iterator  CaloCrystalList_const_iter;
           
	   
	   explicit MakeCaloClusterHack(fhicl::ParameterSet const& pset) :

           // Parameters
           _diagLevel(pset.get<int>("diagLevel",0)),
           _maxFullPrint(pset.get<int>("maxFullPrint",5)),
           _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)),
           _deltaTimePlus(pset.get<double>("deltaTimePlus", 10.)),// ns
           _deltaTimeMinus(pset.get<double>("deltaTimeMinus", 10.)),// ns
           _nCryPerCluster(pset.get<int>("nCryPerCrystal", 0)),
           _EnoiseCut(pset.get<double>("EnoiseCut", 0.090)),//MeV 3 sigma noise
           _ExpandCut(pset.get<double>("ExpandCut", 0.090)),//MeV
           _EminCluster(pset.get<double>("EminCluster", 10)),//MeV
           _MinCalPatRecSeedEnergy(pset.get<double>("MinCalPatRecSeedEnergy", 60.)),//MeV
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

		   fHackData = new THackData("HackData","Hack Data");
		   gROOT->GetRootFolder()->Add(fHackData);

		   // hack = (THackData*) gROOT->GetRootFolder()->FindObject("HackData");
           }
           
	   
	   virtual ~MakeCaloClusterHack() { }
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
           double _EminCluster;

  double _MinCalPatRecSeedEnergy;

           string _g4ModuleLabel;  // Name of the module that made these hits.
           string _caloReadoutModuleLabel;
           string _caloCrystalModuleLabel;
           string _caloClusterAlgorithm;
           string _caloClusterSeeding;
           const string _producerName;

  THackData* fHackData;

           // A category for the error logger.
           const std::string _messageCategory;

           // Give some informationation messages only on the first event.
           bool _firstEvent;

 	   void makeCaloClusters(CaloClusterCollection& caloClusters, 
                	         art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);

	   std::list<CaloCrystalHit const*> formCluster(CaloCrystalList& hitList, CaloCrystalList_iter crystalSeed, Calorimeter const & cal);
	   std::list<CaloCrystalHit const*>::iterator findCrystalSeed(CaloCrystalList& hitList);
           int findMainCluster(double time, std::vector<double>& protoClusterTimingList);
           CLHEP::Hep3Vector calculateCog(std::list<CaloCrystalHit const*> cluster, Calorimeter const & cal, int mode=1); 

   double closestDistance(std::list<CaloCrystalHit const*> cluster, std::list<CaloCrystalHit const*> cluster2, Calorimeter const & cal) ;

   };


//must include a check of the parameters here
   void MakeCaloClusterHack::beginJob(){
           
	if ( (_caloClusterSeeding.compare("TIME") != 0 ) &&  (_caloClusterSeeding.compare("ENERGY") != 0 ) )	     
	   throw cet::exception("CaloCluster") << "caloClusterSeeding must be either \"energy\" or \"time\" \n";

	   cout<<"selected clustering algorithm--> "<<_caloClusterAlgorithm <<", seeded by "<< _caloClusterSeeding<<endl;
	   //        cout << "Diaglevel: "
           //                        << _diagLevel << " "
           //                        << _maxFullPrint
           //                        << endl;

   }



   void MakeCaloClusterHack::produce(art::Event& event) {
     int           ncl;
     //     double        max_energy(70.);
     CaloCluster   *cl, *cl_max(NULL);

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
     
     ncl = caloClusters->size();
     if (ncl > 0) {
//-----------------------------------------------------------------------------
// 2013-10-13 P.Murat: make sure clusters are sorted in energy, then store time 
//                     of the most energetic cluster in a TDataHack object
//-----------------------------------------------------------------------------
       std::sort(caloClusters->begin(),caloClusters->end(),mu2e::caloClusterEnergyPredicate);

       float emax = _MinCalPatRecSeedEnergy;

       for (int i=0; i<ncl; i++) {
	 cl =  &caloClusters->at(i);
	 if ((cl->time() > 400.) && (cl->energyDep() > emax)) {
	   cl_max     = cl;
	   emax       = cl->energyDep();
	   break;
	 }
       }

       if (cl_max != NULL) {
	 CaloClusterTools cluTool(*cl_max);
	 fHackData->SetClusterT0(cluTool.timeFasterCrystal());//cl_max->time());
	 // Store x,y,z postion of the cluster in the tracker system
	 fHackData->SetClusterZ(cl_max->cog3Vector().z() - 10200.);
	 fHackData->fData[0] = cl_max->cog3Vector().x()+3904.;
	 fHackData->fData[1] = cl_max->cog3Vector().y();
       }
       else {
	 fHackData->SetClusterT0(-1.);
	 fHackData->SetClusterZ(-9999.);
       }
     }
     
     event.put(std::move(caloClusters), _producerName);

   }
   
   
		     
   
   
   void MakeCaloClusterHack::makeCaloClusters(CaloClusterCollection& caloClusters,
                                      art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle) {


     
          CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
          Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  
  
          //Get a working copy of the CaloCrystalHits 
	  CaloCrystalList caloCrystalHitsWork;
	  for ( CaloCrystalHitCollection::const_iterator i=caloCrystalHits.begin(); i!=caloCrystalHits.end(); ++i )
	    caloCrystalHitsWork.push_back( &(*i));




	  // Sort crystals by energy/time -> seed of new cluster is always the first of the current list
	  if (_caloClusterSeeding.compare("TIME") == 0)  caloCrystalHitsWork.sort(mu2e::caloCrystalHitTimePredicate);
          else                                           caloCrystalHitsWork.sort(mu2e::caloCrystalHitEnergyPredicate);
	  



          // First, find clusters with seed energy above EminCluster
	  std::vector< CaloCrystalList > protoClusterList;
	  std::vector<double> protoClusterTimingList;
	  	 	      
	  while( ! caloCrystalHitsWork.empty() ){       
	    
	      //CaloCrystalList_iter crystalSeed = findCrystalSeed(caloCrystalHitsWork);	    
	      CaloCrystalList_iter crystalSeed = caloCrystalHitsWork.begin();
	      if ( (*crystalSeed)->energyDep()<_EminCluster ) break;
	      CaloCrystalList  cluster = formCluster(caloCrystalHitsWork,crystalSeed,cal);

	      cluster.sort(mu2e::caloCrystalHitEnergyPredicate);
	      protoClusterList.push_back(cluster); 
	      protoClusterTimingList.push_back( (*cluster.begin())->time() );

	      //std::cout<<"Cluster "<<protoClusterList.size()<<" with "<<cluster.size()<<" crystals "<<std::endl;
	      //for(CaloCrystalList_const_iter i = cluster.begin(); i!= cluster.end();++i)
	      //std::cout<<(*i)->id()<<" "<<(*i)->energyDep()<<" "<<(*i)->time()<<" "<<cal.crystalOrigin((*i)->id())<<std::endl;

	  }



	  // Find clusters split-off, keep only those with time compatible with big cluster
	  std::vector< CaloCrystalList > protoSplitList;
	  while( ! caloCrystalHitsWork.empty() ){       
	    
	      CaloCrystalList_iter crystalSeed = caloCrystalHitsWork.begin();
	      if ( (*crystalSeed)->energyDep() < 2. ) break;
	      double crystalTime = (*crystalSeed)->time();

	      // find cluster if the crystal timing is compatible with the list of energetic clusters
	      if ( findMainCluster(crystalTime, protoClusterTimingList) > -1) {

		  CaloCrystalList cluster = formCluster(caloCrystalHitsWork,crystalSeed,cal);
		  cluster.sort(mu2e::caloCrystalHitEnergyPredicate);
		  protoSplitList.push_back(cluster); 

  	          //std::cout<<"Cluster sec "<<protoClusterList.size()<<std::endl;
	          //for(CaloCrystalList_const_iter i = cluster.begin(); i!= cluster.end();++i)
	          //std::cout<<(*i)->id()<<" "<<(*i)->energyDep()<<" "<<(*i)->time()<<" "<<cal.crystalOrigin((*i)->id())<<std::endl;

	      } 
	      else 
	      { 
		caloCrystalHitsWork.erase(crystalSeed); 
	      }

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
	      
	      CLHEP::Hep3Vector cog = calculateCog(thisList,cal,1);
	      CaloCluster caloCluster(isection,averageTime,totalEnergy,caloCrystalHitPtrVector);	      
      	      caloCluster.SetCog3Vector(cog);
	      caloClusters.push_back(caloCluster);

	      //std::cout<<" Cluster cog "<<cog)<<std::endl;
	      
	  }


          //Very approximative method to produce the cluster, enough for now
          for (std::vector< CaloCrystalList >::const_iterator it = protoSplitList.begin(); it!=protoSplitList.end();++it){
	     
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
              
	      int mainIndex = findMainCluster((*thisList.begin())->time(), protoClusterTimingList);

	      double minDist;
	      if (mainIndex >= 0) minDist = closestDistance(thisList, protoClusterList[mainIndex],cal);
	      else                minDist = 1.e6;
   	      
	      CLHEP::Hep3Vector cog = calculateCog(thisList,cal,1);

	      CaloCluster caloCluster(isection,averageTime,totalEnergy,caloCrystalHitPtrVector);	      
      	      caloCluster.SetCog3Vector(cog);
      	      caloCluster.SetDistance(minDist);	      
      	      caloCluster.SetParentId(mainIndex);	      
	      caloClusters.push_back(caloCluster);

	      //std::cout<<" Cluster cog "<<cog<<std::endl;
	      
	  }


   }





   
   
   
   std::list<CaloCrystalHit const*> MakeCaloClusterHack::formCluster(CaloCrystalList& hitList, CaloCrystalList_iter crystalSeed, Calorimeter const & cal)
   {
       
       CaloCrystalList   crystalsInCluster;
       std::queue<CaloCrystalHit const*>  crystalToVisit;

       double seedTime   = (*crystalSeed)->time();
       crystalToVisit.push(*crystalSeed);  
       crystalsInCluster.push_front(*crystalSeed);
       hitList.erase(crystalSeed);


       while (!crystalToVisit.empty()) {
	    
	    std::vector<int> neighborsId = cal.neighborsByLevel(crystalToVisit.front()->id(),1);
	    crystalToVisit.pop();

	    //check if there are crystals in the hit list corresponding to the neighbours with consistent time
	    for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
	    {	       

	       if ( find(neighborsId.begin(),neighborsId.end(),(*it)->id()) == neighborsId.end() ) continue;
	       if ( (*it)->time() - seedTime  > _deltaTimePlus) continue;
               if (  seedTime - (*it)->time() > _deltaTimeMinus) continue;
	       
	       if ((*it)->energyDep() >= _ExpandCut) crystalToVisit.push(*it);
	       if ((*it)->energyDep() >= _EnoiseCut) crystalsInCluster.push_front(*it);

	       //remove hit from hitList (the -- is needed because erase increments the iterator)
	       hitList.erase(it--);
	    }
       }  

       return crystalsInCluster;
   }  

   std::list<CaloCrystalHit const*>::iterator  MakeCaloClusterHack::findCrystalSeed(CaloCrystalList& hitList)
   {
       
	 CaloCrystalList_iter iseed = hitList.begin();
	 double Emax((*iseed)->energyDep());       
	 double Tmin((*iseed)->time());

	 if (_caloClusterSeeding.compare("TIME") == 0) 
	 {       	  
	    for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
               if ((*it)->time() > Tmin) {iseed=it, Tmin=(*it)->time();}	    
	 } 
	 else 
	 {	    
	    for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
	       if ((*it)->energyDep() > Emax) {iseed=it, Emax=(*it)->energyDep();}
	 }

	 return iseed;      
   }

   int MakeCaloClusterHack::findMainCluster(double time, std::vector<double>& protoClusterTimingList)
   {   
      for (unsigned int it=0;it<protoClusterTimingList.size();++it)
	if ( std::abs( time-protoClusterTimingList[it] ) < _deltaTimePlus) return it;
      return -1;     
   }



   CLHEP::Hep3Vector MakeCaloClusterHack::calculateCog(std::list<CaloCrystalHit const*> cluster, Calorimeter const & cal, int mode) 
   {
       //CaloClusterCogCalculator calc;
       //CLHEP::Hep3Vector newResult = calc.calculateCog(cluster,CaloClusterCogCalculator::log);
       
       
       
       const double Offset(4.0);

       CLHEP::Hep3Vector aVector(0,0,0);
       double sumWeights(0);    

       for (std::list<CaloCrystalHit const*>::const_iterator it = cluster.begin(); it != cluster.end(); ++it)
       {        
	 double energy = (*it)->energyDep();
	 CLHEP::Hep3Vector crystalPos = cal.crystalOrigin((*it)->id());

	 if (energy > 1e-6) {
             double weight = energy;
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




   double MakeCaloClusterHack::closestDistance(std::list<CaloCrystalHit const*> cluster, std::list<CaloCrystalHit const*> cluster2, Calorimeter const & cal) 
   {
       
       double minDistance(1e6);    
       for (CaloCrystalList_const_iter it = cluster.begin(); it != cluster.end(); ++it) {        
	 
	 CLHEP::Hep3Vector crystalPos = cal.crystalOrigin((*it)->id());

	 for (CaloCrystalList_const_iter it2 = cluster2.begin(); it2 != cluster2.end(); ++it2)
	 {
		 CLHEP::Hep3Vector crystalPos2 = cal.crystalOrigin((*it2)->id());
		 double dist = (crystalPos-crystalPos2).mag();
		 if (dist<minDistance) minDistance = dist;
	 }	 
       }
       
       return minDistance;
   }

}// end namespace mu2e



using mu2e::MakeCaloClusterHack;
DEFINE_ART_MODULE(MakeCaloClusterHack);
