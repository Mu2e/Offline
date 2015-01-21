// To do list
// 2014-12-244 P.Murat: back to git commit ID 6310127

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
#include "CaloCluster/inc/CaloClusterUtilities.hh"

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
  
  class MakeCaloClusterNew : public art::EDProducer {


  public:
           
    typedef std::list<CaloCrystalHit const*>                  CaloCrystalList;
    typedef std::list<CaloCrystalHit const*>::iterator        CaloCrystalList_iter;
    typedef std::list<CaloCrystalHit const*>::const_iterator  CaloCrystalList_const_iter;
           
	   
    explicit MakeCaloClusterNew(fhicl::ParameterSet const& pset) :
	     
      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      //           _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)),
      _deltaTimePlus (pset.get<double>("deltaTimePlus" ,  5.)),// ns
      _deltaTimeMinus(pset.get<double>("deltaTimeMinus",  5.)),// ns
      _nCryPerCluster(pset.get<int>("nCryPerCrystal", 0)),
      _EnoiseCut(pset.get<double>("EnoiseCut", 0.090)),//MeV 3 sigma noise
      _ExpandCut(pset.get<double>("ExpandCut", 0.090)),//MeV
      _EminCluster  (pset.get<double>("EminCluster"  , 5.)),//MeV
      _EminSplitSeed(pset.get<double>("EminSplitSeed", 2.)),//MeV
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
      _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
      _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "ENERGY")),
      //      _producerName("Algo"+TOUpper(_caloClusterAlgorithm)+"SeededBy"+TOUpper(_caloClusterSeeding)),
      _producerName(""),
      _messageCategory("HITS"),
      _firstEvent(true)
    {
      TOUpper(_caloClusterSeeding);
					// Tell the framework what we make.
      produces<CaloClusterCollection>(_producerName);
    }
           
	   
    virtual ~MakeCaloClusterNew() { }

    virtual void beginJob();
    virtual void beginRun(art::Run&   aRun   );
    virtual void produce (art::Event& anEvent);

  private:
   
 
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the calorimeter StepPoint collection
    std::string _caloStepPoints;

    // Parameters
    //           double _minimumEnergy;  // minimum energy deposition of G4 step
    double _deltaTimePlus;
    double _deltaTimeMinus;
    double _nCryPerCluster;
    double _EnoiseCut;
    double _ExpandCut;
    double _EminCluster;
    double _EminSplitSeed;
    string _g4ModuleLabel;  // Name of the module that made these hits.
    string _caloReadoutModuleLabel;
    string _caloCrystalModuleLabel;
    string _caloClusterAlgorithm;
    string _caloClusterSeeding;
    const string _producerName;

    const Calorimeter* _cal;

    // A category for the error logger.
    const std::string _messageCategory;
    
    // Give some informationation messages only on the first event.
    bool _firstEvent;
    
    void makeCaloClusters(CaloClusterCollection& caloClusters, 
			  art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle);
    
    std::list<CaloCrystalHit const*> formCluster(CaloCrystalList& hitList, CaloCrystalList_iter crystalSeed);
    
    std::list<CaloCrystalHit const*>::iterator findCrystalSeed(CaloCrystalList& hitList);
    
    int findMainCluster(double time, std::vector<double>& protoClusterTimingList);
    
    int findMainClusterA(int DiskID, double Time, CaloClusterCollection& ClusterList, int NClusters);
    
    CLHEP::Hep3Vector calculateCog(std::list<CaloCrystalHit const*> cluster, int mode=1); 
    
    double closestDistance(std::list<CaloCrystalHit const*> cluster, std::list<CaloCrystalHit const*> cluster2) ;
    
  };


  // must include a check of the parameters here
  void MakeCaloClusterNew::beginJob() {
           
    if ( (_caloClusterSeeding.compare("TIME") != 0 ) &&  (_caloClusterSeeding.compare("ENERGY") != 0 ) )	     
      throw cet::exception("CaloCluster") << "caloClusterSeeding must be either \"energy\" or \"time\" \n";
    
    cout<<"selected clustering algorithm--> "<<_caloClusterAlgorithm <<", seeded by "<< _caloClusterSeeding<<endl;
    //        cout << "Diaglevel: "
    //                        << _diagLevel << " "
    //                        << _maxFullPrint
    //                        << endl;
    
  }


//-----------------------------------------------------------------------------
// update calorimeter handle at begin run - calibrations etc may change....
//-----------------------------------------------------------------------------
  void MakeCaloClusterNew::beginRun(art::Run& aRun) {
    GeomHandle<Calorimeter>   ch;
    _cal = ch.get();
  }


//-----------------------------------------------------------------------------
   void MakeCaloClusterNew::produce(art::Event& event) {

     //Get handles to calorimeter crystal hits
     art::Handle<CaloCrystalHitCollection> caloCrystalHitsHandle;
     event.getByLabel(_caloCrystalModuleLabel, caloCrystalHitsHandle);
     if ( !caloCrystalHitsHandle.isValid()) return;
     
     //Create a new CaloCluster collection and fill it
     unique_ptr<CaloClusterCollection> caloClusters(new CaloClusterCollection);
     makeCaloClusters(*caloClusters,caloCrystalHitsHandle);
     
     event.put(std::move(caloClusters), _producerName);
   }
   
//-----------------------------------------------------------------------------
   void MakeCaloClusterNew::makeCaloClusters(CaloClusterCollection&                       caloClusters,
					     art::Handle<CaloCrystalHitCollection> const& caloCrystalHitsHandle) {

     int                            ncl, isection, seed_disk_id;
     double                         totalEnergy,  averageTime, seed_time;
     std::vector< CaloCrystalList > protoSplitList;
     const CaloCrystalHit*          seed;
     
     CaloCrystalHitCollection const& caloCrystalHits(*caloCrystalHitsHandle);
     CaloCrystalHit const* caloCrystalHitBase = &caloCrystalHits.front();

       //Get a working copy of the CaloCrystalHits 
     CaloCrystalList caloCrystalHitsWork;
     for ( CaloCrystalHitCollection::const_iterator i=caloCrystalHits.begin(); i!=caloCrystalHits.end(); ++i )
       caloCrystalHitsWork.push_back( &(*i));

     // Sort crystals by energy/time -> seed of new cluster is always the first of the current list
     if (_caloClusterSeeding.compare("TIME") == 0)  caloCrystalHitsWork.sort(mu2e::caloCrystalHitTimePredicate);
     else                                           caloCrystalHitsWork.sort(mu2e::caloCrystalHitEnergyPredicate);
	  
     // First, find clusters with seed energy above _EminCluster
     std::vector< CaloCrystalList > protoClusterList;
     //     std::vector<double> protoClusterTimingList;
     
     while( ! caloCrystalHitsWork.empty() ){       
       CaloCrystalList_iter crystalSeed = caloCrystalHitsWork.begin();
       if ( (*crystalSeed)->energyDep() < _EminCluster ) break;
       CaloCrystalList  cluster = formCluster(caloCrystalHitsWork,crystalSeed);

       cluster.sort(mu2e::caloCrystalHitEnergyPredicate);
       protoClusterList.push_back(cluster); 
       //       protoClusterTimingList.push_back( (*cluster.begin())->time() );
       
       //std::cout<<"Cluster "<<protoClusterList.size()<<" with "<<cluster.size()<<" crystals "<<std::endl;
       //for(CaloCrystalList_const_iter i = cluster.begin(); i!= cluster.end();++i)
       //std::cout<<(*i)->id()<<" "<<(*i)->energyDep()<<" "<<(*i)->time()<<" "<<_cal->crystalOrigin((*i)->id())<<std::endl;
     }
//-----------------------------------------------------------------------------
// make clusters from the main proto-list
// P.Murat: this loop  needs to be merged with the previous one
//-----------------------------------------------------------------------------
     for (std::vector< CaloCrystalList >::const_iterator it = protoClusterList.begin(); it!=protoClusterList.end();++it){
       isection    = -1;
       totalEnergy = 0;
       averageTime = 0;

       CaloCrystalList thisList = (*it);
              
       std::vector< art::Ptr< CaloCrystalHit> > caloCrystalHitPtrVector;

       for (CaloCrystalList_const_iter il = thisList.begin(); il !=thisList.end(); ++il){
	      
	 CaloCrystalHit const* hit = *il;
	 totalEnergy += hit->energyDep();
	 averageTime += hit->time();
	 if (isection==-1) isection = _cal->caloSectionId(hit->id());
              
	 size_t idx = ( hit - caloCrystalHitBase );
	 caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
		
       }
       averageTime /= float(thisList.size());
	      
       CLHEP::Hep3Vector cog = calculateCog(thisList,1);
       CaloCluster caloCluster(isection,averageTime,totalEnergy,caloCrystalHitPtrVector);	      
       caloCluster.SetCog3Vector(cog);
       caloClusters.push_back(caloCluster);
     }
//-----------------------------------------------------------------------------
// 2014-04-15 P.Murat: make sure main clusters are sorted in energy
//            I think 'sort' breaks if given an empty vector
//-----------------------------------------------------------------------------
     ncl = caloClusters.size();
     if (ncl > 0) {
       std::sort(caloClusters.begin(),caloClusters.end(),mu2e::caloClusterEnergyPredicate);
     }
//-----------------------------------------------------------------------------
// next: find clusters split-off, keep only those with time compatible with big clusters
//-----------------------------------------------------------------------------
     while( ! caloCrystalHitsWork.empty() ) {       
	    
       CaloCrystalList_iter crystalSeed = caloCrystalHitsWork.begin();
       seed = (*crystalSeed);

       if (seed->energyDep() < _EminSplitSeed) break;
       seed_time    = seed->time();
       seed_disk_id = _cal->caloSectionId(seed->id());

       // find cluster if the crystal timing is compatible with the list of energetic clusters
       if ( findMainClusterA(seed_disk_id,seed_time,caloClusters,ncl) > -1) {
	 CaloCrystalList cluster = formCluster(caloCrystalHitsWork,crystalSeed);
	 cluster.sort(mu2e::caloCrystalHitEnergyPredicate);
	 protoSplitList.push_back(cluster); 
       } 
       else { 
	 caloCrystalHitsWork.erase(crystalSeed); 
       }
     }
//-----------------------------------------------------------------------------
// finally, make "split" clusters from the proto-list to their final destination
// P.Murat: this loop  needs to be merged with the previous one
//-----------------------------------------------------------------------------
     for (std::vector< CaloCrystalList >::const_iterator it = protoSplitList.begin(); it!=protoSplitList.end();++it){
       isection    = -1;
       totalEnergy = 0;
       averageTime = 0;

       CaloCrystalList thisList = (*it);

       seed          = (*thisList.begin());
       seed_time     = seed->time();
       seed_disk_id  = _cal->caloSectionId(seed->id());
              
       std::vector< art::Ptr< CaloCrystalHit> > caloCrystalHitPtrVector;
       
       for (CaloCrystalList_const_iter il = thisList.begin(); il !=thisList.end(); ++il) {
	 CaloCrystalHit const* hit = *il;
	 totalEnergy += hit->energyDep();
	 averageTime += hit->time();
	 if (isection==-1) isection = _cal->caloSectionId(hit->id());
         
	 size_t idx = (hit - caloCrystalHitBase);
	 caloCrystalHitPtrVector.push_back( art::Ptr<CaloCrystalHit>(caloCrystalHitsHandle,idx) );
       }
       averageTime /= float(thisList.size());
//-----------------------------------------------------------------------------
// use timing of the seed crystal to determine the 'main' cluster
//-----------------------------------------------------------------------------
//	      int mainIndex = findMainCluster((*thisList.begin())->time(), protoClusterTimingList);

       int mainIndex = findMainClusterA(seed_disk_id,seed_time,caloClusters,ncl);

       double minDist;
       if (mainIndex >= 0) minDist = closestDistance(thisList, protoClusterList[mainIndex]);
       else                minDist = 1.e6;

       CLHEP::Hep3Vector cog = calculateCog(thisList,1);
       
       CaloCluster caloCluster(isection,averageTime,totalEnergy,caloCrystalHitPtrVector);	      
       caloCluster.SetCog3Vector(cog);
       caloCluster.SetDistance(minDist);	      
       caloCluster.SetParentId(mainIndex);	      
       caloClusters.push_back(caloCluster);
     }
   }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  std::list<CaloCrystalHit const*> MakeCaloClusterNew::formCluster(CaloCrystalList&     hitList    , 
								   CaloCrystalList_iter crystalSeed) {
       
    CaloCrystalList   crystalsInCluster;
    std::queue<CaloCrystalHit const*>  crystalToVisit;

    double seedTime   = (*crystalSeed)->time();
    crystalToVisit.push(*crystalSeed);  
    crystalsInCluster.push_front(*crystalSeed);
    hitList.erase(crystalSeed);

    while (!crystalToVisit.empty()) {
	    
      std::vector<int> neighborsId = _cal->neighbors(crystalToVisit.front()->id());
      crystalToVisit.pop();

      //check if there are crystals in the hit list corresponding to the neighbours with consistent time
      for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it) {

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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  std::list<CaloCrystalHit const*>::iterator  MakeCaloClusterNew::findCrystalSeed(CaloCrystalList& hitList) {
    CaloCrystalList_iter iseed = hitList.begin();
    double Emax((*iseed)->energyDep());       
    double Tmin((*iseed)->time());
    
    if (_caloClusterSeeding.compare("TIME") == 0) {       	  
      for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
	if ((*it)->time() > Tmin) {iseed=it, Tmin=(*it)->time();}	    
    } 
    else {	    
      for (CaloCrystalList_iter it = hitList.begin(); it != hitList.end(); ++it)
	if ((*it)->energyDep() > Emax) {iseed=it, Emax=(*it)->energyDep();}
    }

    return iseed;      
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  int MakeCaloClusterNew::findMainCluster(double time, std::vector<double>& protoClusterTimingList) {   
    for (unsigned int it=0;it<protoClusterTimingList.size();++it)
      if ( std::abs( time-protoClusterTimingList[it] ) < _deltaTimePlus) return it;
    return -1;     
  }

//-----------------------------------------------------------------------------
// find "main" cluster on 'Disk' . closest to the given "split" cluster in time
// 'Time' is the time of the "split" cluster
//-----------------------------------------------------------------------------
  int MakeCaloClusterNew::findMainClusterA(int Disk, double Time, CaloClusterCollection& ClusterList, int N) { 

    double dt, dt_min;
    int    closest = -1;

    dt_min  = _deltaTimePlus;

    for (int i=0; i<N; i++) {
      CaloCluster* cl = &ClusterList.at(i); 
      if (Disk == cl->vaneId()) {
	dt = Time-cl->time();
	if (fabs(dt) < dt_min) {
	  closest = i;
	  dt_min  = fabs(dt);
	} 
      }
    }

    return closest;     
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  CLHEP::Hep3Vector MakeCaloClusterNew::calculateCog(std::list<CaloCrystalHit const*> cluster, int mode) {
    const double Offset(4.0);

    CLHEP::Hep3Vector aVector(0,0,0);
    double sumWeights(0);    

    for (std::list<CaloCrystalHit const*>::const_iterator it = cluster.begin(); it != cluster.end(); ++it) {        
      double energy = (*it)->energyDep();
      CLHEP::Hep3Vector crystalPos = _cal->crystalOrigin((*it)->id());
      
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

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  double MakeCaloClusterNew::closestDistance(std::list<CaloCrystalHit const*> cluster, std::list<CaloCrystalHit const*> cluster2) {
    double minDistance(1e6);    
    for (CaloCrystalList_const_iter it = cluster.begin(); it != cluster.end(); ++it) {        
	 
      CLHEP::Hep3Vector crystalPos = _cal->crystalOrigin((*it)->id());
      
      for (CaloCrystalList_const_iter it2 = cluster2.begin(); it2 != cluster2.end(); ++it2) {
	CLHEP::Hep3Vector crystalPos2 = _cal->crystalOrigin((*it2)->id());
	double dist = (crystalPos-crystalPos2).mag();
	if (dist<minDistance) minDistance = dist;
      }	 
    }
       
    return minDistance;
  }

}// end namespace mu2e



using mu2e::MakeCaloClusterNew;
DEFINE_ART_MODULE(MakeCaloClusterNew);
