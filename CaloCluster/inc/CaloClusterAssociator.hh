//
// Original author B.Echenard
//

#ifndef CaloCluster_CaloClusterAssociator_HH_
#define CaloCluster_CaloClusterAssociator_HH_


// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"


// C++ includes
#include <unordered_map>
#include <map>
#include <queue>



namespace mu2e {


    class CaloClusterAssociator {


	 public:

 	     typedef art::Ptr<CaloCrystalHit>        CaloCrystalHitPtr;                     
 	     typedef std::vector<CaloCrystalHitPtr > CaloCrystalHitPtrVector;                     

	     
	     
	     CaloClusterAssociator(Calorimeter const& cal): _cal(cal),_associatedMainId(), _associatedSplitId() {};              

	     ~CaloClusterAssociator(){};
	     	     


             void      associateSplitOff(CaloProtoClusterCollection const& mainClusterColl, 
	                                 CaloProtoClusterCollection const& splitClusterColl, 
	                                 double deltaTimePlus, double deltaTimeMinus,double maxDist);
             void      associateMain(CaloProtoClusterCollection const& clusterColl, 
				     double deltaTimePlus, double deltaTimeMinus, double maxDist);  
					 
	     double    closestDistance(CaloCrystalHitPtrVector const& cluster, CaloCrystalHitPtrVector const& cluster2);

             unsigned  associatedMainId(unsigned int i)    const {return _associatedMainId.find(i)->second;}
             unsigned  associatedSplitId(unsigned int i)   const {return _associatedSplitId.find(i)->second;}
	     

	 private:
             
	     Calorimeter const&                   _cal;
	     
	     std::map<unsigned int,int>  _associatedMainId;
	     std::map<unsigned int,int>  _associatedSplitId;
	     
    };


} 

#endif
