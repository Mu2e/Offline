//
// Original author B.Echenard
//

#ifndef CaloCluster_ClusterAssociator_HH_
#define CaloCluster_ClusterAssociator_HH_


// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"


// C++ includes
#include <unordered_map>
#include <map>
#include <queue>



namespace mu2e {


    class ClusterAssociator {


	 public:

 	     typedef art::Ptr<CaloCrystalHit>        CaloCrystalHitPtr;                     
 	     typedef std::vector<CaloCrystalHitPtr > CaloCrystalHitPtrVector;                     

	     
	     
	     ClusterAssociator(Calorimeter const& cal): _cal(cal), _associatedSplitId(), _associatedMainId() {};              
	     ~ClusterAssociator(){};
	     	     

             void      associateSplitOff(CaloProtoClusterCollection const& mainClusterColl, 
	                                 CaloProtoClusterCollection const& splitClusterColl, 
	                                 double deltaTimePlus, double deltaTimeMinus,double maxDist);
             void      associateMain(CaloProtoClusterCollection const& clusterColl, double deltaTime, double maxDist);  
					 
	     double    closestDistance(CaloCrystalHitPtrVector const& cluster, CaloCrystalHitPtrVector const& cluster2);

             std::vector<unsigned int> const& associatedMainId(unsigned int i)  const {return _associatedMainId.find(i)->second;}
             unsigned                         associatedSplitId(unsigned int i) const {return _associatedSplitId.find(i)->second;}
	     

	 private:
             
	     Calorimeter const& _cal;
	     
	     std::unordered_map<unsigned int, unsigned int>               _associatedSplitId;
	     std::unordered_map<unsigned int,std::vector<unsigned int> >  _associatedMainId;
	     
    };


} 

#endif
