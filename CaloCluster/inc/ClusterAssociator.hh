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

             
             
             ClusterAssociator(Calorimeter const& cal): cal_(cal), associatedSplitId_(), associatedMainId_() {};              
             ~ClusterAssociator(){};
                          

             void      associateSplitOff(CaloProtoClusterCollection const& mainClusterColl, 
                                         CaloProtoClusterCollection const& splitClusterColl, 
                                         double deltaTime,double maxDist);
             void      associateMain(CaloProtoClusterCollection const& clusterColl, double deltaTime, double maxDist, int strategy);  
                                         
             double    closestDistance(CaloCrystalHitPtrVector const& cluster, CaloCrystalHitPtrVector const& cluster2);

             std::vector<unsigned int> const& associatedMainId(unsigned int i)  const {return associatedMainId_.find(i)->second;}
             unsigned                         associatedSplitId(unsigned int i) const {return associatedSplitId_.find(i)->second;}
             

         private:
             
             Calorimeter const& cal_;
             
             std::unordered_map<unsigned int, unsigned int>               associatedSplitId_;
             std::unordered_map<unsigned int,std::vector<unsigned int> >  associatedMainId_;
             
    };


} 

#endif
