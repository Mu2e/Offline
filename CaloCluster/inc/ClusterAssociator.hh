#ifndef CaloCluster_ClusterAssociator_HH_
#define CaloCluster_ClusterAssociator_HH_
//
// Utility class to associate proto-clusters together  
//  - low energy clusters (split-off) are associated to a single high energy cluster
//  - several high-energy clusters canbe associated together. Currently, all clusters 
//    that are linked together are regrouped into one cluster (i.e. A is linked to B, 
//    and B to C, then ABC are grouped together)

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

#include <unordered_map>
#include <map>
#include <queue>


namespace mu2e {


    class ClusterAssociator 
    {
        public:
            typedef art::Ptr<CaloCrystalHit>        CaloCrystalHitPtr;                     
            typedef std::vector<CaloCrystalHitPtr > CaloCrystalHitPtrVector;                     

            ClusterAssociator(const Calorimeter& cal): cal_(cal), associatedSplitId_(), associatedMainId_() {};              

            
            void      associateSplitOff(const CaloProtoClusterCollection&, const CaloProtoClusterCollection&,double,double);           
            void      associateMain    (const CaloProtoClusterCollection&, double, double, int);  
            double    closestDistance  (const CaloCrystalHitPtrVector&, const CaloCrystalHitPtrVector&);

            const std::vector<unsigned>& associatedMainId (unsigned i) const {return associatedMainId_.find(i)->second;}
                   unsigned              associatedSplitId(unsigned i) const {return associatedSplitId_.find(i)->second;}


        private:
            const Calorimeter& cal_;
            std::unordered_map<unsigned, unsigned>             associatedSplitId_;
            std::unordered_map<unsigned,std::vector<unsigned>> associatedMainId_;
             
    };


} 

#endif
