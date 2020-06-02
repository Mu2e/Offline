#ifndef CaloCluster_ClusterFinder_HH_
#define CaloCluster_ClusterFinder_HH_
//
// Class to find cluster of simply connected crystals
// 
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include <vector>
#include <queue>
#include <list>

namespace mu2e {


    class ClusterFinder 
    {
         public:
             typedef std::list<const CaloCrystalHit*>    CaloCrystalList;
             typedef std::vector<const CaloCrystalHit*>  CaloCrystalVec;

             ClusterFinder(const Calorimeter&, const CaloCrystalHit*, double, double, bool isOnline = false);  
             
             void                   formCluster(std::vector<CaloCrystalList>&);
             const CaloCrystalList& clusterList() const {return clusterList_;}             


         private:
             const Calorimeter*     cal_;
             const CaloCrystalHit*  crystalSeed_;
             double                 seedTime_;
             CaloCrystalList        clusterList_;
             std::queue<int>        crystalToVisit_;
             std::vector<bool>      isVisited_; 
             double                 deltaTime_; 
             double                 ExpandCut_;
 	     bool 		    isOnline_;
    };


}

#endif
