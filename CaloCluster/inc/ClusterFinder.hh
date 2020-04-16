//
// Original author B.Echenard
//

#ifndef CaloCluster_ClusterFinder_HH_
#define CaloCluster_ClusterFinder_HH_


// Mu2e includes
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"


// C++ includes
#include <vector>
#include <queue>
#include <list>



namespace mu2e {


    class ClusterFinder {


         public:
             
             typedef std::list<CaloCrystalHit const*>    CaloCrystalList;
             typedef std::vector<CaloCrystalHit const*>  CaloCrystalVec;


             ClusterFinder(Calorimeter const& cal, CaloCrystalHit const* crystalSeed, double deltaTime, double ExpandCut);              
        
             
             ClusterFinder(Calorimeter const& cal, CaloCrystalHit const* crystalSeed, double deltaTime, double ExpandCut, bool isOnline);  

	     ~ClusterFinder(){};            
             
             CaloCrystalList const& clusterList()  const {return clusterList_;}
             
             void formCluster(std::vector<CaloCrystalList>& idHitVec);



         private:
             
             Calorimeter const*     cal_;
             CaloCrystalHit const*  crystalSeed_;
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
