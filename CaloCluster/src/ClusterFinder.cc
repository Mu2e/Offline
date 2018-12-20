//
// Class to find cluster of simply connected crystals
// 
// Original author B. Echenard
//
// Note: there are few places where a continue could be replaced by a break if the crystal are time ordered, but 
//       the performance gain is so low that it outweighs the risk of forgeting to time order the crystal hits.
//       

#include "CaloCluster/inc/ClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include <iostream>
#include <vector>
#include <algorithm>


namespace mu2e {


       ClusterFinder::ClusterFinder(Calorimeter const& cal, CaloCrystalHit const* crystalSeed, double deltaTime, double ExpandCut) : 
          cal_(&cal), crystalSeed_(crystalSeed),seedTime_(crystalSeed->time()), clusterList_(), crystalToVisit_(), isVisited_(cal.nCrystal()),
          deltaTime_(deltaTime), ExpandCut_(ExpandCut)
       {}
       
       


       void ClusterFinder::formCluster(std::vector<CaloCrystalList>& idHitVec)  
       { 

            clusterList_.clear();            
            clusterList_.push_front(crystalSeed_);
            crystalToVisit_.push(crystalSeed_->id());  

                        
            CaloCrystalList& liste = idHitVec[crystalSeed_->id()];
            liste.erase(std::find(liste.begin(), liste.end(), crystalSeed_));


            while (!crystalToVisit_.empty())
            {            
                 int visitId         = crystalToVisit_.front();
                 isVisited_[visitId] = 1;

                 std::vector<int> const& neighborsId = cal_->crystal(visitId).neighbors();
                 for (auto& iId : neighborsId)
                 {               
                     if (isVisited_[iId]) continue;
                     isVisited_[iId]=1;


                     CaloCrystalList& list = idHitVec[iId];
                     auto it=list.begin();
                     while(it != list.end())
                     {
                         CaloCrystalHit const* hit = *it;
                         if (std::abs(hit->time() - seedTime_) < deltaTime_)
                         { 
                            if (hit->energyDep() > ExpandCut_) crystalToVisit_.push(iId);
                            clusterList_.push_front(hit);
                            it = list.erase(it);   
                         } 
                         else {++it;}
                     } 
                     
                 }
                                       
                 crystalToVisit_.pop();                 
            }
            
           // sort proto-clustres, even if they are sorted in the cluster module, in case somebody 
           // uses the proto-clusters instead of clusters he will get the same behaviour
           clusterList_.sort([] (CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();} );
                        
       } 

}








