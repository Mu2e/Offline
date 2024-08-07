#include "Offline/CaloCluster/inc/ClusterFinder.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include <iostream>
#include <vector>
#include <algorithm>


namespace mu2e {

        ClusterFinder::ClusterFinder(const Calorimeter& cal, const CaloHit* crystalSeed, double deltaTime, double ExpandCut, bool addSecondRing) :
          cal_(&cal), crystalSeed_(crystalSeed), seedTime_(crystalSeed->time()), clusterList_(), crystalToVisit_(), isVisited_(cal.nCrystals()),
          deltaTime_(deltaTime), ExpandCut_(ExpandCut), addSecondRing_(addSecondRing)
        {}


        void ClusterFinder::formCluster(std::vector<CaloCrystalList>& idHitVec)
        {
            std::fill(isVisited_.begin(), isVisited_.end(), false);

            clusterList_.clear();
            clusterList_.push_front(crystalSeed_);
            crystalToVisit_.push(crystalSeed_->crystalID());

            CaloCrystalList& liste = idHitVec[crystalSeed_->crystalID()];
            liste.erase(std::find(liste.begin(), liste.end(), crystalSeed_));

            while (!crystalToVisit_.empty())
            {
                 int visitId         = crystalToVisit_.front();
                 isVisited_[visitId] = true;

                 std::vector<int>  neighborsId = cal_->crystal(visitId).neighbors();
                 if (addSecondRing_) neighborsId.insert(neighborsId.end(), cal_->nextNeighbors(visitId).begin(), cal_->nextNeighbors(visitId).end());

                 for (auto& iId : neighborsId)
                 {
                     if (isVisited_[iId]) continue;
                     isVisited_[iId] = true;
                     CaloCrystalList& list = idHitVec[iId];

                     auto it=list.begin();
                     while(it != list.end())
                     {
                         CaloHit const* hit = *it;
                         if (std::abs(hit->time() - seedTime_) < deltaTime_)
                         {
                             if (hit->energyDep() > ExpandCut_) crystalToVisit_.push(iId);
                             clusterList_.push_front(hit);
                             it = list.erase(it);
                         }
                         else ++it;
                     }
                 }

                 crystalToVisit_.pop();
            }

            // make sure to sort proto0cluster by energy
            clusterList_.sort([](const CaloHit* lhs, const CaloHit* rhs) {return lhs->energyDep() > rhs->energyDep();});
       }

}








