#include "CaloCluster/inc/ClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"

#include <iostream>
#include <vector>
#include <algorithm>


namespace mu2e {

	ClusterFinder::ClusterFinder(const Calorimeter& cal, const CaloCrystalHit* crystalSeed, double deltaTime, double ExpandCut, bool isOnline) : 
	  cal_(&cal), crystalSeed_(crystalSeed), seedTime_(crystalSeed->time()), clusterList_(), crystalToVisit_(), isVisited_(cal.nCrystal()),
	  deltaTime_(deltaTime), ExpandCut_(ExpandCut), isOnline_(isOnline)
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

		 std::vector<int>  neighborsId = cal_->crystal(visitId).neighbors();
                 
                 //Why??? That breaks the connectivity
                 if (isOnline_) neighborsId.insert(neighborsId.end(), cal_->nextNeighbors(visitId).begin(), cal_->nextNeighbors(visitId).end());
		 
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
			 else ++it;
		     } 
                 }

                 crystalToVisit_.pop();                 
            }

	    // make sure to sort proto0cluster by energy
	    clusterList_.sort([](const CaloCrystalHit* lhs, const CaloCrystalHit* rhs) {return lhs->energyDep() > rhs->energyDep();});
       } 

}








