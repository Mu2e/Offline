#include "CaloCluster/inc/ClusterAssociator.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

#include <iostream>


namespace mu2e {

       
     void ClusterAssociator::associateSplitOff(const CaloProtoClusterCollection& mainClusterColl, 
                                               const CaloProtoClusterCollection& splitClusterColl,
                                               double deltaTime, double maxDist)  
     { 

        associatedSplitId_.clear();

        int jStart(0);
        for (unsigned i=0;i<splitClusterColl.size(); ++i)
        {                        
            double minDist(1e6);
            int jmin(-1);

            double timeSmall                            = splitClusterColl[i].time();
            const CaloCrystalHitPtrVector& splitHitsPtr = splitClusterColl[i].caloCrystalHitsPtrVector();             
            const CaloCrystalHitPtr&       hitSmall     = *(splitHitsPtr.begin());             


            for (unsigned j=jStart;j<mainClusterColl.size(); ++j)
            {
                 if (timeSmall - mainClusterColl[j].time() > deltaTime) {jStart = j; continue;}
                 if (mainClusterColl[j].time() - timeSmall > deltaTime) break;

                 const CaloCrystalHitPtrVector& mainHitsPtr = mainClusterColl[j].caloCrystalHitsPtrVector();             
                 const CaloCrystalHitPtr&       hitMain     = *(mainHitsPtr.begin());

                 CLHEP::Hep3Vector crystalPos1 = cal_.crystal(hitMain->id()).position();
                 CLHEP::Hep3Vector crystalPos2 = cal_.crystal(hitSmall->id()).position();
                 double dist0 = (crystalPos1-crystalPos2).mag();

                 if (dist0 > 2.5*maxDist) continue;

                 double dist =  closestDistance(mainHitsPtr,splitHitsPtr);
                 if (dist < maxDist && dist < minDist) {minDist=dist; jmin=j;} 
            }          
            associatedSplitId_[i] = jmin;
        }

     } 


     //------------------------------------------------------------------------------------------------------------------------
     //if you want to turn-off chain linking (see note above), uncomment line 1 and comment after line 2
     void ClusterAssociator::associateMain(const CaloProtoClusterCollection & clusterColl, double deltaTime, double maxDist, int strategy)  
     { 

        std::vector<int>              isAssociatedTo(clusterColl.size(),-1);
        std::vector<std::vector<int>> associatedId(clusterColl.size());

        associatedMainId_.clear();

        for (unsigned i=0;i<clusterColl.size(); ++i)
        {                                       
            associatedMainId_[i].clear();
            CaloCrystalHitPtrVector const& FirstHitsPtr = clusterColl[i].caloCrystalHitsPtrVector();             

            for (unsigned j=i+1;j<clusterColl.size();++j)
            {
                 if (isAssociatedTo[j] > -1) continue;

		 double dtime = clusterColl[j].time() - clusterColl[i].time();                  
		 if (dtime > deltaTime) break;

                 CaloCrystalHitPtrVector const& SecondHitsPtr = clusterColl[j].caloCrystalHitsPtrVector();             
                 double dist        = closestDistance(FirstHitsPtr,SecondHitsPtr);
		 double estimatedDt = dist / 300.0;
		 bool   sameDisk    = cal_.crystal(FirstHitsPtr[0]->id()).diskId() == cal_.crystal(SecondHitsPtr[0]->id()).diskId();

		 if (strategy == 2 && std::abs(estimatedDt-dtime) > deltaTime) continue;

		 if (strategy == 1 && !sameDisk ) continue;
		 if (strategy == 1 && sameDisk  && std::abs(estimatedDt-dtime) > deltaTime ) continue;		 

		 if (strategy == 0 &&  (dtime > deltaTime || dist > maxDist))  continue;

                 isAssociatedTo[j] = i;
                 associatedId[i].push_back(j);
                 //associatedMainId_[i].push_back(j); //line 1
            }          
        }

        //line 2
        for (unsigned i=0;i<associatedId.size();++i)
        {               
             std::set<int> neighbors;
             std::queue<int> list;
             for (int id : associatedId[i]){neighbors.insert(id); list.push(id);}

             while (!list.empty())
             {
                 int nextId = list.front();
                 for (auto id : associatedId[nextId]) {neighbors.insert(id); list.push(id);}
                 associatedId[nextId].clear();  
                 list.pop();               
             }
             for (auto id : neighbors) associatedMainId_[i].push_back(id);
        }


     } 





     //----------------------------------------------------------------------------------------------------------
     double ClusterAssociator::closestDistance(const CaloCrystalHitPtrVector& cluster, const CaloCrystalHitPtrVector& cluster2)
     {
        double minDistance(1e6);    
        for (const auto& hit : cluster)
        {        
            CLHEP::Hep3Vector crystalPos = cal_.crystal(hit->id()).position();

            for (const auto& hit2 : cluster2)
            {
               double dist = (crystalPos - cal_.crystal(hit2->id()).position()).mag();
               if (dist<minDistance) minDistance = dist;
            }         
        }

       return minDistance;
     }


}
