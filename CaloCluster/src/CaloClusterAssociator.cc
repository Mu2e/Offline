//
// Utility class to associate proto-clusters together, separated into 
// two medium energy clusters or a cluster and a split-off (might want to treat them differently)
// 
// Original author B. Echenard
//

// Mu2e includes
#include "CaloCluster/inc/CaloClusterAssociator.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"


// C++ includes
#include <iostream>


namespace mu2e {

       
       void CaloClusterAssociator::associateSplitOff(CaloProtoClusterCollection const& mainClusterColl, 
	                                             CaloProtoClusterCollection const& splitClusterColl,
						     double deltaTimePlus, double deltaTimeMinus, double maxDist)  
       { 
          
	  _associatedSplitId.clear();

          int jStart(0);
	  for (unsigned int i=0;i<splitClusterColl.size(); ++i)
	  {	   	     
	      double minDist(1e6);
	      int jmin(-1);

	      double timeSmall                            = splitClusterColl[i].time();
	      CaloCrystalHitPtrVector const& splitHitsPtr = splitClusterColl[i].caloCrystalHitsPtrVector();	     
	      CaloCrystalHitPtr       const& hitSmall     = *(splitHitsPtr.begin());	     


	      for (unsigned int j=jStart;j<mainClusterColl.size(); ++j)
	      {
		   if (timeSmall - mainClusterColl[j].time() > deltaTimePlus)   {jStart = j; continue;}
		   if (mainClusterColl[j].time() - timeSmall > deltaTimeMinus)  break;

		   CaloCrystalHitPtrVector const& mainHitsPtr = mainClusterColl[j].caloCrystalHitsPtrVector();	     
	           CaloCrystalHitPtr       const& hitMain     = *(mainHitsPtr.begin());

		   CLHEP::Hep3Vector crystalPos1 = _cal.crystalOrigin(hitMain->id());
		   CLHEP::Hep3Vector crystalPos2 = _cal.crystalOrigin(hitSmall->id());
		   double dist0 = (crystalPos1-crystalPos2).mag();

		   if (dist0 > 2.5*maxDist) continue;

		   double dist =  closestDistance(mainHitsPtr,splitHitsPtr);
		   if (dist < maxDist && dist < minDist) { minDist=dist; jmin=j; } 
	      }	  
	                    
	      _associatedSplitId[i]   = jmin;
	  }
	
       } 
 



       //----------------------------------------------------------------------------------------------------------
       void CaloClusterAssociator::associateMain(CaloProtoClusterCollection const& clusterColl, double deltaTime, double maxDist)  
       { 
          
	  _associatedMainId.clear();

          for (unsigned int i=0;i<clusterColl.size(); ++i)
	  {	   	     
	      double minDist(1e6);
	      int jmin(-1);

	      CaloCrystalHitPtrVector const& FirstHitsPtr = clusterColl[i].caloCrystalHitsPtrVector();	     

	      for (unsigned int j=i+1;j<clusterColl.size();++j)
	      {
		 if (clusterColl[j].time() - clusterColl[i].time() > deltaTime) break;
		 
		 CaloCrystalHitPtrVector const& SecondHitsPtr = clusterColl[j].caloCrystalHitsPtrVector();	     

		 double dist = closestDistance(FirstHitsPtr,SecondHitsPtr);
		 if (dist < minDist) {minDist=dist; jmin=j;} 
	      }	  
              
	      _associatedMainId[i] = jmin;
	  }
	    	    
       } 
 


 
       //----------------------------------------------------------------------------------------------------------
       double CaloClusterAssociator::closestDistance(CaloCrystalHitPtrVector const& cluster, CaloCrystalHitPtrVector const& cluster2)
       {
	  double minDistance(1e6);    
	  for (auto const& hit : cluster)
	  {        
	      CLHEP::Hep3Vector crystalPos = _cal.crystalOrigin(hit->id());

	      for (auto const& hit2 : cluster2)
	      {
		 double dist = (crystalPos - _cal.crystalOrigin(hit2->id())).mag();
		 if (dist<minDistance) minDistance = dist;
	      }	 
	  }

	 return minDistance;
       }


}
