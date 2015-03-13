//
// Utility to study the MC content of a calo cluster
// 
// If the CaloHitSimPartMC information is not available, recompute it
// 
// Original author B. Echenard
//

// Mu2e includes
#include "CaloCluster/inc/CaloClusterAssociator.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"


// C++ includes
#include <iostream>
#include <list>
#include <map>


namespace mu2e {

       
       void CaloClusterAssociator::associateSplitOff(CaloProtoClusterCollection const& mainClusterColl, 
	                                             CaloProtoClusterCollection const& splitClusterColl,
						     double deltaTimePlus, double deltaTimeMinus, double maxDist)  
       { 
          
	  _associatedSplitId.clear();

          for (unsigned int i=0;i<splitClusterColl.size(); ++i)
	  {	   	     
	      double minDist(1e6);
	      int jmin(-1);

	      CaloCrystalHitPtrVector const& splitHitsPtr = splitClusterColl[i].caloCrystalHitsPtrVector();	     
	      CaloCrystalHitPtr       const& hitSmall     = *(splitHitsPtr.begin());	     

	      for (unsigned int j=0;j<mainClusterColl.size(); ++j)
	      {
	          CaloCrystalHitPtrVector const& mainHitsPtr = mainClusterColl[j].caloCrystalHitsPtrVector();	     
	          CaloCrystalHitPtr       const& hitMain     = *(mainHitsPtr.begin());

	          if (hitSmall->time() -  hitMain->time()  > deltaTimePlus)  continue;
        	  if (hitMain->time()  -  hitSmall->time() > deltaTimeMinus) continue;

		  double dist =  closestDistance(mainHitsPtr,splitHitsPtr);

		  if (dist < maxDist && dist < minDist) {minDist=dist; jmin=j;} 
	      }	  
              
	      _associatedSplitId[i]   = jmin;
	  }
	    	    
       } 
 

       void CaloClusterAssociator::associateMain(CaloProtoClusterCollection const& clusterColl, 
						 double deltaTimePlus, double deltaTimeMinus, double maxDist)  
       { 
          
	  _associatedMainId.clear();

          for (unsigned int i=0;i<clusterColl.size(); ++i)
	  {	   	     
	      double minDist(1e6);
	      int jmin(-1);

	      CaloCrystalHitPtrVector const& FirstHitsPtr = clusterColl[i].caloCrystalHitsPtrVector();	     
	      CaloCrystalHitPtr       const& hitFirst     = *(FirstHitsPtr.begin());	     

	      for (unsigned int j=i+1;j<clusterColl.size();++j)
	      {

	         if (i==j) continue;
		 CaloCrystalHitPtrVector const& SecondHitsPtr = clusterColl[j].caloCrystalHitsPtrVector();	     
	         CaloCrystalHitPtr       const& hitSecond     = *(SecondHitsPtr.begin());

	         if (hitFirst->time()  -  hitSecond->time() > deltaTimePlus)  continue;
        	 if (hitSecond->time() -  hitFirst->time()  > deltaTimeMinus) continue;

		 double dist =  closestDistance(FirstHitsPtr,SecondHitsPtr);
		 if (dist > maxDist) continue;

		 if (dist < minDist) {minDist=dist; jmin=j;} 
	      }	  
              
	      _associatedMainId[i]   = jmin;
	  }
	    	    
       } 
 
 
       double CaloClusterAssociator::closestDistance(CaloCrystalHitPtrVector const& cluster, CaloCrystalHitPtrVector const& cluster2)
       {
	  double minDistance(1e6);    
	  for (auto const& hit : cluster)
	  {        
	      CLHEP::Hep3Vector crystalPos = _cal.crystalOrigin(hit->id());

	      for (auto const& hit2 : cluster2)
	      {
		 CLHEP::Hep3Vector crystalPos2 = _cal.crystalOrigin(hit2->id());
		 double dist = (crystalPos-crystalPos2).mag();
		 if (dist<minDistance) minDistance = dist;
	      }	 
	  }

	 return minDistance;
       }


}
