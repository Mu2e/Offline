//
// Class to find cluster of simply connected crystals
// 
// Original author B. Echenard
//
// Note: there are few places where a continue could be replaced by a break if the crystal are time ordered, but 
//       the performance gain is so low that it outweigh the risk of forgeting to time order the crystal hits.
//       

// Mu2e includes
#include "CaloCluster/inc/CaloClusterFinder.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

// C++ includes
#include <iostream>
#include <vector>
#include <algorithm>


namespace mu2e {


       CaloClusterFinder::CaloClusterFinder(Calorimeter const& cal, CaloCrystalHit const* crystalSeed, double deltaTimePlus, double deltaTimeMinus, double ExpandCut) : 
          _cal(&cal), _crystalSeed(crystalSeed),_seedTime(crystalSeed->time()), _clusterList(), _inspected(), _crystalToVisit(), _isVisited(cal.nCrystal()),
	  _deltaTimePlus(deltaTimePlus), _deltaTimeMinus(deltaTimeMinus), _ExpandCut(ExpandCut)
       {}
       
       


       void CaloClusterFinder::formCluster(std::vector<CaloCrystalVec>& idHitVec)  
       { 

	    _inspected.clear();
	    
	    _clusterList.push_front(_crystalSeed);
	    _crystalToVisit.push(_crystalSeed->id());  
            _inspected.push_front(_crystalSeed);


	    
	    CaloCrystalVec& vecSeed = idHitVec[_crystalSeed->id()];
	    for (unsigned int ivec=0; ivec < vecSeed.size(); ++ivec)
	       if (vecSeed[ivec]==_crystalSeed) {vecSeed[ivec]=0; break;}


	    while (!_crystalToVisit.empty())
	    {            
        	 int visitId         = _crystalToVisit.front();
		 _isVisited[visitId] = 1;

		 std::vector<int> const& neighborsId = _cal->neighbors(visitId);
         	 for (auto& iId : neighborsId)
        	 {               
		     if (_isVisited[iId]) continue;
		     _isVisited[iId]=1;

		     CaloCrystalVec& vec = idHitVec[iId];
		     for (unsigned int ivec=0; ivec < vec.size(); ++ivec)
		     {
			 if (vec[ivec]==0) continue;
			 CaloCrystalHit const* hit = vec[ivec];
			 

			 if (hit->time() - _seedTime > _deltaTimePlus)  continue; //see note
			 if (_seedTime - hit->time() > _deltaTimeMinus) continue;

         		 if (hit->energyDep() > _ExpandCut) _crystalToVisit.push(iId);
                         _clusterList.push_front(hit);

 			 _inspected.push_front(hit);			 
        		 vec[ivec]=0;
		     } 
		     
		 }
		      
		 _crystalToVisit.pop();		 
	    }
	    
	   // keep sorting proto-clustres, even if they are sorted in the cluster module, in case somebody 
	   // uses the proto-clusters instead of clusters he will get the same behaviour
	   _clusterList.sort([] (CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();} );
	    	    
       } 

}



/*
       void CaloClusterFinder::cutterCluster(std::vector<CaloCrystalVec>& idHitVec)
       {  
	    _inspected.clear();
	    
	    _clusterList.push_front(_crystalSeed);
	    _crystalToVisit.push(_crystalSeed->id());  
            _inspected.push_front(_crystalSeed);

            std::vector<int> listToVisit = _cal->neighbors(_crystalSeed->id());
            listToVisit.insert( listToVisit.end(), _cal->nextNeighbors(_crystalSeed->id()).begin(), _cal->nextNeighbors(_crystalSeed->id()).end() );
 
 	    CaloCrystalVec& vecSeed = idHitVec[_crystalSeed->id()];
	    for (unsigned int ivec=0; ivec < vecSeed.size(); ++ivec)
	       if (vecSeed[ivec]==_crystalSeed) {vecSeed[ivec]=0; break;}
           
	
	    for (int iId : listToVisit)
            {               
		CaloCrystalVec& vec = idHitVec[iId];
		for (unsigned int ivec=0; ivec < vec.size(); ++ivec)
		{
		    if (vec[ivec]==0) continue;
		    CaloCrystalHit const* hit = vec[ivec];

		    if (hit->time() - _seedTime > _deltaTimePlus)  continue;
		    if (_seedTime - hit->time() > _deltaTimeMinus) continue;

                    _clusterList.push_front(hit);
 		    _inspected.push_front(hit);
                    vec[ivec]=0;			 
		} 
		     
            }
	   _clusterList.sort([] (CaloCrystalHit const* lhs, CaloCrystalHit const* rhs) {return lhs->energyDep() > rhs->energyDep();} );	    	    
       } 
*/









