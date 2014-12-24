//
// Utility to study the MC content of a calo cluster
// 
// If the CaloHitSimPartMC information is not available, recompute it
// 
// Original author B. Echenard
//

// Mu2e includes
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "CaloCluster/inc/CaloSeedManager.hh"

// C++ includes
#include <iostream>
#include <queue>
#include <list>
#include <unordered_map>



namespace mu2e {


       void CaloSeedManager::add(CaloCrystalHit const& i) 
       { 
	    int id = i.id();

	    if (_seedMap[id]==0) _seedMap[id] = &i; 
	    else
	    {
	      if (_seedMode==SeedType::Energy &&  (i.energyDep() > _seedMap[id]->energyDep() )) _seedMap[id] = &i;
	      if (_seedMode==SeedType::Time   &&  (i.time()      < _seedMap[id]->time()      )) _seedMap[id] = &i;
	    }  
       }   


       CaloCrystalHit const* CaloSeedManager::seed()
       { 
	    if (_seedMap.empty()) return nullptr;

	    CaloCrystalHit const* seed(0);

	    for (unsigned int i=0;i<_seedMap.size();++i)
	    {
		if ( _seedMap[i] == 0 ) continue;
	        if (seed==0) {seed = _seedMap[i]; continue;} 
		if (_seedMode==SeedType::Energy && (_seedMap[i]->energyDep() > seed->energyDep() ) )  seed = _seedMap[i];
	        if (_seedMode==SeedType::Time   && (_seedMap[i]->time()      < seed->time()      ) )  seed = _seedMap[i];
	    }   
	    return seed;
       }
	 	 

       void CaloSeedManager::checkSeedbyList(CaloCrystalList const& crystalsInCluster, std::vector<CaloCrystalVec> const& idHitMap )
       {        

	   for (auto const& hit : crystalsInCluster)
	   {
		int iId = hit->id();
		if (hit == _seedMap[iId]) checkSeedbyId(iId,idHitMap[iId]); 
           }
       }

       
       void CaloSeedManager::checkSeedbyId(int iId, CaloCrystalVec const& hits)
       {
                
           if (hits.empty()) {_seedMap[iId]=0; return;}     

           CaloCrystalHit const* newSeed(0);
           for (unsigned int i=0; i< hits.size(); ++i)
           {
              if (hits[i]==0) continue;
	      if (newSeed==0) {newSeed = hits[i]; continue;} 
	      if ( _seedMode==SeedType::Energy && ( hits[i]->energyDep() > newSeed->energyDep() ) ) newSeed = hits[i];
              if ( _seedMode==SeedType::Time   && ( hits[i]->time()      < newSeed->time()      ) ) newSeed = hits[i];
           }
           
           _seedMap[iId] = newSeed;       
       }         

      
       	 
    
       void CaloSeedManager::dumpSeed()
       {         
	   std::cout<<"Seeds list"<<std::endl;
	   for (unsigned int i=0;i<_seedMap.size();++i) if (_seedMap[i] !=0) std::cout<<i<<" "<<_seedMap[i]->energyDep()<<std::endl;
           std::cout<<"End seeds "<<std::endl;
       }	 
}
