//
// Utility to study the MC content of a calorimeter crystal hit
// Find the SimParticle(s) associated to each hits at the calorimeter boundary
//  - for each StepPointMC, find mother by going up the tree of SimParticles until start or end is outside the calo (jumps between calo sections are ok) 
//  - record the mchit_ptr for the mother particles - momentum/direction is not recorded in SimParticles at the calorimeter boundary
//    so must use StePointMC to estimate the momentum at the entrance of the calo 
// Original author B. Echenard
//

// C++ includes
#include <map>
#include <vector>
#include <iostream>

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "HitMakers/inc/CaloReadoutUtilities.hh"
#include "HitMakers/inc/CaloHitSimUtil.hh"



//-----------------------------------------

namespace mu2e {


    
     void CaloReadoutUtilities::fillSimMother(Calorimeter const& cal, 
                        		      PtrStepPointMCVector const& mcptr,
					      CaloHitSimPartMC& caloHitSimPartMC)
     {
       
       SimStepMap simStep;
       std::map<int, SimPtr > simParents;

       for (unsigned int i=0; i<mcptr.size(); ++i) 
       {
	   
	   StepPtr const& mchit = mcptr[i];
	   int thisSimpartId           = mchit->simParticle()->id().asInt();
	   int parentSimpartId         = (mchit->simParticle()->parent()) ? mchit->simParticle()->parent()->id().asInt() : -1;

	   //SimParticle already in the list
	   std::map<int, SimPtr >::iterator ifind = simParents.find(thisSimpartId); 
	   if ( ifind != simParents.end() )
	   {
	       SimPtr const& mommy = ifind->second;
	       SimStepMap::iterator mfind = simStep.find(mommy);		 
	       if (mfind != simStep.end() ) mfind->second.update(mchit,mchit->totalEDep());
	       continue;		 
	   }	    

           //SimParticle mother already in the list
	   std::map<int, SimPtr>::iterator ifindParent = simParents.find(parentSimpartId);              
	   if ( ifindParent != simParents.end() )
	   {
	       SimPtr const& mommy = ifindParent->second;
	       simParents.insert( std::pair<int,SimPtr >(thisSimpartId, mommy) );

	       SimStepMap::iterator mfind = simStep.find(mommy);		 
	       if ( mfind!=simStep.end() ) mfind->second.update(mchit,mchit->totalEDep());
	       continue;
	   }	    

	   //SimParticle is added to the list of simParticles
	   SimPtr mother = mchit->simParticle();
           while ( mother->hasParent() && cal.isInsideCalorimeter(mother->startPosition()) ) mother = mother->parent();

	   simParents.insert( std::pair<int, SimPtr >(thisSimpartId, mother) );

   	   SimStepMap::iterator mfind = simStep.find(mother);		 
	   if (mfind != simStep.end()) mfind->second.update(mchit,mchit->totalEDep()); 
	   else simStep.insert(std::pair<SimPtr, CaloHitSimUtil>(mother,CaloHitSimUtil(mchit,mchit->totalEDep())) );

	}
		        
	for (SimStepMap::const_iterator it=simStep.begin();it!=simStep.end();++it)
	{        
	   caloHitSimPartMC.add(it->first,it->second.step(),it->second.edepTot());
	}        	

     }


}
