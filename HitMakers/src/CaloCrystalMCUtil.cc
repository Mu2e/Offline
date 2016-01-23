//
// Utility to study the MC content of a calorimeter crystal hit
//
// Definitions:
//   - "SimParticle ancestor" the SimParticle that crossed the calorimeter boundary and produced the StepPointMC hit. Note 
//     that the SimParticle ancestor might produce other SimParticles which will then produce StepPointMCs
//   - "SimParticle" the SimParticle that directly produced the StepPointMC
//   - "SimParticle mother" (aka mommy) the SimParticle that produced this SimParticle. 
//
//  Find the total energy deposited by each SimParticle ancestor in a crystal, and the StepPointMc with the largest momentum 
//  The time/position of the StepPointMC with the largest momentum = origin of the shower in the crystal (in good approximation), 
//  since this information is not stored by the SimParticles.
//
//
//  N.B: We could have used a virtual detector to keep track of the energy/position of particles hitting the calorimeter, 
//       but we lose the information about the shower development in the crystals. Could do so if we further need to 
//       shrink the data size
// 
//       Incoming SimParticles produces usually a few additional SimParticles along the way. For this reason, it is 
//       more efficient to look at the SimParticle and the SimParticle mother associated to a StepPointMC when browsing the 
//       StepPointMC list instead of just looking at the SimParticle
//
// Original author B. Echenard
//


// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <unordered_map>

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "HitMakers/inc/CaloCrystalMCUtil.hh"



//-----------------------------------------

namespace mu2e {


    
     void CaloCrystalMCUtil::fillSimMother(Calorimeter const& cal, PtrStepPointMCVector const& mcptr, CaloHitSimPartMC& caloHitSimPartMC, 
                                           const std::unordered_map<const StepPointMC*, double> &timeMap)
     {
       
	 SimStepMap simStep;
	 std::map<int, SimPtr > simAncestor;

	 for (unsigned int i=0; i<mcptr.size(); ++i) 
	 {

	     StepPtr const& mchit     = mcptr[i];
	     int thisSimpartId        = mchit->simParticle()->id().asInt();
	     int thisParentSimpartId  = (mchit->simParticle()->parent()) ? mchit->simParticle()->parent()->id().asInt() : -1;

	     //SimParticle ancestor already stored, update CaloHitMCUtil object
	     std::map<int, SimPtr >::iterator ifind = simAncestor.find(thisSimpartId); 
	     if ( ifind != simAncestor.end() )
	     {
		 SimPtr const& ancestor     = ifind->second;
		 SimStepMap::iterator afind = simStep.find(ancestor);		 
		 if (afind != simStep.end() ) afind->second.update(mchit,mchit->totalEDep());
		 continue;		 
	     }	    

             //SimParticle parent is already in the ancestor list, update ancestor list and CaloHitMCUtil object
	     std::map<int, SimPtr>::iterator ifindParent = simAncestor.find(thisParentSimpartId);              
	     if ( ifindParent != simAncestor.end() )
	     {
		 SimPtr const& ancestor = ifindParent->second;
		 simAncestor.insert( std::pair<int,SimPtr >(thisSimpartId, ancestor) );

		 SimStepMap::iterator afind = simStep.find(ancestor);		 
		 if ( afind != simStep.end() ) afind->second.update(mchit,mchit->totalEDep());
		 continue;
	     }	    

	     //find SimParticle ancestor, and store the corresponding SimParticle Id -> SimParticle ancestor
	     SimPtr parent = mchit->simParticle();
             while ( parent->hasParent() && cal.isInsideCalorimeter(parent->startPosition()) ) parent = parent->parent();
	     simAncestor.insert( std::pair<int, SimPtr >(thisSimpartId, parent) );

	     //store or update the SimParticle ancestor and the associated CaloHitMCUtil object 
	     SimStepMap::iterator afind = simStep.find(parent);		 
	     if (afind != simStep.end()) afind->second.update(mchit,mchit->totalEDep()); 
	     else simStep.insert(std::pair<SimPtr, CaloHitMCUtil>(parent,CaloHitMCUtil(mchit,mchit->totalEDep())) );
	  }


	  for (SimStepMap::const_iterator it=simStep.begin();it!=simStep.end();++it)
	  {        
	     StepPtr const& mchit = it->second.step();
	     
	     
	     double htime = mchit->time();
	     auto iter = timeMap.find(&(*mchit));
	     if (iter != timeMap.end()) htime = iter->second;
	     
	     caloHitSimPartMC.add(it->first,it->second.edepTot(),htime,mchit->momentum().mag(),mchit->position());
	  }        	

     }


}
