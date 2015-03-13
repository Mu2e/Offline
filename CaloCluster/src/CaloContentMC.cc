//
// Utility to study the MC content of a calo cluster
// 
// If the CaloHitSimPartMC information is not available, recompute it
// 
// Original author B. Echenard
//

// Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CaloCluster/inc/CaloContentMC.hh"

#include "HitMakers/inc/CaloReadoutUtilities.hh"
#include "HitMakers/inc/CaloHitSimUtil.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"


// C++ includes
#include <map>
#include<iostream>



namespace mu2e {


       CaloContentMC::CaloContentMC(CaloHitMCNavigator const& navi, CaloCluster const& cluster) : 
	_navi(&navi), _cluster(&cluster), _simStep(), _simBaseSet(), _simPart(), _simBase(), _steps(), _edep()
       {	   
           fillCluster();
       }




       void CaloContentMC::fillCluster() 
       {
       
            Calorimeter const & cal = *(GeomHandle<Calorimeter>());
	    CaloReadoutUtilities readoutUtil;
	    
	    for (unsigned int i=0;i<_cluster->size();++i) 
	    {	 

		//run only over the first readout, since the other are duplicating the StepPointMc in the crystal
		CaloCrystalHit const& hit = *(_cluster->caloCrystalHitsPtrVector().at(i));
	        CaloHit const& caloHit = *(hit.readouts().at(0));
                CaloHitSimPartMC const& hitSim = _navi->sim(caloHit);

		//if the CaloHitSimPartMC is empty, then recalculate it on the fly
		if (hitSim.simParticles().size()) 
		{
		    fillMaps(hitSim);
		}
		else 
		{
	            CaloHitSimPartMC hitSimNew;
                    PtrStepPointMCVector const& mcptr = _navi->ptrs(caloHit);
	            readoutUtil.fillSimMother(cal,mcptr,hitSimNew);
		    fillMaps(hitSimNew);
		}
	           
	    }



            for (SimStepMap::const_iterator it=_simStep.begin();it!=_simStep.end();++it)
            {
	       _simPart.push_back(it->first);
	       _steps.push_back(it->second.step());
	       _edep.push_back(it->second.edepTot());
	    }         

	    for (std::set<SimParticlePtr >::const_iterator it=_simBaseSet.begin();it!=_simBaseSet.end();++it)
	       _simBase.push_back(*it);


       }


       void CaloContentMC::fillMaps(CaloHitSimPartMC const& hitSim)
       {

	    std::vector<SimParticlePtr > const& sim  = hitSim.simParticles();
	    std::vector<StepPointMCPtr > const& step = hitSim.stepPoints();
            std::vector<double>          const& edep = hitSim.eDep();

	    for (unsigned int i=0;i<sim.size();++i)
	    {	      
		
        	SimStepMap::iterator mfind = _simStep.find(sim[i]);            
        	if (mfind != _simStep.end()) mfind->second.update(step[i],edep[i]); 
        	else _simStep.insert(std::pair<SimParticlePtr, CaloHitSimUtil>(sim[i],CaloHitSimUtil(step[i],edep[i])) );

		SimParticlePtr mother = sim[i];
        	while ( mother->hasParent() ) mother = mother->parent();
  		_simBaseSet.insert(mother);                
            }
       }


       bool CaloContentMC::hasConversion() 
       {
	     for (std::set<SimParticlePtr >::const_iterator it=_simBaseSet.begin();it!=_simBaseSet.end();++it)
	     {	   
		 SimParticlePtr const& simpartPtr = *it;
		 if (simpartPtr->genParticle() && simpartPtr->genParticle()->generatorId()==GenId::conversionGun) return true;
	     }
	     return false;
       }	   


}
