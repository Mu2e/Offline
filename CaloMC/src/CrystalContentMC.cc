//
// Utility to study the MC content of a calo cluster. Browser over SimParticles of each crystal, and keep only distinct
// entries, updating the total energy, time and position
//
//
// Original author B. Echenard
//


#include "CaloMC/inc/CrystalContentMC.hh"
#include "CaloMC/inc/CaloContentSim.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <map>
#include<iostream>



namespace mu2e {


       CrystalContentMC::CrystalContentMC(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit) :
          _simContentMap(), _eDepTot(0), _time(0)
       {
           fillCrystal(cal, caloHitTruth, caloCrystalHit);
       }


       void CrystalContentMC::fillCrystal(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit)
       {           
	    
	    for (auto i=caloHitTruth.begin(), ie = caloHitTruth.end(); i !=ie; ++i)
	    {	       
	        const auto& caloCrystalHitPtr = i->first;
		const auto& sim = i->second;
		const auto& caloShowerPtr = caloHitTruth.data(i);	     

		if (caloCrystalHitPtr.get() != &caloCrystalHit) continue;
		
		double eDep = caloShowerPtr->energy();
		double time = caloShowerPtr->time();

		_eDepTot += eDep;
		if ( std::abs(time-caloCrystalHitPtr->time()) < std::abs(_time-caloCrystalHitPtr->time()) ) _time = time; 
		
		double pIn(0);
		for (const auto& steps : caloShowerPtr->caloShowerSteps()) pIn = std::max(pIn,steps->momentumIn());
		
		auto parent(sim);
                while ( parent->hasParent() && cal.isInsideCalorimeter(parent->startPosition()) ) parent = parent->parent();	             
                
		auto mfind = _simContentMap.find(parent);		     
		if (mfind != _simContentMap.end()) mfind->second.update(eDep,time,pIn);
		else _simContentMap.insert(std::pair<art::Ptr<SimParticle>, CaloContentSim>(parent,CaloContentSim(eDep,time,pIn)) );
	   }

       }

}
