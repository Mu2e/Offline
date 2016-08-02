//
// Utility to study the MC content of a calo cluster. Browser over SimParticles of each crystal, and keep only distinct
// entries, updating the total energy, time and position
//
//
// Original author B. Echenard
//

#include "CaloMC/inc/ClusterContentMC.hh"
#include "CaloMC/inc/CaloContentSim.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <map>
#include<iostream>



namespace mu2e {


       ClusterContentMC::ClusterContentMC(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster) :
          _simContentMap(), _hasConversion(false), _eDepTot(0), _time(0)
       {
           fillCluster(cal, caloClusterTruth, cluster);
       }


       void ClusterContentMC::fillCluster(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster)
       {           
	    for (auto i=caloClusterTruth.begin(), ie = caloClusterTruth.end(); i !=ie; ++i)
	    {	       
	        const auto& caloClusterPtr = i->first;
		const auto& sim = i->second;
		const auto& caloShowerPtr = caloClusterTruth.data(i);	     

		if (caloClusterPtr.get() != &cluster) continue;
		
		double eDep = caloShowerPtr->energy();
		double time = caloShowerPtr->time();

		_eDepTot += eDep;
		if ( std::abs(time-caloClusterPtr->time()) < std::abs(_time-caloClusterPtr->time()) ) _time = time; 

		double pIn(0);
		for (const auto& steps : caloShowerPtr->caloShowerSteps()) pIn = std::max(pIn,steps->momentumIn());
		auto parent(sim);
                while ( parent->hasParent() && cal.isInsideCalorimeter(parent->startPosition()) ) parent = parent->parent();	             
                
		auto mfind = _simContentMap.find(parent);		     
		if (mfind != _simContentMap.end()) mfind->second.update(eDep,time,pIn);
		else _simContentMap.insert(std::pair<art::Ptr<SimParticle>, CaloContentSim>(parent,CaloContentSim(eDep,time,pIn)) );
		
		if (parent->genParticle() && parent->genParticle()->generatorId()==GenId::conversionGun) _hasConversion = true; 		
	    }

       }

}
