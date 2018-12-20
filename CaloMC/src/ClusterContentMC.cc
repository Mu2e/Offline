//
// Utility to study the MC content of a calo cluster. Browser over SimParticles of each crystal, and keep only distinct
// entries, updating the total energy, time and position
//

#include "CaloMC/inc/ClusterContentMC.hh"
#include "CaloMC/inc/CaloContentSim.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include <map>
#include<iostream>



namespace mu2e {


       ClusterContentMC::ClusterContentMC(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster) :
          simContentMap_(), hasConversion_(false), eDepTot_(0), time_(0)
       {
           fillCluster(cal, caloClusterTruth, cluster);
       }


       void ClusterContentMC::fillCluster(const Calorimeter& cal, const CaloClusterMCTruthAssns& caloClusterTruth, const CaloCluster& cluster)
       {           
	   std::set<const CaloShowerSim*> caloShowerSimSeen;

	   for (auto i=caloClusterTruth.begin(), ie = caloClusterTruth.end(); i !=ie; ++i)
	   {	       
	       const auto& caloClusterPtr = i->first;
	       const auto& sim = i->second;
	       const auto& caloShowerSimPtr = caloClusterTruth.data(i);	     

	       if (caloClusterPtr.get() != &cluster) continue;
	       if (caloShowerSimSeen.find(caloShowerSimPtr.get()) != caloShowerSimSeen.end()) continue;
	       caloShowerSimSeen.insert(caloShowerSimPtr.get());		

               double pIn  = caloShowerSimPtr->momentumIn();
               double eDep = caloShowerSimPtr->energy();
	       double time = caloShowerSimPtr->time();

	       eDepTot_ += eDep;
	       if ( std::abs(time-caloClusterPtr->time()) < std::abs(time_-caloClusterPtr->time()) ) time_ = time; 

	       auto parent(sim);
               while ( parent->hasParent() && cal.geomUtil().isInsideCalorimeter(parent->startPosition()) ) parent = parent->parent();	             

	       auto mfind = simContentMap_.find(parent);		     
	       if (mfind != simContentMap_.end()) mfind->second.update(eDep,time,pIn);
	       else simContentMap_.insert(std::pair<art::Ptr<SimParticle>, CaloContentSim>(parent,CaloContentSim(eDep,time,pIn)) );

	       if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) hasConversion_ = true; 		
	   }
       }

}
