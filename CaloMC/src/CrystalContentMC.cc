//
// Utility to study the MC content of a calo cluster. Browser over SimParticles of each crystal, and keep only distinct
// entries, updating the total energy, time and position
//
// Now there is a trick here. The truch matching is between SimParticle and CaloHits with a data payload about the CaloShower.
// This means that there can be several caloShowers for a given SimParticle / CaloHit entry. BUT this also means there can be 
// several caloShowers for a given CaloHit if several SimParticle contibutes to the Shower. In that case, there is a double-counting 
// problem, prevented by the  caloShowerSeen flag
//  

#include "CaloMC/inc/CrystalContentMC.hh"
#include "CaloMC/inc/CaloContentSim.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

#include <map>
#include<iostream>



namespace mu2e {


       CrystalContentMC::CrystalContentMC(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit) :
          simContentMap_(), eDepTot_(0), time_(0)
       {
           fillCrystal(cal, caloHitTruth, caloCrystalHit);
       }


       void CrystalContentMC::fillCrystal(const Calorimeter& cal, const CaloHitMCTruthAssns& caloHitTruth, const CaloCrystalHit& caloCrystalHit)
       {           
           std::set<const CaloShowerSim*> caloShowerSimSeen;

           for (auto i=caloHitTruth.begin(), ie = caloHitTruth.end(); i !=ie; ++i)
           {               
               const auto& caloCrystalHitPtr = i->first;
               const auto& sim = i->second;
               const auto& caloShowerSimPtr = caloHitTruth.data(i);             

               if (caloCrystalHitPtr.get() != &caloCrystalHit) continue;                
               if (caloShowerSimSeen.find(caloShowerSimPtr.get()) != caloShowerSimSeen.end()) continue;
               caloShowerSimSeen.insert(caloShowerSimPtr.get());                

               double pIn = caloShowerSimPtr->momentumIn();
               double eDep = caloShowerSimPtr->energy();
               double time = caloShowerSimPtr->time();

               eDepTot_ += eDep;
               if ( std::abs(time-caloCrystalHitPtr->time()) < std::abs(time_-caloCrystalHitPtr->time()) ) time_ = time; 

               auto parent(sim);
               while ( parent->hasParent() && cal.geomUtil().isInsideCalorimeter(parent->startPosition()) ) parent = parent->parent();                     

               auto mfind = simContentMap_.find(parent);                     
               if (mfind != simContentMap_.end()) mfind->second.update(eDep,time,pIn);
               else simContentMap_.insert(std::pair<art::Ptr<SimParticle>, CaloContentSim>(parent,CaloContentSim(eDep,time,pIn)) );

	       if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) hasConversion_ = true; 		
          }
      }


}
