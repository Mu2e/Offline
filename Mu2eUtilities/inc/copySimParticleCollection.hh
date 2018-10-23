#ifndef Mu2eUtilities_copySimParticleCollection_hh
#define Mu2eUtilities_copySimParticleCollection_hh

//
// Copy a SimParticleCollection and reseat Ptr's
//
// When making a new SimparticleCollection from an existing collection,
// it is not as simple to copy the objects.
// The art::Ptr's in the SimParticles need to be fixed so they point 
// the their parents (and daughters) in the new collection.
// 
// The inputs are the old collection, the new collection,
// and the art product pieces needed for making art::Ptr's
// to the new collection.  Examples of how to get those 
// pieces are in Filters.
// 
// 

#include <memory>
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Common/EDProductGetter.h"

namespace mu2e {

  void copySimParticleCollection (SimParticleCollection  const& oldColl,
				SimParticleCollection& newColl,
				art::ProductID const& SPpid,
				art::EDProductGetter const* SPpg) {

    // loop over the collection
    for(auto const& spPair : oldColl) {
      auto key = spPair.first;
      auto const& oldSP = spPair.second;
      // copy the SimParticle
      newColl[key] = oldSP;
      auto& newSP = newColl[key];
      // point the parent art::Ptr to the parent in the new collection
      newSP.parent() = 
	art::Ptr<SimParticle>(SPpid,oldSP.parent().key(),SPpg);
      // and repoint daughter art::Ptr's
      for (auto& dauPtr : newSP.daughters()) {
	dauPtr = art::Ptr<SimParticle>(SPpid,dauPtr.key(),SPpg);
      }
    } // loop over SimParticles

    return;
  } // end copySimParticleCollection

}


#endif
