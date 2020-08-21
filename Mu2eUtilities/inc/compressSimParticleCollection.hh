#ifndef Mu2eUtilities_compressSimParticleCollection_hh
#define Mu2eUtilities_compressSimParticleCollection_hh

//
// Compress a SimParticleCollection.
//
//
// Contact person Rob Kutschke
//
// Notes:
// 1) This function copies an existing SimParticleCollection to a new one, but removes
//    SimParticles that are flagged as uninteresting.
//
// 2) This requires maintenance on the mother/daughter links to remove references to particles
//    that no longer exist.
//
// 3) It also requires reseating the mother/daughter art::Ptr objects to point into the new collection.
//
// 4) The art::Ptr objects within the output collection are not usable within the module that creates
//    the output collection; they are usable in all subsequent modules.
//
// 5) The arguments are:
//    1 - the art ProductID of the output collection
//    2 - a productGetter object to get items from the output collection
//    3 - the input collection
//    4 - the object that knows whether to keep or delete each item - see note 7.
//    5 - the output collection.
//
// 6) The code will throw if you try to save a secondary particle without also saving its mother.
//
// 7) The SELECTOR template argument can be any class that supports the following:
//       bool operator[]( cet::map_vector_key) const;
//
//    One example of a type that will work is:
//       cet::map_vector<bool>
//    For other examples see Filters/src/HitsInConversionTimeWindow_module.cc
//
//    Note that std::map<cet::map_vector_key,bool> will not work because it does not have
//    a const operator[].
//
// 8) A primary particle does not have a mother and this is represented by a mother Ptr
//    with a null key; such a Ptr needs no maintenance.  A secondary particle must have
//    a mother.  This code with throw if one tries to create an output collection in which
//    the mother of a secondary has been deleted.
//

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Common/EDProductGetter.h"

namespace mu2e {

  typedef std::map<cet::map_vector_key, cet::map_vector_key> KeyRemap;

  // Pass in the old key to check if it's already added to keyRemap, if it hasn't been then use nextNewKey for the next key
  cet::map_vector_key getNewKey(const cet::map_vector_key& oldKey, KeyRemap* keyRemap, const unsigned int& nextNewKey) {
    cet::map_vector_key nextKey;

    if ( keyRemap->find(oldKey) == keyRemap->end() ) { // might have already added the key since parents have a position reserved before they are added to the output
      nextKey = cet::map_vector_key(nextNewKey);
      keyRemap->insert( std::make_pair(oldKey, nextKey) ); // update the map
    }
    else {
      nextKey = keyRemap->at(oldKey);
    }

    return nextKey;
  }


  template<typename SELECTOR, typename OUTCOLL>
  void compressSimParticleCollection ( art::ProductID         const& newProductID,
                                       art::EDProductGetter   const* productGetter,
                                       SimParticleCollection  const& in,
                                       SELECTOR               const& keep,
                                       OUTCOLL&        out,
				       KeyRemap* keyRemap = NULL){

    unsigned int initial_out_size = out.size();
    for ( SimParticleCollection::const_iterator i=in.begin(), e=in.end(); i!=e; ++i ){
      if ( keep[i->first] ){

        // Default construct and replace to avoid multiple searches through the collection.
	cet::map_vector_key oldSimKey = i->first;
	cet::map_vector_key newSimKey;
	if (keyRemap) {
	  newSimKey = getNewKey(oldSimKey, keyRemap, initial_out_size + keyRemap->size());
	}
	else { 
	  newSimKey = oldSimKey;
	}
        SimParticle& sim = out[newSimKey];
        sim = i->second;

	if (keyRemap) {
	  sim.id() = newSimKey; // need to make sure the SimParticle's trackId is the same as its key in the output collection
	}

        // See note 1).
        if ( sim.isSecondary() ){
          cet::map_vector_key parentKey = cet::map_vector_key(sim.parent().key());
          if ( keep[parentKey] ) {
	    art::Ptr<SimParticle> newParentPtr;
	    if (keyRemap) {
	      cet::map_vector_key newParentKey = getNewKey(parentKey, keyRemap, initial_out_size + keyRemap->size());
	      newParentPtr = art::Ptr<SimParticle>( newProductID, newParentKey.asUint(), productGetter);
	    }
	    else {
	      newParentPtr = art::Ptr<SimParticle>( newProductID, sim.parent().key(), productGetter);
	    }
            sim.parent() = newParentPtr;
          } else{
            // This particle is a secondary particle but does not have a mother.
            // Do we wish to throw or warn here?
          }
        }

        // Remove daughters that have been deleted.
        std::vector<art::Ptr<SimParticle> > daughters;
        std::vector<art::Ptr<SimParticle> > const& oldDaughters(sim.daughters());
        for ( std::vector<art::Ptr<SimParticle> >::const_iterator j=oldDaughters.begin(),
                je=oldDaughters.end(); j !=je; ++j ){
          cet::map_vector_key dkey = cet::map_vector_key(j->key());

          if ( keep[dkey] ){
	    art::Ptr<SimParticle> newDPtr;
	    if (keyRemap) {
	      cet::map_vector_key newDKey = getNewKey(dkey, keyRemap, initial_out_size + keyRemap->size());
	      newDPtr = art::Ptr<SimParticle>( newProductID, newDKey.asUint(), productGetter);
	    }
	    else {
	      newDPtr = art::Ptr<SimParticle>( newProductID, j->key(), productGetter);
	    }
            daughters.push_back(newDPtr);
          }
        }
        sim.setDaughterPtrs(daughters);
      }
    }
  } // end compressSimParticleCollection

}
#endif /* Mu2eUtilities_compressSimParticleCollection_hh */
