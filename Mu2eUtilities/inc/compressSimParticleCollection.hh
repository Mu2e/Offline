#ifndef Mu2eUtilities_compressSimParticleCollection_hh
#define Mu2eUtilities_compressSimParticleCollection_hh

//
// Compress a SimParticleCollection.
//
// $Id: compressSimParticleCollection.hh,v 1.2 2013/08/28 05:58:37 gandr Exp $
// $Author: gandr $
// $Date: 2013/08/28 05:58:37 $
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

  // Pass in the old key to check if it's already added to keyRemap, if it ahsn't been then use nextNewKey for the next key
  cet::map_vector_key getNewKey(const cet::map_vector_key& oldKey, KeyRemap& keyRemap, const unsigned int& nextNewKey) {
    cet::map_vector_key nextKey;
    //    std::cout << "oldKey = " << oldKey << " ";
    if ( keyRemap.find(oldKey) == keyRemap.end() ) { // might have already added the key since we add parents earlier when remapping
      nextKey = cet::map_vector_key(nextNewKey);
      keyRemap.insert( std::make_pair(oldKey, nextKey) ); // update the map
      //      std::cout << " not added yet, newKey = " << nextKey << ")" << std::endl;
    }
    else {
      nextKey = keyRemap.at(oldKey);
      //      std::cout << " already added, newKey = " << nextKey << ")" << std::endl;
    }

    return nextKey;
  }


  template<typename SELECTOR, typename OUTCOLL>
  void compressSimParticleCollection ( art::ProductID         const& newProductID,
                                       art::EDProductGetter   const* productGetter,
                                       SimParticleCollection  const& in,
                                       SELECTOR               const& keep,
                                       OUTCOLL&        out,
				       SimParticleRemapping* remap = NULL){

    unsigned int initial_out_size = out.size();
    KeyRemap keyRemap;
    for ( SimParticleCollection::const_iterator i=in.begin(), e=in.end(); i!=e; ++i ){
      if ( keep[i->first] ){

        // Default construct and replace to avoid multiple searches through the collection.
	bool no_kept_parents_or_daughters = true;
	cet::map_vector_key oldSimKey = i->first;
	//	std::cout << "i->first = " << i->first << std::endl;
	cet::map_vector_key newSimKey;
	if (remap) {
	  newSimKey = getNewKey(oldSimKey, keyRemap, initial_out_size + keyRemap.size());
	  //	  std::cout << "Get new i->first key = " << newSimKey << std::endl;
	}
	else { 
	  newSimKey = oldSimKey;
	}
        SimParticle& sim = out[newSimKey];
        sim = i->second;

        // See note 1).
        if ( sim.isSecondary() ){
          cet::map_vector_key parentKey = cet::map_vector_key(sim.parent().key());
	  //	  std::cout << "Should have a parent" << std::endl;
          if ( keep[parentKey] ) {
	    no_kept_parents_or_daughters = false;
	    art::Ptr<SimParticle> newParentPtr;
	    if (remap) {
	      cet::map_vector_key newParentKey = getNewKey(parentKey, keyRemap, initial_out_size + keyRemap.size());
	      //	      std::cout << "Get new parent key = " << newParentKey << std::endl;
	      newParentPtr = art::Ptr<SimParticle>( newProductID, newParentKey.asUint(), productGetter);
	      (*remap)[sim.parent()] = newParentPtr;
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
	  //	  std::cout << "Should have daughters" << std::endl;
          if ( keep[dkey] ){
	    no_kept_parents_or_daughters = false;
	    art::Ptr<SimParticle> newDPtr;
	    if (remap) {
	      cet::map_vector_key newDKey = getNewKey(dkey, keyRemap, initial_out_size + keyRemap.size());
	      //	      std::cout << "Get new daughter key = " << newDKey << std::endl;
	      newDPtr = art::Ptr<SimParticle>( newProductID, newDKey.asUint(), productGetter);
	      (*remap)[*j] = newDPtr;
	    }
	    else {
	      newDPtr = art::Ptr<SimParticle>( newProductID, j->key(), productGetter);
	    }
            daughters.push_back(newDPtr);
          }
        }
        sim.setDaughterPtrs(daughters);

	if (remap) {
	  // Sometimes a SimParticle has no parents or daughters that are being kept
	  // In this case we haven't seen the old SimParticle Ptr to be able to fill the remap
	  // So we have to create the oldSimPtr ourselves and infer the old ProductID and EDProductGetter from other SimParticles in the input collection
	  if (no_kept_parents_or_daughters) {
	    // need to work out the SimParticle Ptr ourselves
	    art::Ptr<SimParticle> aDaughter = art::Ptr<SimParticle>();
	    for (const auto& aSimPartPair : in) {
	      const auto& aSimParticle = aSimPartPair.second;
	      if (aSimParticle.daughters().size()>0) {
		aDaughter = *(aSimParticle.daughters().begin());
		
		art::ProductID oldPID = aDaughter.id();
		const art::EDProductGetter* oldProductGetter = aDaughter.productGetter();
	    
		art::Ptr<SimParticle> oldSimPtr = art::Ptr<SimParticle>(oldPID, oldSimKey.asUint(), oldProductGetter);
		art::Ptr<SimParticle> newSimPtr = art::Ptr<SimParticle>(newProductID, newSimKey.asUint(), productGetter);
		(*remap)[oldSimPtr] = newSimPtr;
		  
		break;
	      }
	    }
	  }
	}
      }
    }

    /*    if (remap) {
      std::cout << "Final Remapping: " << std::endl;
      std::cout << "Output Size = " << out.size() << ", Max Size = " << out.max_size() << ", Remap size = " << remap->size() << std::endl;
      for (const auto& i_simPtrPair : *remap) {
	std::cout << i_simPtrPair.first << " --> " << i_simPtrPair.second << std::endl;
      }
    }
    */
  } // end compressSimParticleCollection

}
#endif /* Mu2eUtilities_compressSimParticleCollection_hh */
