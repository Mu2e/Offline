// Cache values of masses by PDG Id to reduce expensive
// lookups in the particle data table.
//
#include "Offline/GlobalConstantsService/inc/MassCache.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {
  MassCache::MassCache ():
    cache_(),
    lastMass_(0.),
    lastId_(PDGCode::unknown){
  }

  double MassCache::mass( id_type id ){

    // The easiest case.
    if ( id == lastId_ ) return lastMass_;
    lastId_  = id;

    // Next easiest case.
    // See if this particle id is already in the cache.
    //   - if so, get an iterator to it.
    //   - if not, construct a new entry with a placeholder mass.
    const double placeholder(0.);
    std::pair< map_type::iterator, bool> result = cache_.insert( value_type( id, placeholder) );

    // The particle id already existed in the map; return its mass.
    if ( !result.second ){
      lastMass_ = result.first->second;
      return lastMass_;
    }

    // The particle id was not in the map, so figure out its mass.
    lastMass_ = getMassFromPDT(id);

    // Insert the mass into the cache.
    result.first->second = lastMass_;

    return lastMass_;

  }

  // Get the mass from the particle data table
  double MassCache::getMassFromPDT( id_type id ){

    auto pdt = *GlobalConstantsHandle<ParticleDataList>();

    return pdt.particle(id).mass();
  }

}
