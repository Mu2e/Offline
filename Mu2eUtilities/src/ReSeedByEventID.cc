//
// Reseed a random engine using the art::EventID and a user supplied salt.
//
// Notes:
// 1) The current algorithm requires an engine that can accept 4 seeds.
//    We could adapt the algorithm to fewer seeds but have chosen not to
//    unless and until it is needed.  We know that the MixMaxRng engine
//    does support up to 4 seeds and, for the foreseeable future, it is
//    the only engine we plan to use.
//
// 2) We have not checked all CLHEP engines but we know that the following
//    do NOT support 4 seeds:
//       - HepJamesRandom
//       - Ranlux64Engine
//       - MTwistEngine
//
// 3) Because of 1 and 2 the c'tor whitelists only MixMaxRng.
//
// Fixme: The code may require maintenance if we change to a different engine type.
//

#include "Offline/Mu2eUtilities/inc/ReSeedByEventID.hh"
#include "fhiclcpp/exception.h"

#include <iostream>
#include <string>

namespace mu2e {
  ReSeedByEventID::ReSeedByEventID( CLHEP::HepRandomEngine& engine, SeedService::seed_t salt, int verbosity)
    :engine_(engine)
    ,salt_{salt}
    ,verbosity_{verbosity}
  {
    const std::string engineType{engine.name()};

    if ( verbosity_ > 0 ){
      std::cout << "ReSeedbyEventID.  Engine name:             " << engineType  << std::endl;
      std::cout << "                  salt:                    " << salt_       << std::endl;
    }

    // Only accept whitelisted engines.  See notes above.
    if ( engineType != "MixMaxRng" ){
      throw cet::exception("BADCONFIG")
        << "ReSeedByEventID requires that engine::name() be MixMaxRng.  It is: " << engineType << "\n";
    }

  }

  void  ReSeedByEventID::reseed ( art::EventID const& id ) const{

    if ( verbosity_ > 1 ){
      std::cout << "ReseedByEventID: " << id  << std::endl;
    }

    const std::array<long,4> seeds{ id.run(), id.subRun(), id.event(), salt_};
    engine_.setSeeds( seeds.data(), seeds.size() );

  }

} // end namespace mu2e
