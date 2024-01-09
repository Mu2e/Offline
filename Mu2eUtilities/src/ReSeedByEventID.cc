#include "Offline/Mu2eUtilities/inc/ReSeedByEventID.hh"
//#include "art/Framework/EventProcessor/Scheduler.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "fhiclcpp/exception.h"

#include <iostream>
#include <string>

namespace mu2e {
  ReSeedByEventID::ReSeedByEventID( CLHEP::HepRandomEngine& engine, SeedService::seed_t salt):_salt(salt){
    std::string engineType{art::ServiceHandle<art::RandomNumberGenerator>{}->defaultEngineKind()};
    std::cout << "Kind:    " << engineType << std::endl;
    if ( engineType != "MixMaxRng" ){
      throw cet::exception("BADCONFIG")
        << "ReSeedByEventID requires that RandomNumberGenerator::defaultEngineKind be MixMaxRng.  It is: " << engineType << "\n";
    }
    // Ideally only do this when the number of schedules is >1.
    // The next line is not working.
    //std::cout << "Threads: " << art::ServiceHandle<art::Scheduler>{}->num_threads() << std::endl;
    _engine = static_cast<CLHEP::MixMaxRng*>(&engine);
  }

  void  ReSeedByEventID::reseed ( art::EventID const& id ) const{
    std::cout << "Reseed: " << id  << " " << _salt << std::endl;
    std::array<long,4> seeds;
    seeds[0] = id.run();
    seeds[1] = id.subRun();
    seeds[2] = id.event();
    seeds[3] = _salt;
    _engine->setSeeds( seeds.data(), seeds.size() );
  }

} // end namespace mu2e
