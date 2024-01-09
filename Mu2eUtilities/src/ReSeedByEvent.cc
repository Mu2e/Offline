#include "Offline/Mu2eUtilities/inc/ReSeedByEvent.hh"
//#include "art/Framework/EventProcessor/Scheduler.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <iostream>
#include <string>

namespace mu2e {
  ReSeedByEvent::ReSeedByEvent( CLHEP::HepRandomEngine& engine, SeedService::seed_t salt):_salt(salt){
    std::string engineType{art::ServiceHandle<art::RandomNumberGenerator>{}->defaultEngineKind()};
    std::cout << "Kind:    " << engineType << std::endl;
    if ( engineType != "MixMaxRng" ){
      // throw
    }
    //std::cout << "Threads: " << art::ServiceHandle<art::Scheduler>{}->num_threads() << std::endl;
    _engine = static_cast<CLHEP::MixMaxRng*>(&engine);
  }

  void  ReSeedByEvent::reseed ( art::EventID const& id ) const{
    std::cout << "reseed: " << id  << std::endl;
    std::array<long,4> seeds;
    seeds[0] = id.run();
    seeds[1] = id.subRun();
    seeds[2] = id.event();
    seeds[3] = _salt;
    _engine->setSeeds( seeds.data(), seeds.size() );
  }

} // end namespace mu2e
