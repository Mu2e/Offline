//
// A producer module that does nothing.
// Use it to turn an existing module label into a no-op so that
// it is not necessary to remove a module from a path.
//
// Optionally create a specified number of random engines.
// This is needed if the module being dummied out consumes random seeds
// and it is important to maintain seed consistency in other modules.
//
// Original author Rob Kutschke
//

#include "SeedService/inc/SeedService.hh"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "CLHEP/Random/RandomEngine.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class NullProducer : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<unsigned> nEngines{Name("nEngines"),
          Comment("Number of engines to create"),0};
    };
    typedef art::EDProducer::Table<Config> Parameters;


    explicit NullProducer(Parameters const& conf);

    void produce( art::Event& e) override;

  private:

    // Non-owning pointers. Engines are owned by the RandomNumberGenerator service.
    std::vector<CLHEP::HepRandomEngine*> _engines;

  };

  NullProducer::NullProducer( Parameters const& conf ):
    art::EDProducer{conf}
    {
      _engines.reserve(conf().nEngines());
      for ( unsigned i=0; i<conf().nEngines(); ++i){
        _engines.emplace_back( &createEngine(art::ServiceHandle<SeedService>{}->getSeed()) );
      }
  }

  void NullProducer::produce(art::Event& event) {
  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::NullProducer)
