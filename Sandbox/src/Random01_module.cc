//
// Module used to test the random number servce.
//
// Contact person Rob Kutschke
//
// Notes:
// 1) The module must create the engine even if the engine state is begin restored
//    from the event at the start of each event; if it is not created, the
//    service does not know the identity of the engine to which to restore the state.
//

#include "SeedService/inc/SeedService.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class Random01 : public art::EDAnalyzer {
  public:

    explicit Random01(fhicl::ParameterSet const& pset);

    virtual void analyze(art::Event const& e) override;

  private:

    std::string                                    myLabel_;
    int                                            seed_;
    art::RandomNumberGenerator::base_engine_t&     engine_;
    CLHEP::RandFlat                                flat_;


  };

  Random01::Random01(fhicl::ParameterSet const& pset):
    EDAnalyzer(pset),
    myLabel_(pset.get<std::string>("module_label")),
    seed_( art::ServiceHandle<SeedService>()->getSeed() ),
    engine_(createEngine(seed_)),
    flat_(engine_){
    mf::LogVerbatim("TEST") << "Constructor: "
                            << myLabel_
                            << " Seed: " << seed_
                            << endl;
  }

  void Random01::analyze( art::Event const& event ) {
    mf::LogVerbatim("TEST") << "Event: "
                            << myLabel_           << " "
                            << event.id().event() << " "
                            << flat_.fire()
                            << endl;
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::Random01);
