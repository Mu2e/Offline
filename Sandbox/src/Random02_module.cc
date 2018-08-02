//
// A variant of Random01 in which the intermediate objects are not explicit.
//

#include "SeedService/inc/SeedService.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class Random02 : public art::EDAnalyzer {
  public:

    explicit Random02(fhicl::ParameterSet const& pset);

    virtual void analyze(art::Event const& e) override;

  private:

    std::string                                    myLabel_;
    CLHEP::HepRandomEngine&                        engine_;
    CLHEP::RandFlat                                flat1_;
    CLHEP::RandFlat                                flat2_;
  };

  Random02::Random02(fhicl::ParameterSet const& pset):
    EDAnalyzer(pset),
    myLabel_(pset.get<std::string>("module_label")),
    engine_{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    flat1_{engine_},
    flat2_{engine_, 1., 2.}
  {
    mf::LogVerbatim("TEST") << "Random02 Constructor: "
                            << myLabel_
                            << endl;
  }

  void Random02::analyze( art::Event const& event ) {
    mf::LogVerbatim("TEST") << "Random02 Event: "
                            << myLabel_           << " "
                            << event.id().event() << " "
                            << flat1_.fire()      << " "
                            << flat2_.fire()
                            << endl;
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::Random02);
