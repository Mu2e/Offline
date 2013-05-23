//
// Module used to test the random number servce.
//
// $Id: Random01_module.cc,v 1.3 2013/05/23 19:15:56 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/05/23 19:15:56 $
//
// Contact person Rob Kutschke
//
// Notes:
// 1) The module must create the engine even if the engine state is begin restored
//    from the event at the start of each event; if it is not created, the
//    service does not know the identity of the engine to which to restore the state.
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class Random01 : public art::EDProducer {
  public:

    explicit Random01(fhicl::ParameterSet const& pset);

    virtual void produce(art::Event& e);

  private:

    std::string                                    myLabel_;
    int                                            seed_;
    art::RandomNumberGenerator::base_engine_t&     engine_;
    CLHEP::RandFlat                                flat_;


  };

  Random01::Random01(fhicl::ParameterSet const& pset):
    myLabel_(pset.get<std::string>("module_label")),
    seed_(pset.get<int>("seed")),
    engine_(createEngine(seed_)),
    flat_(engine_){
  }

  void Random01::produce( art::Event& event ) {
    mf::LogVerbatim("TEST") << "Event: "
                            << myLabel_           << " "
                            << event.id().event() << " "
                            << flat_.fire()
                            << endl;
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::Random01);
