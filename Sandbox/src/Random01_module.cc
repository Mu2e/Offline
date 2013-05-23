//
// Module used to test the random number servce.
//
// $Id: Random01_module.cc,v 1.1 2013/05/23 16:37:37 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/05/23 16:37:37 $
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

#include "CLHEP/Random/RandFlat.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class Random01 : public art::EDProducer {
  public:

    explicit Random01(fhicl::ParameterSet const& pset);

    virtual void produce(art::Event& e);

  private:

    int                                            seed_;
    art::RandomNumberGenerator::base_engine_t&     engine_;
    CLHEP::RandFlat                                flat_;


  };

  Random01::Random01(fhicl::ParameterSet const& pset):
    seed_(pset.get<int>("seed")),
    engine_(createEngine(seed_)),
    flat_(engine_){
  }

  void Random01::produce( art::Event& event ) {
    cout << "Event: "
         << event.id().event() << " "
         << flat_.fire()
         << endl;
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::Random01);
