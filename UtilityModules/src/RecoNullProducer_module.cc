//
// A producer module that does nothing.  Similar to NullProducer_module.cc
// but with no refernces to the SeedService or the RandomEngine header.
//
// Use it to turn an existing module label into a no-op so that
// it is not necessary to remove a module from a path.
//
// Original author Rob Kutschke
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class RecoNullProducer : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
    };
    typedef art::EDProducer::Table<Config> Parameters;


    explicit RecoNullProducer(Parameters const& conf);

    void produce( art::Event& e) override;

  private:

  };

  RecoNullProducer::RecoNullProducer( Parameters const& conf ):
    art::EDProducer{conf}
    {}

  void RecoNullProducer::produce(art::Event& event) {
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RecoNullProducer)
