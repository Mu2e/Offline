//
//  Do-nothing dummy module, serves as a placeholder in a sequence to allow nulling
//  out a producer or to provide a slot for a real producer without having to redefine a sequence
//  Original author D. Brown (LBNL)
//
#include "art/Framework/Core/EDProducer.h"

namespace mu2e {
  class DummyProducer : public art::EDProducer {
    public:
      struct Config {};
      using Parameters = art::EDProducer::Table<Config>;
      explicit DummyProducer(const Parameters& conf);
      void produce(art::Event& evt) override;
  };
  DummyProducer::DummyProducer(const Parameters& config ) : art::EDProducer{config} {}
  void DummyProducer::produce(art::Event& ) {}
}
DEFINE_ART_MODULE(mu2e::DummyProducer)
