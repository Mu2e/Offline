// Counts the number of (presumably generated) events, and writes the
// result to SubRun.
//
// Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "cetlib_except/exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

#include "MCDataProducts/inc/GenEventCount.hh"

namespace mu2e {

  //================================================================
  class GenEventCounter : public art::EDProducer {
    GenEventCount::count_t seenEvents_;
  public:
    struct Config {};
    using Parameters = art::EDProducer::Table<Config>;
    explicit GenEventCounter(const Parameters& conf);
    virtual void produce(art::Event& event) override;
    virtual void endSubRun(art::SubRun& sr) override;
  };

  //================================================================
  GenEventCounter::GenEventCounter(const Parameters& conf)
    : art::EDProducer{conf},
      seenEvents_(0)
  {
    produces<mu2e::GenEventCount, art::InSubRun>();
  }

  //================================================================
  void GenEventCounter::produce(art::Event&) {
    ++seenEvents_;
  }

  //================================================================
  void GenEventCounter::endSubRun(art::SubRun& sr) {

    // The intention is to have at most one object of type
    // GenEventCount per SubRun.  Throw here if such an
    // object is already in SubRun.
    std::vector<art::Handle<GenEventCount> > hh = sr.getMany<GenEventCount>();
    if(!hh.empty()) {
      std::ostringstream os;
      os<<"GenEventCounter: refusing to write event count in "
        <<sr.id()<<" because conflicting products exist:\n";
      for(const auto& h : hh) {
        os<<"    moduleLabel = "<<h.provenance()->moduleLabel()
          <<", instance = "<<h.provenance()->productInstanceName()
          <<", process = "<<h.provenance()->processName()
          <<"\n";
      }
      os<<"\n";
      throw cet::exception("BADCONFIG")<<os.str();
    }

    mf::LogInfo("Summary")<<"Creating GenEventCount record: "<<seenEvents_
                          <<" events for "<<sr.id()<<"\n";

    sr.put(std::unique_ptr<GenEventCount>(new GenEventCount(seenEvents_)));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GenEventCounter);
