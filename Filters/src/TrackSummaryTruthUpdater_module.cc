// Update TrackSummaryTruthAssns to point to a filtered SimParticle collection
//
// Andrei Gaponenko, 2014

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "MCDataProducts/inc/TrackSummaryTruthAssns.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

namespace mu2e {

  class TrackSummaryTruthUpdater : public art::EDProducer {
  public:
    explicit TrackSummaryTruthUpdater(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt) override;
  private:
    art::InputTag inputTST_;
    art::InputTag remapping_;
  };

  //================================================================
  TrackSummaryTruthUpdater::TrackSummaryTruthUpdater(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , inputTST_(pset.get<std::string>("inputTST"))
    , remapping_(pset.get<std::string>("remapping"))
  {
    produces<TrackSummaryTruthAssns>();
  }

  //================================================================
  void TrackSummaryTruthUpdater::produce(art::Event& event) {

    std::unique_ptr<TrackSummaryTruthAssns> out(new TrackSummaryTruthAssns());

    auto ih = event.getValidHandle<TrackSummaryTruthAssns>(inputTST_);
    auto remap = event.getValidHandle<SimParticleRemapping>(remapping_);

    // redmine #6420:
    // for(typename TrackSummaryTruthAssns::size_type i=0; i<ih->size(); ++i) {
    // use "unsigned" for now
    for(unsigned i=0; i<ih->size(); ++i) {
      const auto iter = remap->find(ih->at(i).first);
      if(iter != remap->end()) {
        out->addSingle(iter->second , ih->at(i).second, ih->data(i));
      }
    }

    event.put(std::move(out));
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackSummaryTruthUpdater);
