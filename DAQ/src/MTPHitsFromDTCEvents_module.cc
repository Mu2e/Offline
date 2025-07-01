///////////////////////////////////////////////////////////////////////////////
// MTPHitsFromDTCEvents : add MTPHitCollection to the event
// M. Stortini
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/Decoders/MTPDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include "Offline/RecoDataProducts/inc/MTPHit.hh"

#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include <iostream>
#include <string>
#include <map>
#include <memory>

namespace art {
  class MTPHitsFromDTCEvents;
}
// ======================================================================

class art::MTPHitsFromDTCEvents : public EDProducer {

public:

  struct Config {
    fhicl::Atom<int>   debugLevel        {fhicl::Name("debugLevel"    ), fhicl::Comment("debug level"           )};
    fhicl::Atom<float> clockFrequency    {fhicl::Name("clockFrequency"), fhicl::Comment("clock frequency in MHz")};
  };

  explicit MTPHitsFromDTCEvents(const art::EDProducer::Table<Config>& config);
  virtual ~MTPHitsFromDTCEvents() {}

  // --- overloaded functions of the art producer
  virtual void produce (art::Event& ArtEvent) override;
  virtual void beginRun(art::Run&   ArtRun  ) override;

  //-----------------------------------------------------------------------------
  // helper functions
  //-----------------------------------------------------------------------------
  artdaq::Fragments getFragments(art::Event& event);

private:

  //-----------------------------------------------------------------------------
  // fcl parameters
  //-----------------------------------------------------------------------------
  int      _debugLevel;
  float    _clockFrequency;

};

// ======================================================================
art::MTPHitsFromDTCEvents::MTPHitsFromDTCEvents(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _debugLevel    (config().debugLevel    ()),
    _clockFrequency(config().clockFrequency())
{
  produces<mu2e::MTPHitCollection>();
}


//-----------------------------------------------------------------------------
void art::MTPHitsFromDTCEvents::beginRun(art::Run&  ArtRun) {}

// ----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
void art::MTPHitsFromDTCEvents::produce(Event& event) {

  // Collection of MTPHits for the event
  std::unique_ptr<mu2e::MTPHitCollection> mtp_hits(new mu2e::MTPHitCollection);

  // get fragments then loop over them
  artdaq::Fragments fragments = getFragments(event);
  for (const auto& frag : fragments) {
    // from fragment get DTCEventFragment
    mu2e::DTCEventFragment eventFragment(frag);
    // from DTCEventFragment get the vector of DTC_SubEvents for the given subsystem (MTP)
    auto dtcSubEvents = eventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_MTP);
    // loop over the DTC_SubEvents
    for (auto& dtcSubEvent : dtcSubEvents) {
      // get the decoder
      mu2e::MTPDataDecoder decoder(dtcSubEvent);
      // now loop over the blockIndex's
      for (size_t blockIndex = 0; blockIndex < decoder.block_count(); blockIndex++) {
        mu2e::MTPDataDecoder::mtp_data_t dataPacketsVec = decoder.GetMTPDataPackets(blockIndex);
        // loop over dataPacketsVec and grab products of interest
        for (size_t vecIndex = 0; vecIndex < dataPacketsVec.size(); vecIndex++) {
          // grap MTPDataPacket and initialize MTPHit, grabbing timeStamp0 and timeStamp1
          const mu2e::MTPDataDecoder::MTPDataPacket* packet = dataPacketsVec.at(vecIndex);
          // grab the two time stamp counters, convert them to ns, and save them
          int channelID = 0; // not in payload yet, will be in future
          uint16_t counter0 = packet->GetTimestamp(0);
          float time0 = counter0*1000.0/_clockFrequency;
          mu2e::MTPHit mtpHit0(time0, channelID);
          mtp_hits->emplace_back(mtpHit0);
          uint16_t counter1 = packet->GetTimestamp(1);
          float time1 = counter1*1000.0/_clockFrequency;
          mu2e::MTPHit mtpHit1(time1, channelID);
          mtp_hits->emplace_back(mtpHit1);
          if (_debugLevel == 1) { std::cout << "time0, time1 = " << time0 << ", " << time1 << std::endl; }
        }
      }
    }
  }

//-----------------------------------------------------------------------------
// Store the mtp hits
//-----------------------------------------------------------------------------
  event.put(std::move(mtp_hits));

}

// ----------------------------------------------------------------------------
// get art fragments from event
//-----------------------------------------------------------------------------
artdaq::Fragments art::MTPHitsFromDTCEvents::getFragments(art::Event& event) {

  artdaq::Fragments    fragments;
  artdaq::FragmentPtrs containerFragments;

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();
  for (const auto& handle : fragmentHandles) {
    if (!handle.isValid() || handle->empty()) {
      continue;
    }

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      for (const auto& cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        if (contf.fragment_type() != mu2e::FragmentType::DTCEVT) {
          break;
        }

        for (size_t ii = 0; ii < contf.block_count(); ++ii) {
          containerFragments.push_back(contf[ii]);
          fragments.push_back(*containerFragments.back());
        }
      }
    } else {
      if (handle->front().type() == mu2e::FragmentType::DTCEVT) {
        for (auto frag : *handle) {
          fragments.emplace_back(frag);
        }
      }
    }
  }
  return fragments;
}

// ======================================================================

DEFINE_ART_MODULE(art::MTPHitsFromDTCEvents)

// ======================================================================
