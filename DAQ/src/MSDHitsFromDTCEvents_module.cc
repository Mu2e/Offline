///////////////////////////////////////////////////////////////////////////////
// MSDHitsFromDTCEvents : add MSDHitCollection to the event
// M. Stortini
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/MobileSyncDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include "Offline/RecoDataProducts/inc/MSDHit.hh"

#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace art {
class MSDHitsFromDTCEvents;
}
// ======================================================================

class art::MSDHitsFromDTCEvents : public EDProducer {

public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> debugLevel{Name("debugLevel"), Comment("debug level")};
  };

  explicit MSDHitsFromDTCEvents(const art::EDProducer::Table<Config>& config);
  virtual ~MSDHitsFromDTCEvents() {}

  // --- overloaded functions of the art producer
  virtual void produce(art::Event& ArtEvent) override;
  virtual void beginRun(art::Run& ArtRun) override;

  //-----------------------------------------------------------------------------
  // helper functions
  //-----------------------------------------------------------------------------
  artdaq::Fragments getFragments(art::Event& event);

private:
  //-----------------------------------------------------------------------------
  // fcl parameters
  //-----------------------------------------------------------------------------
  int _debugLevel;
};

// ======================================================================
art::MSDHitsFromDTCEvents::MSDHitsFromDTCEvents(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, _debugLevel(config().debugLevel()) {
  produces<mu2e::MSDHitCollection>();
}

//-----------------------------------------------------------------------------
void art::MSDHitsFromDTCEvents::beginRun(art::Run& ArtRun) {}

// ----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
void art::MSDHitsFromDTCEvents::produce(Event& event) {

  // Collection of MSDHits for the event
  std::unique_ptr<mu2e::MSDHitCollection> msd_hits(new mu2e::MSDHitCollection);

  // get fragments then loop over them
  artdaq::Fragments fragments = getFragments(event);
  size_t index = 0;
  for (const auto& frag : fragments) {
    // from fragment get DTCEventFragment
    mu2e::DTCEventFragment eventFragment(frag);
    // from DTCEventFragment get the vector of DTC_SubEvents for the given subsystem (MSD)
    auto dtcSubEvents =
        eventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_MobileSync);
    if (_debugLevel > 2) {
      std::cout << "  Fragment " << index << " has " << dtcSubEvents.size()
                << " MSD DTC subevents, total sub events = "
                << eventFragment.getData().GetSubEvents().size() << std::endl;
      ++index;
    }
    // loop over the DTC_SubEvents
    size_t index_subevt = 0;
    for (auto& dtcSubEvent : dtcSubEvents) {
      // get the decoder
      mu2e::MobileSyncDataDecoder decoder(dtcSubEvent);
      // now loop over the blockIndex's
      if (_debugLevel > 3) {
        std::cout << "    Sub event " << index_subevt << " has " << decoder.block_count()
                  << " blocks\n";
        ++index_subevt;
      }
      for (size_t blockIndex = 0; blockIndex < decoder.block_count(); blockIndex++) {
        mu2e::MobileSyncDataDecoder::sync_data_t dataPacketsVec =
            decoder.GetMobileSyncPackets(blockIndex);
        // get number of hits in packet
        const unsigned nHitsInPacket = decoder.GetNHitsPerPacket();
        if (_debugLevel > 4) {
          std::cout << "      block " << blockIndex << " has " << dataPacketsVec.size()
                    << " packets\n";
        }
        // loop over dataPacketsVec and grab products of interest
        for (const auto& packet : dataPacketsVec) {
          // loop over hits in packet and fill msd_hits
          for (unsigned hit = 0; hit < nHitsInPacket; hit++) {
            double hitTime = 0;
            double hitTOT = 0;
            mu2e::MSDHit msdHit;
            if (decoder.GetHitTime(&packet, hit, hitTime)) {
              msdHit.setTime(hitTime);
              if (decoder.GetHitTOT(&packet, hit, hitTOT)) {
                msdHit.setTOT(hitTOT);
              }
              msd_hits->emplace_back(msdHit);
              if (_debugLevel > 0) {
                std::cout << "Hit " << hit << " in packet: time = " << msdHit.time()
                          << " ns, TOT = " << msdHit.tot() << " ns" << std::endl;
              }
            }
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------------
  // Store the msd hits
  //-----------------------------------------------------------------------------
  if (_debugLevel > 1) {
    std::cout << "[MSDHitsFromDTCEvents::" << __func__ << "] Event " << event.id() << " has "
              << msd_hits->size() << " MSD hits\n";
  }
  event.put(std::move(msd_hits));
}

// ----------------------------------------------------------------------------
// get art fragments from event
//-----------------------------------------------------------------------------
artdaq::Fragments art::MSDHitsFromDTCEvents::getFragments(art::Event& event) {

  artdaq::Fragments fragments;
  artdaq::FragmentPtrs containerFragments;

  std::vector<art::Handle<artdaq::Fragments>> fragmentHandles;
  fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();
  if (_debugLevel > 1) {
    std::cout << "[MSDHitsFromDTCEvents::" << __func__ << "] Event " << event.id() << " has "
              << fragmentHandles.size() << " fragment handles\n";
  }
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
  if (_debugLevel > 2) {
    std::cout << "[MSDHitsFromDTCEvents::" << __func__ << "] Found " << fragments.size()
              << " fragments\n";
  }
  return fragments;
}

// ======================================================================

DEFINE_ART_MODULE(art::MSDHitsFromDTCEvents)

// ======================================================================

