// ======================================================================
//
// ArtFragmentsFromDTCEvents_plugin:  Add cal data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CRVFragment.hh"
#include "artdaq-core-mu2e/Data/CalorimeterFragment.hh"
#include "artdaq-core-mu2e/Data/TrackerFragment.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>

#include <string>

#include <array>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

namespace art {
class ArtFragmentsFromDTCEvents;
}

// ======================================================================

class art::ArtFragmentsFromDTCEvents : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<float> makeCaloFrag{fhicl::Name("makeCaloFrag"),
                                    fhicl::Comment("Create ArtCaloFragmets")};
    fhicl::Atom<float> makeTrkFrag{fhicl::Name("makeTrkFrag"),
                                   fhicl::Comment("Create ArtTrkFragmets")};
    fhicl::Atom<float> makeCRVFrag{fhicl::Name("makeCRVFrag"),
                                   fhicl::Comment("Create ArtCRVFragmets")};
    fhicl::Atom<float> makeSTMFrag{fhicl::Name("makeSTMFrag"),
                                   fhicl::Comment("Create ArtSTMFragmets")};
  };

  // --- C'tor/d'tor:
  explicit ArtFragmentsFromDTCEvents(const art::EDProducer::Table<Config>& config);
  virtual ~ArtFragmentsFromDTCEvents() {}

  virtual void beginRun(art::Run&) override;

  // --- Production:
  virtual void produce(Event&);

private:
  int diagLevel_;
  int makeCaloFrag_, makeTrkFrag_, makeCRVFrag_, makeSTMFrag_;
};

// ======================================================================

void art::ArtFragmentsFromDTCEvents::beginRun(art::Run& Run) {}

art::ArtFragmentsFromDTCEvents::ArtFragmentsFromDTCEvents(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_(config().diagLevel()), makeCaloFrag_(config().makeCaloFrag()),
    makeTrkFrag_(config().makeTrkFrag()), makeCRVFrag_(config().makeCRVFrag()),
    makeSTMFrag_(config().makeSTMFrag()) {

  if (config().makeCaloFrag() > 0) {
    produces<std::vector<mu2e::CalorimeterFragment>>();
  }
  if (config().makeTrkFrag() > 0) {
    produces<std::vector<mu2e::TrackerFragment>>();
  }
  if (config().makeCRVFrag() > 0) {
    produces<std::vector<mu2e::CRVFragment>>();
  }
  // if (config().makeSTMFrag() > 0) { produces<std::vector<mu2e::STMFragment>>();}
}

// ----------------------------------------------------------------------

void art::ArtFragmentsFromDTCEvents::produce(Event& event) {

  art::EventNumber_t eventNumber = event.event();

  // Collection of CaloHits for the event
  std::unique_ptr<std::vector<mu2e::CalorimeterFragment>> caloFragColl(
      new std::vector<mu2e::CalorimeterFragment>);
  std::unique_ptr<std::vector<mu2e::TrackerFragment>> trkFragColl(
      new std::vector<mu2e::TrackerFragment>);
  std::unique_ptr<std::vector<mu2e::CRVFragment>> crvFragColl(new std::vector<mu2e::CRVFragment>);

  // std::unique_ptr<std::vector<mu2e::STMFragment>>         stm_frags(new
  // std::vector<mu2e::STMFragment>); if (makeSTMFrag_ > 0) {

  artdaq::Fragments fragments;
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

  if (diagLevel_ > 0) {
    std::cout << "[ArtFragmentsFromDTCEvents::produce] Found nHandlesnFragments  "
              << fragments.size() << std::endl;
  }

  size_t nFrags(0);

  for (const auto& frag : fragments) {
    mu2e::DTCEventFragment bb(frag);

    if (makeTrkFrag_ > 0) { // TRACKER
      auto trkSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);
      for (auto const& subevent : trkSEvents) {
        mu2e::TrackerFragment tf(subevent);
        trkFragColl->emplace_back(tf);
        ++nFrags;
      }
    }

    if (makeCaloFrag_ > 0) { // CALORIMETER
      auto caloSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter);
      for (auto& subevent : caloSEvents) {
        mu2e::CalorimeterFragment cf(subevent);
        caloFragColl->emplace_back(cf);
        ++nFrags;
      }
    }

    if (makeCRVFrag_ > 0) { // CRV
      auto crvSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_CRV);
      for (auto& subevent : crvSEvents) {
        mu2e::CRVFragment cf(subevent);
        crvFragColl->emplace_back(cf);
        ++nFrags;
      }

      //FIXME: Temporary implementation until the DTC header gets fixed.
      //Currently, the DTC header uses the Subsystem ID for the tracker.
      //Checking the TDAQ header of the first data block instead.
      auto crvSEventsT = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);
      for(auto& subevent : crvSEventsT)
      {
        mu2e::CRVFragment cf(subevent);
        cf.setup_event();
        if(cf.block_count()>0)
        {
          auto block = cf.dataAtBlockIndex(0);
          if(block == nullptr) continue;
          auto header = block->GetHeader();
          if(header->GetSubsystemID() != DTCLib::DTC_Subsystem::DTC_Subsystem_CRV) continue;

          crvFragColl->emplace_back(cf);
          ++nFrags;
        }
      }
    }
  }

  if (nFrags == 0) {
    std::cout << "[ArtFragmentsFromDTCEvents::produce] found no fragments!" << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::ArtFragmentsFromDTCEvents::produce exiting eventNumber="
              << (int)(event.event()) << " / timestamp=" << (int)eventNumber << std::endl;
  }

  if (makeCaloFrag_ > 0) {
    event.put(std::move(caloFragColl));
  }
  if (makeTrkFrag_ > 0) {
    event.put(std::move(trkFragColl));
  }
  if (makeCRVFrag_ > 0) {
    event.put(std::move(crvFragColl));
  }
  //  if (makeSTMFrag_ > 0)  {event.put(std::move(stmFragColl));}

} // produce()

// ======================================================================

DEFINE_ART_MODULE(art::ArtFragmentsFromDTCEvents)

// ======================================================================
