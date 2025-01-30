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
#include "artdaq-core-mu2e/Data/CRVDataDecoder.hh"
#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
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
    fhicl::Atom<int> makeCaloFrag{fhicl::Name("makeCaloFrag"),
                                   fhicl::Comment("Create ArtCaloFragmets")};
    fhicl::Atom<int> makeTrkFrag{fhicl::Name("makeTrkFrag"),
                                   fhicl::Comment("Create ArtTrkFragmets")};
    fhicl::Atom<int> makeCRVFrag{fhicl::Name("makeCRVFrag"),
                                   fhicl::Comment("Create ArtCRVFragmets")};
    fhicl::Atom<int> makeSTMFrag{fhicl::Name("makeSTMFrag"),
                                   fhicl::Comment("Create ArtSTMFragmets")};
  };

  // --- C'tor/d'tor:
  explicit ArtFragmentsFromDTCEvents(const art::EDProducer::Table<Config>& config);
  ~ArtFragmentsFromDTCEvents() override {}

  void beginRun(art::Run&) override;

  // --- Production:
  void produce(Event&) override;

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
    produces<std::vector<mu2e::CalorimeterDataDecoder>>();
  }
  if (config().makeTrkFrag() > 0) {
    produces<std::vector<mu2e::TrackerDataDecoder>>();
  }
  if (config().makeCRVFrag() > 0) {
    produces<std::vector<mu2e::CRVDataDecoder>>();
  }
  // if (config().makeSTMFrag() > 0) { produces<std::vector<mu2e::STMFragment>>();}
}

// ----------------------------------------------------------------------

void art::ArtFragmentsFromDTCEvents::produce(Event& event) {

  art::EventNumber_t eventNumber = event.event();

  // Collection of CaloHits for the event
  std::unique_ptr<std::vector<mu2e::CalorimeterDataDecoder>> caloFragColl(
      new std::vector<mu2e::CalorimeterDataDecoder>);
  std::unique_ptr<std::vector<mu2e::TrackerDataDecoder>> trkFragColl(
      new std::vector<mu2e::TrackerDataDecoder>);
  std::unique_ptr<std::vector<mu2e::CRVDataDecoder>> crvFragColl(new std::vector<mu2e::CRVDataDecoder>);

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
        mu2e::TrackerDataDecoder tf(subevent);
        trkFragColl->emplace_back(tf);
        ++nFrags;
      }
    }

    if (makeCaloFrag_ > 0) { // CALORIMETER
      auto caloSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter);
      for (auto& subevent : caloSEvents) {
        mu2e::CalorimeterDataDecoder cf(subevent);
        caloFragColl->emplace_back(cf);
        ++nFrags;
      }
    }

    if (makeCRVFrag_ > 0) { // CRV
      auto crvSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_CRV);
      for (auto& subevent : crvSEvents) {
        mu2e::CRVDataDecoder cf(subevent);
        crvFragColl->emplace_back(cf);
        ++nFrags;
      }
      auto crvSEventsTmp = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);  //currently wrongly encoded in the DTC Subevent header
      for (auto& subevent : crvSEventsTmp) {
        mu2e::CRVDataDecoder cf(subevent);
        crvFragColl->emplace_back(cf);
        ++nFrags;
      }
    }
  }

  if ( (diagLevel_ > 0) && (nFrags == 0)) {
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
