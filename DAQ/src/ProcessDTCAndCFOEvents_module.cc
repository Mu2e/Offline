// ======================================================================
//
// Add data-decoders and EventHeader to the event from DTC_Events
// and CFO_Events
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
#include <artdaq-core-mu2e/Data/EventHeader.hh>
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/RecoDataProducts/inc/DAQerror.hh"

#include <iostream>

#include <string>

#include <array>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

namespace art {
class ProcessDTCAndCFOEvents;
}

// ======================================================================

class art::ProcessDTCAndCFOEvents : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<int> makeCaloFrag{fhicl::Name("makeCaloFrag"),
                                   fhicl::Comment("Create CaloDataDecoders")};
    fhicl::Atom<int> makeTrkFrag{fhicl::Name("makeTrkFrag"),
                                   fhicl::Comment("Create TrackerDataDecoders")};
    fhicl::Atom<int> makeCRVFrag{fhicl::Name("makeCRVFrag"),
                                   fhicl::Comment("Create CRVDataDecoders")};
  };

  // --- C'tor/d'tor:
  explicit ProcessDTCAndCFOEvents(const art::EDProducer::Table<Config>& config);
  ~ProcessDTCAndCFOEvents() override {}

  void beginRun(art::Run&) override;

  // --- Production:
  void produce(Event&) override;

private:
  int diagLevel_;
  int makeCaloFrag_, makeTrkFrag_, makeCRVFrag_;
};

// ======================================================================

void art::ProcessDTCAndCFOEvents::beginRun(art::Run& Run) {}

art::ProcessDTCAndCFOEvents::ProcessDTCAndCFOEvents(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_(config().diagLevel()), makeCaloFrag_(config().makeCaloFrag()),
    makeTrkFrag_(config().makeTrkFrag()), makeCRVFrag_(config().makeCRVFrag()) {

  if (config().makeCaloFrag() > 0) {
    produces<std::vector<mu2e::CalorimeterDataDecoder>>();
  }
  if (config().makeTrkFrag() > 0) {
    produces<std::vector<mu2e::TrackerDataDecoder>>();
  }
  if (config().makeCRVFrag() > 0) {
    produces<std::vector<mu2e::CRVDataDecoder>>();
  }

  produces<mu2e::EventHeader>();
  produces<mu2e::DAQerrorCollection>();
}

// ----------------------------------------------------------------------

void art::ProcessDTCAndCFOEvents::produce(Event& event) {

  art::EventNumber_t eventNumber = event.event();

  // Collection of CaloHits for the event
  std::unique_ptr<std::vector<mu2e::CalorimeterDataDecoder>> caloFragColl(
      new std::vector<mu2e::CalorimeterDataDecoder>);
  std::unique_ptr<std::vector<mu2e::TrackerDataDecoder>> trkFragColl(
      new std::vector<mu2e::TrackerDataDecoder>);
  std::unique_ptr<std::vector<mu2e::CRVDataDecoder>> crvFragColl(new std::vector<mu2e::CRVDataDecoder>);
  std::unique_ptr<mu2e::EventHeader> evtHdr(new mu2e::EventHeader);
  std::unique_ptr<mu2e::DAQerrorCollection> daqErrors(new mu2e::DAQerrorCollection);

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
    std::cout << "[ProcessDTCAndCFOEvents::produce] Found nHandlesnFragments  "
              << fragments.size() << std::endl;
  }

  size_t nFrags(0);

  for (size_t fragIdx=0; const auto& frag : fragments) {
    mu2e::DTCEventFragment bb(frag);

    if (makeTrkFrag_ > 0) { // TRACKER
      auto trkSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);
      for (auto const& subevent : trkSEvents) {
        trkFragColl->emplace_back(subevent);
        ++nFrags;
      }
    }

    if (makeCaloFrag_ > 0) { // CALORIMETER
      auto caloSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter);
      for (auto& subevent : caloSEvents) {
        caloFragColl->emplace_back(subevent);
        ++nFrags;
      }
    }

    if (makeCRVFrag_ > 0) { // CRV
      auto crvSEvents = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_CRV);
      for (auto& subevent : crvSEvents) {
        crvFragColl->emplace_back(subevent);
        ++nFrags;
      }
      auto crvSEventsTmp = bb.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker);  //currently wrongly encoded in the DTC Subevent header
      for (auto& subevent : crvSEventsTmp) {
        crvFragColl->emplace_back(subevent);
        ++nFrags;
      }
    }

    const DTCLib::DTC_Event &dtcEvent =  bb.getData();
    const std::vector<DTCLib::DTC_SubEvent> &dtcSubEvents = dtcEvent.GetSubEvents();
    size_t expectedSize = dtcEvent.GetEventByteCount();
    size_t actualSize = sizeof(DTCLib::DTC_EventHeader);
    for(size_t iSubEvent=0; iSubEvent<dtcSubEvents.size(); ++iSubEvent) actualSize+=dtcSubEvents.at(iSubEvent).GetSubEventByteCount();
    if(diagLevel_ > 0)
    {
      std::cout << "[ProcessDTCAndCFOEvents::produce] expected event size: " << expectedSize << ", actual event size: " << actualSize << std::endl;
    }
    if(expectedSize!=actualSize)
    {
      std::cerr << "[ProcessDTCAndCFOEvents::produce] mismatch between expected event size and actual event size!" << std::endl;
      daqErrors->emplace_back(mu2e::DAQerrorCode::byteCountMismatch,fragIdx);
    }
  }

  if ( (diagLevel_ > 0) && (nFrags == 0)) {
    std::cout << "[ProcessDTCAndCFOEvents::produce] found no fragments!" << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::ProcessDTCAndCFOEvents::produce exiting eventNumber="
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
  event.put(std::move(evtHdr));
  event.put(std::move(daqErrors));
} // produce()

// ======================================================================

DEFINE_ART_MODULE(art::ProcessDTCAndCFOEvents)

// ======================================================================
