// ======================================================================
//
// StrawRecoFromFragments_plugin:  Add tracker data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include <artdaq-core/Data/Fragment.hh>

#include <iostream>

#include <string>

#include <memory>

namespace art {
class StrawDigisFromDataDecoders;
}

// ======================================================================

class art::StrawDigisFromDataDecoders : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<int> useTrkADC{fhicl::Name("useTrkADC"),
                               fhicl::Comment("parse tracker ADC waveforms")};
    fhicl::Atom<art::InputTag> trkTag{fhicl::Name("trkTag"), fhicl::Comment("trkTag")};
  };

  // --- C'tor/d'tor:
  explicit StrawDigisFromDataDecoders(const art::EDProducer::Table<Config>& config);
  virtual ~StrawDigisFromDataDecoders() {}

  // --- Production:
  virtual void produce(Event&);

private:
  void
  analyze_tracker_(const mu2e::TrackerDataDecoder& cc,
                   std::unique_ptr<mu2e::StrawDigiCollection> const& straw_digis,
                   std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> const& straw_digi_adcs);
  int diagLevel_;
  int useTrkADC_;

  art::InputTag trkFragmentsTag_;

  const int hexShiftPrint = 7;

}; // StrawDigisFromDataDecoders

// ======================================================================

art::StrawDigisFromDataDecoders::StrawDigisFromDataDecoders(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, diagLevel_(config().diagLevel()), useTrkADC_(config().useTrkADC()),
    trkFragmentsTag_(config().trkTag()) {
  produces<mu2e::StrawDigiCollection>();
  if (useTrkADC_) {
    produces<mu2e::StrawDigiADCWaveformCollection>();
  }
  produces<mu2e::IntensityInfoTrackerHits>();
  // FIXME!
  produces<mu2e::ProtonBunchTime>();
}

// ----------------------------------------------------------------------

void art::StrawDigisFromDataDecoders::produce(Event& event) {
  art::EventNumber_t eventNumber = event.event();

  // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);
  std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> straw_digi_adcs(
      new mu2e::StrawDigiADCWaveformCollection);
  // IntensityInfoTrackerHits
  std::unique_ptr<mu2e::IntensityInfoTrackerHits> intInfo(new mu2e::IntensityInfoTrackerHits);

  // FIXME! this is temporary
  std::unique_ptr<mu2e::ProtonBunchTime> pbt(new mu2e::ProtonBunchTime);
  pbt->pbtime_ = 0;
  pbt->pbterr_ = 0;
  event.put(std::move(pbt));

  size_t totalSize = 0;
  size_t numTrkFrags = 0;
  auto fragmentHandle = event.getValidHandle<std::vector<mu2e::TrackerDataDecoder> >(trkFragmentsTag_);

  for (auto frag : *fragmentHandle) {
    analyze_tracker_(frag, straw_digis, straw_digi_adcs);
    for (size_t i = 0; i < frag.block_count(); ++i) {
      totalSize += frag.blockSizeBytes(i);
    }
    numTrkFrags++;
  }

  if (numTrkFrags == 0) {
    std::cout << "[StrawDigisFromDataDecoders::produce] found no Tracker fragments!" << std::endl;
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numTrkFrags << " TRK fragments. " << std::endl;

    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::StrawDigisFromDataDecoders::produce exiting eventNumber="
              << (int)(event.event()) << " / timestamp=" << (int)eventNumber << std::endl;
  }

  intInfo->setNTrackerHits(straw_digis->size());
  event.put(std::move(intInfo));
  // Store the straw digis in the event
  event.put(std::move(straw_digis));
  if (useTrkADC_) {
    event.put(std::move(straw_digi_adcs));
  }

} // produce()

void art::StrawDigisFromDataDecoders::analyze_tracker_(
    const mu2e::TrackerDataDecoder& cc, std::unique_ptr<mu2e::StrawDigiCollection> const& straw_digis,
    std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> const& straw_digi_adcs) {

  if (diagLevel_ > 1) {
    std::cout << std::endl;
    std::cout << "TrackerDataDecoder: ";
    std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
    std::cout << std::endl;
    std::cout << "\t"
              << "====== Example Block Sizes ======" << std::endl;
    for (size_t i = 0; i < 10; i++) {
      if (i < cc.block_count()) {
        std::cout << "\t" << i << "\t" << cc.blockSizeBytes(i) << std::endl;
      }
    }
    std::cout << "\t"
              << "=========================" << std::endl;
  }

  for (size_t curBlockIdx = 0; curBlockIdx < cc.block_count(); curBlockIdx++) {
#if 0 // TODO: Review this code and update as necessary
    if (diagLevel_ > 1) {
      // Print binary contents the first 3 packets starting at the current position
      // In the case of the tracker simulation, this will be the whole tracker
      // DataBlock. In the case of the calorimeter, the number of data packets
      // following the header packet is variable.
      cc.printPacketAtByte(cc.blockIndexBytes(0) + 16 * (0 + 3 * curBlockIdx));
      cc.printPacketAtByte(cc.blockIndexBytes(0) + 16 * (1 + 3 * curBlockIdx));
      cc.printPacketAtByte(cc.blockIndexBytes(0) + 16 * (2 + 3 * curBlockIdx));

      // Print out decimal values of 16 bit chunks of packet data
      for (int i = hexShiftPrint; i >= 0; i--) {
        std::cout << "0x" << std::hex << std::setw(4) << std::setfill('0') << (adc_t) * (pos + i)
                  << std::dec << std::setw(0);
        std::cout << " ";
      }
      std::cout << std::endl;
    }
#endif

    auto block = cc.dataAtBlockIndex(curBlockIdx);
    if (block == nullptr) {
      mf::LogError("StrawDigisFromDataDecoders")
          << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (diagLevel_ > 1) {

      std::cout << "timestamp: "
                << static_cast<int>(hdr->GetEventWindowTag().GetEventWindowTag(true)) << std::endl;
      std::cout << "hdr->SubsystemID: " << static_cast<int>(hdr->GetSubsystemID()) << std::endl;
      std::cout << "dtcID: " << static_cast<int>(hdr->GetID()) << std::endl;
      std::cout << "rocID: " << static_cast<int>(hdr->GetLinkID()) << std::endl;
      std::cout << "packetCount: " << static_cast<int>(hdr->GetPacketCount()) << std::endl;
      std::cout << "EVB mode: " << static_cast<int>(hdr->GetEVBMode()) << std::endl;

      std::cout << std::endl;
    }

    // Parse phyiscs information from TRK packets
    if (hdr->GetPacketCount() > 0) {

      // Create the StrawDigi data products
      auto trkDataVec = cc.GetTrackerData(curBlockIdx, useTrkADC_);
      if (trkDataVec.empty()) {
        mf::LogError("StrawDigisFromDataDecoders")
            << "Error retrieving Tracker data from DataBlock " << curBlockIdx
            << "! Aborting processing of this block!";
        continue;
      }

      for (auto& trkDataPair : trkDataVec) {

        mu2e::StrawId sid(trkDataPair.first->StrawIndex);
        mu2e::TrkTypes::TDCValues tdc = {trkDataPair.first->TDC0(), trkDataPair.first->TDC1()};
        mu2e::TrkTypes::TOTValues tot = {trkDataPair.first->TOT0, trkDataPair.first->TOT1};
        mu2e::TrkTypes::ADCValue pmp = trkDataPair.first->PMP;

        // Fill the StrawDigiCollection
        straw_digis->emplace_back(sid, tdc, tot, pmp);
        if (useTrkADC_) {
          straw_digi_adcs->emplace_back(trkDataPair.second);
        }

        if (diagLevel_ > 1) {
          std::cout << "MAKEDIGI: " << sid.asUint16() << " " << tdc[0] << " " << tdc[1] << " "
                    << tot[0] << " " << tot[1] << " ";
          for (size_t i = 0; i < trkDataPair.second.size(); i++) {
            std::cout << trkDataPair.second[i];
            if (i < trkDataPair.second.size() - 1) {
              std::cout << " ";
            }
          }
          std::cout << std::endl;

          std::cout << std::endl;

          std::cout << "strawIdx: " << sid.asUint16() << std::endl;
          std::cout << "TDC0: " << tdc[0] << std::endl;
          std::cout << "TDC1: " << tdc[1] << std::endl;
          std::cout << "TOT0: " << tot[0] << std::endl;
          std::cout << "TOT1: " << tot[1] << std::endl;
          std::cout << "PMP:  " << pmp << std::endl;
          std::cout << "Waveform: {";
          for (size_t i = 0; i < trkDataPair.second.size(); i++) {
            std::cout << trkDataPair.second[i];
            if (i < trkDataPair.second.size() - 1) {
              std::cout << ",";
            }
          }
          std::cout << "}" << std::endl;

          std::cout << "FPGA Flags: ";
          for (size_t i = 8; i < 16; i++) {
            if (((0x0001 << (15 - i)) & trkDataPair.first->ErrorFlags) > 0) {
              std::cout << "1";
            } else {
              std::cout << "0";
            }
          }
          std::cout << std::endl;

          std::cout << "LOOP: " << hdr->GetEventWindowTag().GetEventWindowTag(true) << " "
                    << curBlockIdx << std::endl;

          // Text format: timestamp strawidx tdc0 tdc1 nsamples sample0-11
          // Example: 1 1113 36978 36829 12 1423 1390 1411 1354 2373 2392 2342 2254 1909 1611 1525
          // 1438
          std::cout << "GREPMETRK: " << hdr->GetEventWindowTag().GetEventWindowTag(true) << " ";
          std::cout << sid.asUint16() << " ";
          std::cout << tdc[0] << " ";
          std::cout << tdc[1] << " ";
          std::cout << tot[0] << " ";
          std::cout << tot[1] << " ";
          std::cout << pmp << " ";
          std::cout << trkDataPair.second.size() << " ";
          for (size_t i = 0; i < trkDataPair.second.size(); i++) {
            std::cout << trkDataPair.second[i];
            if (i < trkDataPair.second.size() - 1) {
              std::cout << " ";
            }
          }
          std::cout << std::endl;
        } // End debug output
      }
    }
  }

  // cc.ClearUpgradedPackets();
}

// ======================================================================

DEFINE_ART_MODULE(art::StrawDigisFromDataDecoders)

// ======================================================================
