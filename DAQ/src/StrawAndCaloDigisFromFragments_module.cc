// ======================================================================
//
// StrawAndCaloDigisFromFragments_plugin:  Add tracker/cal data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "mu2e-artdaq-core/Overlays/CalorimeterFragment.hh"
#include "mu2e-artdaq-core/Overlays/FragmentType.hh"
#include "mu2e-artdaq-core/Overlays/TrackerFragment.hh"

#include "DataProducts/inc/TrkTypes.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/ProtonBunchTime.hh"

#include <artdaq-core/Data/Fragment.hh>

#include <iostream>

#include <string>

#include <memory>

namespace art {
class StrawAndCaloDigisFromFragments;
}

// ======================================================================

class art::StrawAndCaloDigisFromFragments : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<int> parseCAL{fhicl::Name("parseCAL"), fhicl::Comment("parseCAL")};
    fhicl::Atom<int> parseTRK{fhicl::Name("parseTRK"), fhicl::Comment("parseTRK")};
    fhicl::Atom<int> useTrkADC{fhicl::Name("useTrkADC"), fhicl::Comment("parse tracker ADC waveforms")};
    fhicl::Atom<art::InputTag> caloTag{fhicl::Name("caloTag"), fhicl::Comment("caloTag")};
    fhicl::Atom<art::InputTag> trkTag{fhicl::Name("trkTag"), fhicl::Comment("trkTag")};
  };

  // --- C'tor/d'tor:
  explicit StrawAndCaloDigisFromFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawAndCaloDigisFromFragments() {}

  // --- Production:
  virtual void produce(Event&);

private:
  void analyze_tracker_(const artdaq::Fragment& f,
                        std::unique_ptr<mu2e::StrawDigiCollection> const& straw_digis,
                        std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> const& straw_digi_adcs);
  void analyze_calorimeter_(const artdaq::Fragment& f,
                            std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis);

  int diagLevel_;

  int parseCAL_;
  int parseTRK_;
  int useTrkADC_;

  art::InputTag trkFragmentsTag_;
  art::InputTag caloFragmentsTag_;

  const int hexShiftPrint = 7;

}; // StrawAndCaloDigisFromFragments

// ======================================================================

art::StrawAndCaloDigisFromFragments::StrawAndCaloDigisFromFragments(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_(config().diagLevel()), 
    parseCAL_(config().parseCAL()),
    parseTRK_(config().parseTRK()),
    useTrkADC_(config().useTrkADC()),
    trkFragmentsTag_(config().trkTag()),
    caloFragmentsTag_(config().caloTag()) {
  if (parseTRK_) {
    produces<mu2e::StrawDigiCollection>();
    if (useTrkADC_) {
      produces<mu2e::StrawDigiADCWaveformCollection>();
    }
  }
  if (parseCAL_) {
    produces<mu2e::CaloDigiCollection>();
  }
  //FIXME!
  produces<mu2e::ProtonBunchTime>();
}

// ----------------------------------------------------------------------

void art::StrawAndCaloDigisFromFragments::produce(Event& event) {
  art::EventNumber_t eventNumber = event.event();

  // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);
  std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> straw_digi_adcs(new mu2e::StrawDigiADCWaveformCollection);

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);

  //FIXME! this is temporary
  std::unique_ptr<mu2e::ProtonBunchTime> pbt(new mu2e::ProtonBunchTime);
  pbt->pbtime_ = 0;
  pbt->pbterr_ = 0;
  event.put(std::move(pbt));

  art::Handle<artdaq::Fragments> trkFragments, calFragments;
  size_t numTrkFrags(0), numCalFrags(0);
  size_t totalSize = 0;
  if (parseTRK_) {
    event.getByLabel(trkFragmentsTag_, trkFragments);
    if (!trkFragments.isValid()) {
      std::cout << "[StrawAndCaloDigisFromFragments::produce] found no Tracker fragments!"
                << std::endl;
      event.put(std::move(straw_digis));
      return;
    }
    numTrkFrags = trkFragments->size();
    for (size_t idx = 0; idx < numTrkFrags; ++idx) {
      auto size = ((*trkFragments)[idx]).sizeBytes(); // * sizeof(artdaq::RawDataType);
      totalSize += size;
      analyze_tracker_((*trkFragments)[idx], straw_digis, straw_digi_adcs);
      //      std::cout << "\tTRK Fragment " << idx << " has size " << size << std::endl;
    }
  }
  if (parseCAL_) {
    event.getByLabel(caloFragmentsTag_, calFragments);
    if (!calFragments.isValid()) {
      std::cout << "[StrawAndCaloDigisFromFragments::produce] found no Calorimeter fragments!"
                << std::endl;
      event.put(std::move(calo_digis));
      return;
    }
    numCalFrags = calFragments->size();
    for (size_t idx = 0; idx < numCalFrags; ++idx) {
      auto size = ((*calFragments)[idx]).sizeBytes(); // * sizeof(artdaq::RawDataType);
      totalSize += size;
      analyze_calorimeter_((*calFragments)[idx], calo_digis);
      //      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
    }
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numTrkFrags << " TRK fragments, and ";
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::StrawAndCaloDigisFromFragments::produce exiting eventNumber="
              << (int)(event.event()) << " / timestamp=" << (int)eventNumber << std::endl;
  }

  // Store the straw digis and calo digis in the event
  if (parseTRK_) {
    event.put(std::move(straw_digis));
    if (useTrkADC_) {
      event.put(std::move(straw_digi_adcs));
    }
  }
  if (parseCAL_) {
    event.put(std::move(calo_digis));
  }

} // produce()

void art::StrawAndCaloDigisFromFragments::analyze_tracker_(
    const artdaq::Fragment& f, std::unique_ptr<mu2e::StrawDigiCollection> const& straw_digis,
    std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> const& straw_digi_adcs) {

  mu2e::TrackerFragment cc(f);

  if (diagLevel_ > 1) {
    std::cout << std::endl;
    std::cout << "TrackerFragment: ";
    std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
    std::cout << "\tByte Count: " << f.dataSizeBytes() << std::endl;
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
      mf::LogError("StrawAndCaloDigisFromFragments")
          << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (diagLevel_ > 1) {

      std::cout << "timestamp: " << static_cast<int>(hdr.GetEventWindowTag().GetEventWindowTag(true))
                << std::endl;
      std::cout << "hdr->SubsystemID: " << static_cast<int>(hdr.GetSubsystemID()) << std::endl;
      std::cout << "dtcID: " << static_cast<int>(hdr.GetID()) << std::endl;
      std::cout << "rocID: " << static_cast<int>(hdr.GetLinkID()) << std::endl;
      std::cout << "packetCount: " << static_cast<int>(hdr.GetPacketCount()) << std::endl;
      std::cout << "EVB mode: " << static_cast<int>(hdr.GetEVBMode()) << std::endl;

      std::cout << std::endl;
    }

    // Parse phyiscs information from TRK packets
    if (hdr.GetPacketCount() > 0 && parseTRK_ > 0) {

      // Create the StrawDigi data products
      auto trkDataVec = cc.GetTrackerData(curBlockIdx);
      if (trkDataVec.empty()) {
        mf::LogError("StrawAndCaloDigisFromFragments")
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
        if (useTrkADC_){ 
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
          std::cout << "PMP:  " << pmp    << std::endl;
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

          std::cout << "LOOP: " << hdr.GetEventWindowTag().GetEventWindowTag(true) << " " << curBlockIdx
                    << std::endl;

          // Text format: timestamp strawidx tdc0 tdc1 nsamples sample0-11
          // Example: 1 1113 36978 36829 12 1423 1390 1411 1354 2373 2392 2342 2254 1909 1611 1525
          // 1438
          std::cout << "GREPMETRK: " << hdr.GetEventWindowTag().GetEventWindowTag(true) << " ";
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

  cc.ClearUpgradedPackets();
}

void art::StrawAndCaloDigisFromFragments::analyze_calorimeter_(
    const artdaq::Fragment& f, std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis) {

  mu2e::CalorimeterFragment cc(f);

  if (diagLevel_ > 1) {
    std::cout << std::endl;
    std::cout << "ArtFragmentReader: ";
    std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
    std::cout << "\tByte Count: " << f.dataSizeBytes() << std::endl;
    std::cout << std::endl;
    std::cout << "\t"
              << "====== Example Block Sizes ======" << std::endl;
    for (size_t i = 0; i < 10; i++) {
      if (i < cc.block_count()) {
        std::cout << "\t" << i  << "\t" << cc.blockSizeBytes(i)
                  << std::endl;
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
      mf::LogError("StrawAndCaloDigisFromFragments")
          << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (diagLevel_ > 1) {

      std::cout << "timestamp: " << static_cast<int>(hdr.GetEventWindowTag().GetEventWindowTag(true)) << std::endl;
      std::cout << "hdr->SubsystemID: " << static_cast<int>(hdr.GetSubsystemID()) << std::endl;
      std::cout << "dtcID: " << static_cast<int>(hdr.GetID()) << std::endl;
      std::cout << "rocID: " << static_cast<int>(hdr.GetLinkID()) << std::endl;
      std::cout << "packetCount: " << static_cast<int>(hdr.GetPacketCount()) << std::endl;
      std::cout << "EVB mode: " << static_cast<int>(hdr.GetEVBMode()) << std::endl;

      std::cout << std::endl;
    }

    if (hdr.GetPacketCount() > 0 &&
             parseCAL_ > 0) { // Parse phyiscs information from CAL packets

      auto calData = cc.GetCalorimeterData(curBlockIdx);
      if (calData == nullptr) {
        mf::LogError("StrawAndCaloDigisFromFragments")
            << "Error retrieving Calorimeter data from block " << curBlockIdx
            << "! Aborting processing of this block!";
        continue;
      }

      if (diagLevel_ > 0) {
        std::cout << "[StrawAndCaloDigiFromFragments] NEW CALDATA: NumberOfHits "
                  << calData->NumberOfHits << std::endl;
      }

      auto hits = cc.GetCalorimeterHits(curBlockIdx);
      bool err = false;
      for (size_t hitIdx = 0; hitIdx < calData->NumberOfHits; hitIdx++) {

        // Fill the CaloDigiCollection
        if (hitIdx > hits.size()) {
          mf::LogError("StrawAndCaloDigisFromFragments")
              << "Error retrieving Calorimeter data from block " << curBlockIdx << " for hit "
              << hitIdx << "! Aborting processing of this block!";
          err = true;
          break;
        }

        if (diagLevel_ > 0) {
          std::cout << "[StrawAndCaloDigiFromFragments] calo hit " << hitIdx << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tChNumber   "
                    << (int)hits[hitIdx].first.ChannelNumber
                    << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tDIRACA     "
                    << (int)hits[hitIdx].first.DIRACA
                    << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tDIRACB     "
                    << (int)hits[hitIdx].first.DIRACB
                    << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tErrorFlags "
                    << (int)hits[hitIdx].first.ErrorFlags
                    << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tTime	      "
                    << (int)hits[hitIdx].first.Time
                    << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tNSamples   "
                    << (int)hits[hitIdx].first.NumberOfSamples << std::endl;
          std::cout << "[StrawAndCaloDigiFromFragments] \tIndexMax   "
                    << (int)hits[hitIdx].first.IndexOfMaxDigitizerSample << std::endl;
        }

        // IMPORTANT NOTE: we don't have a final
        // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
        // placeholder. Because we still need to know which crystal a hit belongs to, we are
        // temporarily storing the 4-bit apdID and 12-bit crystalID in the Reserved DIRAC A slot.
        // Also, note that until we have an actual map, channel index does not actually correspond
        // to the physical readout channel on a ROC.
        uint16_t crystalID = hits[hitIdx].first.DIRACB & 0x0FFF;
        uint16_t apdID = hits[hitIdx].first.DIRACB >> 12;

        // FIXME: Can we match vector types here?
        std::vector<int> caloHits;
	caloHits.reserve(hits[hitIdx].second.size());
        std::copy(hits[hitIdx].second.begin(), hits[hitIdx].second.end(), std::back_inserter(caloHits));

        calo_digis->emplace_back((crystalID * 2 + apdID), hits[hitIdx].first.Time, caloHits,
                                 hits[hitIdx].first.IndexOfMaxDigitizerSample);

        if (diagLevel_ > 1) {
          // Until we have the final mapping, the BoardID is just a placeholder
          // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);

          std::cout << "Crystal ID: " << (int)crystalID << std::endl;
          std::cout << "APD ID: " << (int)apdID << std::endl;
          std::cout << "Time: " << (int)hits[hitIdx].first.Time << std::endl;
          std::cout << "NumSamples: " << (int)hits[hitIdx].first.NumberOfSamples << std::endl;
          std::cout << "Waveform: {";
          for (size_t i = 0; i < hits[hitIdx].second.size(); i++) {
            std::cout << hits[hitIdx].second[i];
            if (i < hits[hitIdx].second.size() - 1) {
              std::cout << ",";
            }
          }
          std::cout << "}" << std::endl;

          // Text format: timestamp crystalID roID time nsamples samples...
          // Example: 1 201 402 660 18 0 0 0 0 1 17 51 81 91 83 68 60 58 52 42 33 23 16
          std::cout << "GREPMECAL: " << hdr.GetEventWindowTag().GetEventWindowTag(true) << " ";
          std::cout << crystalID << " ";
          std::cout << apdID << " ";
          std::cout << hits[hitIdx].first.Time << " ";
          std::cout << hits[hitIdx].second.size() << " ";
          for (size_t i = 0; i < hits[hitIdx].second.size(); i++) {
            std::cout << hits[hitIdx].second[i];
            if (i < hits[hitIdx].second.size() - 1) {
              std::cout << " ";
            }
          }
          std::cout << std::endl;
        } // End debug output

      } // End loop over readout channels in DataBlock
      if (err)
        continue;

    } // End Cal Mode
  }
}
// ======================================================================

DEFINE_ART_MODULE(art::StrawAndCaloDigisFromFragments)

// ======================================================================
