// ======================================================================
//
// CrvDigisFromFragments_plugin:  Add CRV data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "mu2e-artdaq-core/Overlays/CRVFragment.hh"

#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <iostream>

#include <string>

#include <memory>

#define LONG_FROM_CRV 0 // Whether to print long-form data or single-line summary

namespace art {
class CrvDigisFromFragments;
}

using art::CrvDigisFromFragments;

// ======================================================================

class art::CrvDigisFromFragments : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level")};
    fhicl::Atom<art::InputTag> crvFragmentsTag{fhicl::Name("crvTag"),
                                               fhicl::Comment("crv Fragments Tag")};
  };

  // --- C'tor/d'tor:
  explicit CrvDigisFromFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CrvDigisFromFragments() {}

  // --- Production:
  virtual void produce(Event&);

private:
  int decompressCrvDigi(uint8_t adc);

  int diagLevel_;

  art::InputTag crvFragmentsTag_;

}; // CrvDigisFromFragments

// ======================================================================

CrvDigisFromFragments::CrvDigisFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, diagLevel_(config().diagLevel()),
    crvFragmentsTag_(config().crvFragmentsTag()) {
  produces<EventNumber_t>();
  produces<mu2e::CrvDigiCollection>();
}

// ----------------------------------------------------------------------

int CrvDigisFromFragments::decompressCrvDigi(uint8_t adc) {
  // TODO: Temporary implementation until we have the real compression used at the FEBs
  int toReturn = adc;
  if (adc >= 50 && adc < 75)
    toReturn = (adc - 50) * 2 + 50;
  if (adc >= 75 && adc < 100)
    toReturn = (adc - 75) * 4 + 100;
  if (adc >= 100 && adc < 125)
    toReturn = (adc - 100) * 8 + 200;
  if (adc >= 125)
    toReturn = (adc - 125) * 16 + 400;
  toReturn += 95;
  return toReturn;
}

void CrvDigisFromFragments::produce(Event& event) {

  art::EventNumber_t eventNumber = event.event();

  auto crvFragments = event.getValidHandle<artdaq::Fragments>(crvFragmentsTag_);
  size_t numCrvFrags = crvFragments->size();

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << crvFragments->size() << " CRV fragments." << std::endl;

    size_t totalSize = 0;
    for (size_t idx = 0; idx < crvFragments->size(); ++idx) {
      auto size = ((*crvFragments)[idx]).size() * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tCRV Fragment " << idx << " has size " << size << std::endl;
    }

    std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;
  }

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CrvDigiCollection> crv_digis(new mu2e::CrvDigiCollection);

  // Loop over the CRV fragments
  for (size_t idx = 0; idx < numCrvFrags; ++idx) {

    const auto& fragment((*crvFragments)[idx]);
    mu2e::CRVFragment cc(fragment);

    if (diagLevel_ > 1) {
      std::cout << std::endl;
      std::cout << "ArtFragmentReader: ";
      std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
      std::cout << "\tByte Count: " << fragment.dataSizeBytes() << std::endl;
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

#if 0 // TODO: Look into restoring some form of this later...
      if (diagLevel_ > 1) {
        // Print binary contents the first 3 packets starting at the current position
        // In the case of the tracker simulation, this will be the whole tracker
        // DataBlock. In the case of the calorimeter, the number of data packets
        // following the header packet is variable.
        cc.printPacketAtByte(cc.blockIndexBytes(0) + 16 * (0 + 3 * curBlockIdx));
        cc.printPacketAtByte(cc.blockIndexBytes(0) + 16 * (1 + 3 * curBlockIdx));
        cc.printPacketAtByte(cc.blockIndexBytes(0) + 16 * (2 + 3 * curBlockIdx));

        // Print out decimal values of 16 bit chunks of packet data
        for (int i = 7; i >= 0; i--) {
          std::cout << (adc_t) * (pos + i);
          std::cout << " ";
        }
        std::cout << std::endl;
      }
#endif

      auto block = cc.dataAtBlockIndex(curBlockIdx);
      if (block == nullptr) {
        std::cerr << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
        continue;
      }
      auto hdr = block->GetHeader();

      if (hdr.GetSubsystemID() != 2) {
        throw cet::exception("DATA") << " CRV packet does not have system ID 2";
      }

      // Parse phyiscs information from the CRV packets
      if (hdr.GetPacketCount() > 0) {
        auto crvRocHdr = cc.GetCRVROCStatusPacket(curBlockIdx);
        if (crvRocHdr == nullptr) {
          std::cerr << "Error retrieving CRV ROC Status Packet from DataBlock " << curBlockIdx
                    << "!" << std::endl;
          continue;
        }

        auto crvHits = cc.GetCRVHitReadoutPackets(curBlockIdx);
        for (auto const& crvHit : crvHits) {

          // Fill the CrvDigiCollection
          // CrvDigi(const std::array<unsigned int, NSamples> &ADCs, unsigned int startTDC,
          //         mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
          // TODO: This is a temporary implementation.
          // There will be a major change on the barIndex+SiPMNumber system,
          // which will be replaced by a channel ID system
          // Only a toy model is used here. The real implementation will follow.
          int channel = crvHit.SiPMID & 0x7F; // right 7 bits
          int FEB = crvHit.SiPMID >> 7;
          int crvBarIndex = (FEB * 64 + channel) / 4;
          int SiPMNumber = (FEB * 64 + channel) % 4;

          std::array<unsigned int, 8> adc;
          for (int j = 0; j < 8; j++)
            adc[j] = decompressCrvDigi(crvHit.Waveform().at(j));
          crv_digis->emplace_back(adc, crvHit.HitTime, mu2e::CRSScintillatorBarIndex(crvBarIndex),
                                  SiPMNumber);
        }

        if (diagLevel_ > 1) {

          for (auto const& crvHit : crvHits) {

              #if LONG_FORM_CRV
            // TODO: This is a temporary implementation.
            // There will be a major change on the barIndex+SiPMNumber system,
            // which will be replaced by a channel ID system
            // Only a toy model is used here. The real implementation will follow.
            int channel = crvHit.SiPMID & 0x7F; // right 7 bits
            int FEB = crvHit.SiPMID >> 7;
            int crvBarIndex = (FEB * 64 + channel) / 4;
            int SiPMNumber = (FEB * 64 + channel) % 4;

            std::cout << "MAKEDIGI: " << SiPMNumber << " " << crvBarIndex << " " << crvHit.HitTime
                      << " " << crvHits.size() << " ";

            auto hits = crvHit.Waveform();
            for (size_t j = 0; j < hits.size(); j++) {
              std::cout << decompressCrvDigi(hits[j]);
              if (j < hits.size() - 1) {
                std::cout << " ";
              }
            }
            std::cout << std::endl;

            std::cout << "timestamp: " << hdr.GetTimestamp().GetTimestamp(true) << std::endl;
            std::cout << "hdr->SubsystemID: " << hdr.GetSubsystemID() << std::endl;
            std::cout << "hdr->DTCID: " << hdr.GetID() << std::endl;
            std::cout << "rocID: " << hdr.GetRingID() << std::endl;
            std::cout << "packetCount: " << hdr.GetPacketCount() << std::endl;
            std::cout << "EVB mode: " << hdr.GetEVBMode() << std::endl;

            //	    for(int i=7; i>=0; i--) {
            //	      std::cout << (adc_t) *(pos+8+i);
            //	      std::cout << " ";
            //	    }
            //	    std::cout << std::endl;
            //
            //	    for(int i=7; i>=0; i--) {
            //	      std::cout << (adc_t) *(pos+8*2+i);
            //	      std::cout << " ";
            //	    }
            //	    std::cout << std::endl;

            std::cout << "SiPMNumber: " << crvHit.SiPMID % 4 << std::endl;
            std::cout << "scintillatorBarIndex: " << crvHit.SiPMID / 4 << std::endl;
            std::cout << "TDC: " << crvHit.HitTime << std::endl;
            std::cout << "Waveform: {";
            for (size_t j = 0; j < hits.size(); j++) {
              std::cout << hits[j];
              if (j < hits.size() - 1) {
                std::cout << " ";
              }
            }
            std::cout << "}" << std::endl;
            #else
            // Text format: timestamp sipmID tdc nsamples sample_list
            std::cout << "GREPMECRV: " << hdr.GetTimestamp().GetTimestamp(true) << " ";
            std::cout << crvHit.SiPMID << " ";
            std::cout << crvHit.HitTime << " ";
            auto hits = crvHit.Waveform();
            for (size_t j = 0; j < hits.size(); j++) {
              std::cout << hits[j];
              if (j < hits.size() - 1) {
                std::cout << " ";
              }
            }
            std::cout << std::endl;
            #endif
          }

          std::cout << "LOOP: " << eventNumber << " " << curBlockIdx << " "
                    << "(" << hdr.GetTimestamp().GetTimestamp(true) << ")" << std::endl;

        } // End debug output
      }   // End parsing CRV packets
    }     // End loop over DataBlocks within fragment
  }       // Close loop over fragments

  if (diagLevel_ > 0) {
    std::cout << "mu2e::CrvDigisFromFragments::produce exiting eventNumber=" << (int)(event.event())
              << " / timestamp=" << (int)eventNumber << std::endl;
  }

  event.put(std::unique_ptr<EventNumber_t>(new EventNumber_t(eventNumber)));

  // Store the straw digis and calo digis in the event
  event.put(std::move(crv_digis));

} // produce()

// ======================================================================

DEFINE_ART_MODULE(CrvDigisFromFragments)

// ======================================================================
