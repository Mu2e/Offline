// ======================================================================
//
// CaloRecoFromFragments_plugin:  Add cal data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

//-- insert calls to proditions ..for calodmap-----
#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
//-------------------------------------------------

#include "Offline/DAQ/inc/CaloDAQUtilities.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

#include <artdaq-core/Data/Fragment.hh>

#include <iostream>

#include <string>

#include <memory>

namespace art {
class CaloRecoFromFragments;
}

// ======================================================================

class art::CaloRecoFromFragments : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> data_type {fhicl::Name("dataType" ) , fhicl::Comment("Data type (0:standard, 1:debug, 2:counters)"), 0};
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<art::InputTag> caloTag{fhicl::Name("caloTag"), fhicl::Comment("caloTag")};
  };

  // --- C'tor/d'tor:
  explicit CaloRecoFromFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CaloRecoFromFragments() {}

  // --- Production:
  virtual void produce(Event&);

private:
  mu2e::ProditionsHandle<mu2e::CaloDAQMap> _calodaqconds_h;

  void analyze_calorimeter_(mu2e::CaloDAQMap const& calodaqconds,
                            const mu2e::CalorimeterDataDecoder& cc,
                            std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis);

  int data_type_;
  int diagLevel_;

  art::InputTag caloFragmentsTag_;
  mu2e::CaloDAQUtilities caloDAQUtil_;

  const int hexShiftPrint = 7;

}; // CaloRecoFromFragments
// ======================================================================
art::CaloRecoFromFragments::CaloRecoFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, 
    data_type_(config().data_type()),
    diagLevel_(config().diagLevel()),
    caloFragmentsTag_(config().caloTag()), caloDAQUtil_("CaloRecoFromFragments") {
  produces<mu2e::CaloDigiCollection>();
}
// ----------------------------------------------------------------------
void art::CaloRecoFromFragments::produce(Event& event) {
  art::EventNumber_t eventNumber = event.event();

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);

  mu2e::CaloDAQMap const& calodaqconds = _calodaqconds_h.get(event.id()); // Get calo daq cond

  size_t totalSize = 0;
  size_t numCalFrags = 0;
  auto fragmentHandle =
      event.getValidHandle<std::vector<mu2e::CalorimeterDataDecoder> >(caloFragmentsTag_);

  for (auto frag : *fragmentHandle) {
    analyze_calorimeter_(calodaqconds, frag, calo_digis);
    for (size_t i = 0; i < frag.block_count(); ++i) {
      totalSize += frag.blockSizeBytes(i);
    }
    numCalFrags++;
  }

  if (numCalFrags == 0) {
    std::cout << "[CaloRecoFromFragments::produce] found no Calorimeter fragments!" << std::endl;
    event.put(std::move(calo_digis));
    return;
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::CaloRecoFromFragments::produce exiting eventNumber=" << (int)(event.event())
              << " / timestamp=" << (int)eventNumber << std::endl;
  }

  // Store the calo digis in the event
  event.put(std::move(calo_digis));

} // produce()

void art::CaloRecoFromFragments::analyze_calorimeter_(
    mu2e::CaloDAQMap const& calodaqconds, const mu2e::CalorimeterDataDecoder& cc,
    std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis) {

  if (diagLevel_ > 1) {
    std::cout << std::endl;
    std::cout << "ArtFragmentReader: ";
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
      mf::LogError("CaloRecoFromFragments")
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

    if (hdr->GetPacketCount() > 0) { // Parse phyiscs information from CAL packets

      if (data_type_ == 0){ //Standard calo data

        auto calHitDataVec = cc.GetCalorimeterHitData(curBlockIdx);
        if (calHitDataVec == nullptr) {
          mf::LogError("CaloRecoFromFragments")
              << "Error retrieving Calorimeter data from block " << curBlockIdx
              << "! Aborting processing of this block!";
          continue;
        }
  
        bool err = false;
        for (unsigned int hitIdx = 0; hitIdx < calHitDataVec->size(); hitIdx++) {
          std::pair<mu2e::CalorimeterDataDecoder::CalorimeterHitDataPacket, std::vector<uint16_t>>  calData = calHitDataVec->at(hitIdx);
            // Fill the CaloDigiCollection
  
            if (diagLevel_ > 0) {
              std::cout << "[CaloRecoFromFragments] calo hit " << hitIdx << std::endl;
              std::cout << "[CaloRecoFromFragments] \tChNumber   "
                        << (int)calHitDataVec->at(hitIdx).first.ChannelNumber << std::endl;
              std::cout << "[CaloRecoFromFragments] \tDIRACA     " << (int)calHitDataVec->at(hitIdx).first.DIRACA
                        << std::endl;
              std::cout << "[CaloRecoFromFragments] \tDIRACB     " << (int)calHitDataVec->at(hitIdx).first.DIRACB
                        << std::endl;
              std::cout << "[CaloRecoFromFragments] \tErrorFlags " << (int)calHitDataVec->at(hitIdx).first.ErrorFlags
                        << std::endl;
              std::cout << "[CaloRecoFromFragments] \tTime              "
                        << (int)calHitDataVec->at(hitIdx).first.Time << std::endl;
              std::cout << "[CaloRecoFromFragments] \tNSamples   "
                        << (int)calHitDataVec->at(hitIdx).first.NumberOfSamples << std::endl;
              std::cout << "[CaloRecoFromFragments] \tIndexMax   "
                        << (int)calHitDataVec->at(hitIdx).first.IndexOfMaxDigitizerSample << std::endl;
            }
  
            uint16_t packetid = calHitDataVec->at(hitIdx).first.DIRACA;
            uint16_t dirac = packetid & 0xFF;
            uint16_t diracChannel = (packetid >>8) & 0x1F;
            mu2e::CaloRawSiPMId rawId(dirac,diracChannel);
            uint16_t roID = calodaqconds.offlineId(rawId).id();
            // uint16_t dettype = (packetId & 0x7000) >> 13;
  
            // FIXME: Can we match vector types here?
            std::vector<int> caloHits;
            caloHits.reserve(calHitDataVec->at(hitIdx).second.size());
            std::copy(calHitDataVec->at(hitIdx).second.begin(), calHitDataVec->at(hitIdx).second.end(),
                      std::back_inserter(caloHits));
  
            calo_digis->emplace_back(roID, calHitDataVec->at(hitIdx).first.Time, caloHits,
                                     calHitDataVec->at(hitIdx).first.IndexOfMaxDigitizerSample);
  
            if (diagLevel_ > 1) {
              // Until we have the final mapping, the BoardID is just a placeholder
              // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);
              uint16_t crystalID = roID / 2;
              std::cout << "Crystal ID: " << (int)crystalID << std::endl;
              std::cout << "SiPM ID: " << (int)roID << std::endl;
              std::cout << "Time: " << (int)calHitDataVec->at(hitIdx).first.Time << std::endl;
              std::cout << "NumSamples: " << (int)calHitDataVec->at(hitIdx).first.NumberOfSamples << std::endl;
              std::cout << "Waveform: {";
              for (size_t i = 0; i < calHitDataVec->at(hitIdx).second.size(); i++) {
                std::cout << calHitDataVec->at(hitIdx).second[i];
                if (i < calHitDataVec->at(hitIdx).second.size() - 1) {
                  std::cout << ",";
                }
              }
              std::cout << "}" << std::endl;
  
              // Text format: timestamp crystalID roID time nsamples samples...
              // Example: 1 201 402 660 18 0 0 0 0 1 17 51 81 91 83 68 60 58 52 42 33 23 16
              std::cout << "GREPMECAL: " << hdr->GetEventWindowTag().GetEventWindowTag(true) << " ";
              std::cout << crystalID << " ";
              std::cout << roID << " ";
              std::cout << calHitDataVec->at(hitIdx).first.Time << " ";
              std::cout << calHitDataVec->at(hitIdx).second.size() << " ";
              for (size_t i = 0; i < calHitDataVec->at(hitIdx).second.size(); i++) {
                std::cout << calHitDataVec->at(hitIdx).second[i];
                if (i < calHitDataVec->at(hitIdx).second.size() - 1) {
                  std::cout << " ";
                }
              }
              std::cout << std::endl;
            } // End debug output
  
          //} // End loop over readout channels in DataBlock
          if (err)
            continue;
        }
      } else if (data_type_ == 1){ //debug calo data

        auto calHitTestDataVec = cc.GetCalorimeterHitTestData(curBlockIdx);
        if (calHitTestDataVec == nullptr) {
          mf::LogError("CaloRecoFromFragments")
              << "Error retrieving Calorimeter test data from block " << curBlockIdx
              << "! Aborting processing of this block!";
          continue;
        }
  
        bool err = false;
        for (unsigned int hitIdx = 0; hitIdx < calHitTestDataVec->size(); hitIdx++) {
          std::pair<mu2e::CalorimeterDataDecoder::CalorimeterHitTestDataPacket, std::vector<uint16_t>>  calData = calHitTestDataVec->at(hitIdx);
            // Fill the CaloDigiCollection
  
            if (diagLevel_ > 0) {
              std::cout << "[CaloRecoFromFragments] calo hit " << hitIdx << std::endl;
              std::cout << "[CaloRecoFromFragments] \tBoardID   "
                        << (int)calHitTestDataVec->at(hitIdx).first.BoardID << std::endl;
              std::cout << "[CaloRecoFromFragments] \tChNumber   "
                        << (int)calHitTestDataVec->at(hitIdx).first.ChannelID << std::endl;
              std::cout << "[CaloRecoFromFragments] \tInPayloadEventWindowTag   "
                        << (int)calHitTestDataVec->at(hitIdx).first.InPayloadEventWindowTag << std::endl;

              std::cout << "[CaloRecoFromFragments] \tErrorFlags " << (int)calHitTestDataVec->at(hitIdx).first.ErrorFlags
                        << std::endl;
              std::cout << "[CaloRecoFromFragments] \tTime              "
                        << (int)calHitTestDataVec->at(hitIdx).first.Time << std::endl;
              std::cout << "[CaloRecoFromFragments] \tNSamples   "
                        << (int)calHitTestDataVec->at(hitIdx).first.NumberOfSamples << std::endl;
              std::cout << "[CaloRecoFromFragments] \tIndexMax   "
                        << (int)calHitTestDataVec->at(hitIdx).first.IndexOfMaxDigitizerSample << std::endl;
            }
  
            uint16_t packetid = 0; //calHitTestDataVec->at(hitIdx).first.DIRACA;
            uint16_t dirac = packetid & 0xFF;
            uint16_t diracChannel = (packetid >>8) & 0x1F;
            mu2e::CaloRawSiPMId rawId(dirac,diracChannel);
            uint16_t roID = calodaqconds.offlineId(rawId).id();
            // uint16_t dettype = (packetId & 0x7000) >> 13;
  
            std::vector<int> caloHits;
            caloHits.reserve(calHitTestDataVec->at(hitIdx).second.size());
            std::copy(calHitTestDataVec->at(hitIdx).second.begin(), calHitTestDataVec->at(hitIdx).second.end(),
                      std::back_inserter(caloHits));
  
            calo_digis->emplace_back(roID, uint(calHitTestDataVec->at(hitIdx).first.Time), caloHits,
                                     uint(calHitTestDataVec->at(hitIdx).first.IndexOfMaxDigitizerSample));
  
            if (diagLevel_ > 1) {
              // Until we have the final mapping, the BoardID is just a placeholder
              // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);
              uint16_t crystalID = roID / 2;
              std::cout << "Crystal ID: " << (int)crystalID << std::endl;
              std::cout << "SiPM ID: " << (int)roID << std::endl;
              std::cout << "Time: " << (int)calHitTestDataVec->at(hitIdx).first.Time << std::endl;
              std::cout << "NumSamples: " << (int)calHitTestDataVec->at(hitIdx).first.NumberOfSamples << std::endl;
              std::cout << "Waveform: {";
              for (size_t i = 0; i < calHitTestDataVec->at(hitIdx).second.size(); i++) {
                std::cout << calHitTestDataVec->at(hitIdx).second[i];
                if (i < calHitTestDataVec->at(hitIdx).second.size() - 1) {
                  std::cout << ",";
                }
              }
              std::cout << "}" << std::endl;
  
              // Text format: timestamp crystalID roID time nsamples samples...
              // Example: 1 201 402 660 18 0 0 0 0 1 17 51 81 91 83 68 60 58 52 42 33 23 16
              std::cout << "GREPMECAL: " << hdr->GetEventWindowTag().GetEventWindowTag(true) << " ";
              std::cout << crystalID << " ";
              std::cout << roID << " ";
              std::cout << calHitTestDataVec->at(hitIdx).first.Time << " ";
              std::cout << calHitTestDataVec->at(hitIdx).second.size() << " ";
              for (size_t i = 0; i < calHitTestDataVec->at(hitIdx).second.size(); i++) {
                std::cout << calHitTestDataVec->at(hitIdx).second[i];
                if (i < calHitTestDataVec->at(hitIdx).second.size() - 1) {
                  std::cout << " ";
                }
              }
              std::cout << std::endl;
            } // End debug output
  
          //} // End loop over readout channels in DataBlock
          if (err)
            continue;
        }

      
      }
    } // End Cal Mode
  }
}
// ======================================================================

DEFINE_ART_MODULE(art::CaloRecoFromFragments)

// ======================================================================
