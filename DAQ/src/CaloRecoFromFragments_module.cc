// ======================================================================
//
// CaloRecoFromFragments_plugin:  Add cal data products to the event
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

//-- insert calls to proditions ..for calodmap-----
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "CaloConditions/inc/CaloDAQMap.hh"
//-------------------------------------------------

#include "RecoDataProducts/inc/CaloDigi.hh"
#include "DAQ/inc/CaloDAQUtilities.hh"

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
			    const artdaq::Fragment& f,
                            std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis);

  int                        diagLevel_;

  art::InputTag             caloFragmentsTag_;
  mu2e::CaloDAQUtilities    caloDAQUtil_;

  const int                 hexShiftPrint = 7;

}; // CaloRecoFromFragments
// ======================================================================
art::CaloRecoFromFragments::CaloRecoFromFragments(const art::EDProducer::Table<Config>& config) :
  art::EDProducer{config},
  diagLevel_(config().diagLevel()), 
  caloFragmentsTag_(config().caloTag()),
  caloDAQUtil_("CaloRecoFromFragments") {
    produces<mu2e::CaloDigiCollection>();
  }
// ----------------------------------------------------------------------
void art::CaloRecoFromFragments::produce(Event& event) {
  art::EventNumber_t eventNumber = event.event();
  
  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);

  mu2e::CaloDAQMap const& calodaqconds = _calodaqconds_h.get(event.id()); // Get calo daq cond
  
  art::Handle<artdaq::Fragments> calFragments;
  size_t numCalFrags(0);
  size_t totalSize = 0;
  event.getByLabel(caloFragmentsTag_, calFragments);
  if (!calFragments.isValid()) {
    std::cout << "[CaloRecoFromFragments::produce] found no Calorimeter fragments!"
	      << std::endl;
    event.put(std::move(calo_digis));
    return;
  }
  
  numCalFrags = calFragments->size();
  if (diagLevel_ > 1) {
    std::cout << "[CaloRecoFromFragments::produce] found "<< numCalFrags <<" Calorimeter fragments"
	      << std::endl;
  }
  for (size_t idx = 0; idx < numCalFrags; ++idx) {
    auto size = ((*calFragments)[idx]).sizeBytes(); // * sizeof(artdaq::RawDataType);
    totalSize += size;
    analyze_calorimeter_(calodaqconds, (*calFragments)[idx], calo_digis);
    //      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
  }
  

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::CaloRecoFromFragments::produce exiting eventNumber="
              << (int)(event.event()) << " / timestamp=" << (int)eventNumber << std::endl;
  }

  // Store the calo digis in the event
  event.put(std::move(calo_digis));


} // produce()

void art::CaloRecoFromFragments::analyze_calorimeter_(mu2e::CaloDAQMap const& calodaqconds,
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
      mf::LogError("CaloRecoFromFragments")
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

    if (hdr.GetPacketCount() > 0) { // Parse phyiscs information from CAL packets

      auto calData = cc.GetCalorimeterData(curBlockIdx);
      if (calData == nullptr) {
        mf::LogError("CaloRecoFromFragments")
            << "Error retrieving Calorimeter data from block " << curBlockIdx
            << "! Aborting processing of this block!";
        continue;
      }

      if (diagLevel_ > 0) {
        std::cout << "[CaloRecoFromFragments] NEW CALDATA: NumberOfHits "
                  << calData->NumberOfHits << std::endl;
      }

      auto hits = cc.GetCalorimeterHits(curBlockIdx);
      bool err = false;
      for (size_t hitIdx = 0; hitIdx < calData->NumberOfHits; hitIdx++) {

        // Fill the CaloDigiCollection
        if (hitIdx > hits.size()) {
          mf::LogError("CaloRecoFromFragments")
              << "Error retrieving Calorimeter data from block " << curBlockIdx << " for hit "
              << hitIdx << "! Aborting processing of this block!";
          err = true;
          break;
        }

        if (diagLevel_ > 0) {
          std::cout << "[CaloRecoFromFragments] calo hit " << hitIdx << std::endl;
          std::cout << "[CaloRecoFromFragments] \tChNumber   "
                    << (int)hits[hitIdx].first.ChannelNumber
                    << std::endl;
          std::cout << "[CaloRecoFromFragments] \tDIRACA     "
                    << (int)hits[hitIdx].first.DIRACA
                    << std::endl;
          std::cout << "[CaloRecoFromFragments] \tDIRACB     "
                    << (int)hits[hitIdx].first.DIRACB
                    << std::endl;
          std::cout << "[CaloRecoFromFragments] \tErrorFlags "
                    << (int)hits[hitIdx].first.ErrorFlags
                    << std::endl;
          std::cout << "[CaloRecoFromFragments] \tTime	      "
                    << (int)hits[hitIdx].first.Time
                    << std::endl;
          std::cout << "[CaloRecoFromFragments] \tNSamples   "
                    << (int)hits[hitIdx].first.NumberOfSamples << std::endl;
          std::cout << "[CaloRecoFromFragments] \tIndexMax   "
                    << (int)hits[hitIdx].first.IndexOfMaxDigitizerSample << std::endl;
        }

        // IMPORTANT NOTE: we don't have a final
        // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
        // placeholder. Because we still need to know which crystal a hit belongs to, we are
        // temporarily storing the 4-bit roID and 12-bit crystalID in the Reserved DIRAC A slot.
        // Also, note that until we have an actual map, channel index does not actually correspond
        // to the physical readout channel on a ROC.
        //uint16_t crystalID = hits[hitIdx].first.DIRACB & 0x0FFF;
        //uint16_t roID = hits[hitIdx].first.DIRACB >> 12;

	uint16_t packetid = hits[hitIdx].first.DIRACA;
	uint16_t roID = calodaqconds.packetIdTocaloRoId(packetid);	
	//uint16_t dettype = (packetId & 0x7000) >> 13;
							
        // FIXME: Can we match vector types here?
        std::vector<int> caloHits;
	caloHits.reserve(hits[hitIdx].second.size());
        std::copy(hits[hitIdx].second.begin(), hits[hitIdx].second.end(), std::back_inserter(caloHits));

        calo_digis->emplace_back(roID, hits[hitIdx].first.Time, caloHits,
                                 hits[hitIdx].first.IndexOfMaxDigitizerSample);

        if (diagLevel_ > 1) {
          // Until we have the final mapping, the BoardID is just a placeholder
          // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);
	  uint16_t crystalID = roID/2;
          std::cout << "Crystal ID: " << (int)crystalID << std::endl;
          std::cout << "SiPM ID: " << (int)roID << std::endl;
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
          std::cout << roID << " ";
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

DEFINE_ART_MODULE(art::CaloRecoFromFragments)

// ======================================================================
