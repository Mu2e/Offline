// ======================================================================
//
// CaloDigiFromFragments_plugin:  Add cal data products to the event
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
#include <map>

#include <memory>

namespace art {
class CaloDigiFromFragments;
}

// ======================================================================

class art::CaloDigiFromFragments : public EDProducer {

public:
  struct Config {
    fhicl::Atom<int> data_type {fhicl::Name("dataType" ) , fhicl::Comment("Data type (0:standard, 1:debug, 2:counters)"), 0};
    fhicl::Atom<int> diagLevel {fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level"), 0};
    fhicl::Atom<art::InputTag> caloTag {fhicl::Name("caloTag"), fhicl::Comment("Input module")};
    fhicl::Atom<bool> useOfflineID {fhicl::Name("useOfflineID"), fhicl::Comment("Use calo disk mapping for SiPM IDs (Default: true)"), true};
    fhicl::Atom<bool> useDTCROCID {fhicl::Name("useDTCROCID"), fhicl::Comment("Use DTC and ROC numbers instead of boardID for SiPM IDs (Default: false)"), false};
  };

  // --- C'tor/d'tor:
  explicit CaloDigiFromFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CaloDigiFromFragments() {}

  // --- Production:
  virtual void produce(Event&);
  virtual void endJob();

private:
  mu2e::ProditionsHandle<mu2e::CaloDAQMap> _calodaqconds_h;

  void analyze_calorimeter_(mu2e::CaloDAQMap const& calodaqconds,
                            const mu2e::CalorimeterDataDecoder& cc,
                            std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis);

  int data_type_;
  int diagLevel_;

  art::InputTag caloFragmentsTag_;
  mu2e::CaloDAQUtilities caloDAQUtil_;
  bool useOfflineID_;
  bool useDTCROCID_;

  const int hexShiftPrint = 7;

  long int total_events;
  long int total_hits;
  long int total_hits_good;
  long int total_hits_bad;
  std::map<mu2e::CaloDAQUtilities::CaloHitError, uint> failure_counter;

}; // CaloDigiFromFragments
// ======================================================================
art::CaloDigiFromFragments::CaloDigiFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, 
    data_type_(config().data_type()),
    diagLevel_(config().diagLevel()),
    caloFragmentsTag_(config().caloTag()),
    caloDAQUtil_("CaloDigiFromFragments"),
    useOfflineID_(config().useOfflineID()),
    useDTCROCID_(config().useDTCROCID()) {
  produces<mu2e::CaloDigiCollection>();
  total_events = 0;
  total_hits = 0;
  total_hits_good = 0;
  total_hits_bad = 0;
}
// ----------------------------------------------------------------------
void art::CaloDigiFromFragments::produce(Event& event) {
  art::EventNumber_t eventNumber = event.event();
  total_events++;

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);

  mu2e::CaloDAQMap const& calodaqconds = _calodaqconds_h.get(event.id()); // Get calo daq cond

  size_t totalSize = 0;
  size_t numCalDecoders = 0;
  auto decoderColl = event.getValidHandle<std::vector<mu2e::CalorimeterDataDecoder> >(caloFragmentsTag_);

  for (auto decoder : *decoderColl) {
    analyze_calorimeter_(calodaqconds, decoder, calo_digis);
    for (size_t i = 0; i < decoder.block_count(); ++i) {
      totalSize += decoder.blockSizeBytes(i);
    }
    numCalDecoders++;
  }

  if (numCalDecoders == 0) {
    std::cout << "[CaloDigiFromFragments::produce] found no Calorimeter decoders!" << std::endl;
    event.put(std::move(calo_digis));
    return;
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "[CaloDigiFromFragments::produce] Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << numCalDecoders << " CALO decoders." << std::endl;
    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  // Store the calo digis in the event
  event.put(std::move(calo_digis));

} // produce()

void art::CaloDigiFromFragments::analyze_calorimeter_(
    mu2e::CaloDAQMap const& calodaqconds, const mu2e::CalorimeterDataDecoder& cc,
    std::unique_ptr<mu2e::CaloDigiCollection> const& calo_digis) {

  auto dtcID = cc.event_.GetDTCID();

  //Loop through the ROCs of this caloDecoder
  for (size_t iROC = 0; iROC < cc.block_count(); iROC++) {

    if (data_type_ == 0){ //Standard calo data

      auto calHitDataVec = cc.GetCalorimeterHitData(iROC);
      if (calHitDataVec == nullptr) {
        mf::LogError("CaloDigiFromFragments")
            << "Error retrieving Calorimeter data from ROC " << iROC
            << "! Aborting processing of this ROC!";
        continue;
      }
  
      //Loop through the hits of this ROC
      for (unsigned int hitIdx = 0; hitIdx < calHitDataVec->size(); hitIdx++) {
        mu2e::CalorimeterDataDecoder::CalorimeterHitDataPacket& thisHitPacket = calHitDataVec->at(hitIdx).first;
        std::vector<uint16_t>& thisHitWaveform = calHitDataVec->at(hitIdx).second;

        // Fill the CaloDigiCollection

        if (diagLevel_ > 1) {
          std::cout << "[CaloDigiFromFragments] calo hit " << hitIdx << std::endl;
          std::cout << "[CaloDigiFromFragments] \tChNumber   "
                    << (int)thisHitPacket.ChannelNumber << std::endl;
          std::cout << "[CaloDigiFromFragments] \tDIRACA     " << (int)thisHitPacket.DIRACA
                    << std::endl;
          std::cout << "[CaloDigiFromFragments] \tDIRACB     " << (int)thisHitPacket.DIRACB
                    << std::endl;
          std::cout << "[CaloDigiFromFragments] \tErrorFlags " << (int)thisHitPacket.ErrorFlags
                    << std::endl;
          std::cout << "[CaloDigiFromFragments] \tTime              "
                    << (int)thisHitPacket.Time << std::endl;
          std::cout << "[CaloDigiFromFragments] \tNSamples   "
                    << (int)thisHitPacket.NumberOfSamples << std::endl;
          std::cout << "[CaloDigiFromFragments] \tIndexMax   "
                    << (int)thisHitPacket.IndexOfMaxDigitizerSample << std::endl;
        }

        uint16_t packetid = thisHitPacket.DIRACA;
        uint16_t dirac = packetid & 0xFF;
        uint16_t diracChannel = (packetid >>8) & 0x1F;
        mu2e::CaloRawSiPMId rawId(dirac,diracChannel);
        uint16_t roID = calodaqconds.offlineId(rawId).id();
        // uint16_t dettype = (packetId & 0x7000) >> 13;
  
        calo_digis->emplace_back(roID, thisHitPacket.Time, std::vector<int>(thisHitWaveform.begin(), thisHitWaveform.end()),
                                   thisHitPacket.IndexOfMaxDigitizerSample);

        if (diagLevel_ > 1) {
          // Until we have the final mapping, the BoardID is just a placeholder
          // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);
          uint16_t crystalID = roID / 2;
          std::cout << "Crystal ID: " << (int)crystalID << std::endl;
          std::cout << "SiPM ID: " << (int)roID << std::endl;
          std::cout << "Time: " << (int)thisHitPacket.Time << std::endl;
          std::cout << "NumSamples: " << (int)thisHitPacket.NumberOfSamples << std::endl;
          std::cout << "Waveform: {";
          for (auto sample : thisHitWaveform) {
            std::cout << sample << " ";
          }
          std::cout << "}" << std::endl;
        } // End debug output

      }
    } else if (data_type_ == 1){ //debug calo data

      auto calHitTestDataVec = cc.GetCalorimeterHitTestData(iROC);
      if (calHitTestDataVec == nullptr) {
        mf::LogError("CaloDigiFromFragments")
            << "Error retrieving Calorimeter test data from block " << iROC
            << "! Aborting processing of this block!";
        continue;
      }

      //Loop through the hits of this ROC
      total_hits += calHitTestDataVec->size();
      for (unsigned int hitIdx = 0; hitIdx < calHitTestDataVec->size(); hitIdx++) {

        mu2e::CalorimeterDataDecoder::CalorimeterHitTestDataPacket& thisHitPacket = calHitTestDataVec->at(hitIdx).first;
        std::vector<uint16_t>& thisHitWaveform = calHitTestDataVec->at(hitIdx).second;

        //Check if hit was decoded correctly
        auto errorCode = caloDAQUtil_.isHitGood(calHitTestDataVec->at(hitIdx));
        if (errorCode){
          failure_counter[errorCode]++;
          total_hits_bad++;
          if (diagLevel_ > 0){
            std::cout << "[CaloDigiFromFragments] BAD calo hit! DTC: " << dtcID << ", ROC: " << iROC << ", hit number: " << hitIdx << " [failure code: " << errorCode << "]" << std::endl;
            caloDAQUtil_.printCaloPulse(thisHitPacket);
            std::cout << "[CaloDigiFromFragments] \twaveform size \t" << thisHitWaveform.size() << std::endl;
          }
          continue;
        }

        if (diagLevel_ > 1) {
          std::cout << "[CaloDigiFromFragments] calo hit " << hitIdx << std::endl;
          caloDAQUtil_.printCaloPulse(thisHitPacket);
        }
        
        // Fill the CaloDigiCollection
        mu2e::CaloRawSiPMId rawId(thisHitPacket.BoardID,thisHitPacket.ChannelID);
        uint16_t SiPMID = ( useOfflineID_ ? calodaqconds.offlineId(rawId).id() : rawId.id() );
        if (useDTCROCID_) SiPMID = dtcID*120 + iROC*20 + thisHitPacket.ChannelID;
        //Constructor: CaloDigi(int SiPMID, int t0, const std::vector<int>& waveform, size_t peakpos)
        total_hits_good++;
        calo_digis->emplace_back(SiPMID, int(thisHitPacket.Time), std::vector<int>(thisHitWaveform.begin(), thisHitWaveform.end()),
                                 uint(thisHitPacket.IndexOfMaxDigitizerSample));
        if (diagLevel_ > 1) {
          uint16_t crystalID = SiPMID / 2;
          std::cout << "Crystal ID: " << (int)crystalID << std::endl;
          std::cout << "SiPM ID: " << (int)SiPMID << std::endl;
          std::cout << "Time: " << (int)thisHitPacket.Time << std::endl;
          std::cout << "NumSamples: " << (int)thisHitPacket.NumberOfSamples << std::endl;
          std::cout << "Waveform: {";
          for (auto sample : thisHitWaveform) {
            std::cout << sample << " ";
          }
          std::cout << "}" << std::endl;
        } // End debug output
      }

    } // End loop over hits
  } // End loop over ROCs
}


void art::CaloDigiFromFragments::endJob(){

  std::cout << "\n ----- [CaloDigiFromFragments] Decoding errors summary ----- " << std::endl;
  std::cout << "Total events: " << total_events << std::endl;
  std::cout << "Total hits: " << total_hits << std::endl;
  std::cout << "Total good hits: " << total_hits_good << std::endl;
  std::cout << "Total bad hits: " << total_hits_bad << std::endl;
  for (auto fail : failure_counter){
    std::cout << "Failure mode " << fail.first
      << " [bad " << caloDAQUtil_.getCaloHitErrorName(fail.first)
      << "], count: " << fail.second
      << " (" << int(100.*fail.second/total_hits) << "%)\n";
  }

}

// ======================================================================

DEFINE_ART_MODULE(art::CaloDigiFromFragments)

// ======================================================================
