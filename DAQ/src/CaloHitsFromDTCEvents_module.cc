// ======================================================================
//
// Make fast CaloHits directly from events
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/Decoders/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include <artdaq-core/Data/Fragment.hh>

#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"

#include "Offline/DAQ/inc/CaloDAQUtilities.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"

#include <iostream>

#include <string>

#include <array>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

namespace art {
class CaloHitsFromDataDTCEvents;
}

// ======================================================================

class art::CaloHitsFromDataDTCEvents : public EDProducer {

  struct CrystalInfo {
    CrystalInfo(int crID, int a, float& b, float& c) :
        _crystalID(crID), _nSiPM(a), _time(b), _eDep(c) {}
    size_t _crystalID;
    size_t _nSiPM;
    float _time;
    float _eDep;
  };

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<float> digiSampling{fhicl::Name("digiSampling"),
                                    fhicl::Comment("calorimeter sampling period in ns")};
    fhicl::Atom<float> deltaTPulses{
        fhicl::Name("deltaTPulses"),
        fhicl::Comment(
            "time-gate between two signal from different SiPMs coupled with the same crystal")};
    fhicl::Atom<float> nPEperMeV{fhicl::Name("nPEperMeV"),
                                 fhicl::Comment("number of photo-electrons per MeV")};
    fhicl::Atom<float> noiseLevelMeV{fhicl::Name("noiseLevelMeV"),
                                     fhicl::Comment("Noise level in MeV")};
    fhicl::Atom<float> nSigmaNoise{fhicl::Name("nSigmaNoise"),
                                   fhicl::Comment("Maxnumber of sigma Noise to combine digi")};
    fhicl::Atom<float> hitEDepMax{
        fhicl::Name("hitEDepMax"),
        fhicl::Comment("Maximum hit energy in MeV (to reject saturated pulses)")};
    fhicl::Atom<float> hitEDepMin{fhicl::Name("hitEDepMin"),
                                  fhicl::Comment("Minimum hit energy in MeV")};
    fhicl::Atom<float> caphriEDepMax{fhicl::Name("caphriEDepMax"),
                                     fhicl::Comment("Maximum CAPHRI hit energy in MeV")};
    fhicl::Atom<float> caphriEDepMin{fhicl::Name("caphriEDepMin"),
                                     fhicl::Comment("Minimum CAPHRI hit energy in MeV")};
  };

  // --- C'tor/d'tor:
  explicit CaloHitsFromDataDTCEvents(const art::EDProducer::Table<Config>& config);
  virtual ~CaloHitsFromDataDTCEvents() {}

  virtual void beginRun(art::Run&) override;

  // --- Production:
  virtual void produce(Event&);
  virtual void endJob();

private:
  mu2e::ProditionsHandle<mu2e::CaloDAQMap> _calodaqconds_h;

  void analyze_calorimeter_(mu2e::CaloDAQMap const& calodaqconds,
                            const mu2e::CalorimeterDataDecoder& cc,
                            std::unique_ptr<mu2e::CaloHitCollection> const& calo_hits,
                            std::unique_ptr<mu2e::CaloHitCollection> const& caphri_hits,
                            unsigned short& evtEnergy);

  void addPulse(uint16_t& crystalID, float& time, float& eDep,
                std::unique_ptr<mu2e::CaloHitCollection> const& hits_calo,
                std::unique_ptr<mu2e::CaloHitCollection> const& hits_caphri);

  int diagLevel_;

  art::InputTag caloFragmentsTag_;
  float digiSampling_;
  float deltaTPulses_, hitEDepMax_, hitEDepMin_, caphriEDepMax_, caphriEDepMin_, nPEperMeV_,
      noise2_, nSigmaNoise_;

  const int hexShiftPrint = 7;

  std::unordered_map<uint16_t, std::list<mu2e::CaloHit>>
      pulseMap_; // Temporary hack until the Calorimeter channel map is finalized
  mu2e::CaloDAQUtilities caloDAQUtil_;

  std::array<float, mu2e::CaloConst::_nCrystalChannel> peakADC2MeV_;
  std::array<float, mu2e::CaloConst::_nCrystalChannel> timeCalib_;

  long int total_events;
  long int total_hits;
  long int total_hits_good;
  long int total_hits_bad;
  std::map<mu2e::CaloDAQUtilities::CaloHitError, uint> failure_counter;
};

// ======================================================================

void art::CaloHitsFromDataDTCEvents::beginRun(art::Run& Run) {

  // FIX ME!
  // here we need to load the prodition-service with the calibrations
  // for converting teh peakADC into MeV and sync times
  for (size_t i = 0; i < peakADC2MeV_.size(); ++i) {
    peakADC2MeV_[i] = 0.0461333;
    timeCalib_[i] = 0.;
  }
}

void art::CaloHitsFromDataDTCEvents::addPulse(
    uint16_t& crystalID, float& time, float& eDep,
    std::unique_ptr<mu2e::CaloHitCollection> const& hits_calo,
    std::unique_ptr<mu2e::CaloHitCollection> const& hits_caphri) {

  bool addNewHit(true);
  bool isCaphri = mu2e::CrystalId(crystalID).isCaphri();
  size_t counter(0);
  for (auto& pulse : pulseMap_[crystalID]) {
    ++counter;
    if (std::fabs(pulse.time() - time) < deltaTPulses_) {

      float ratio = (eDep - pulse.energyDep()) / (eDep + pulse.energyDep());
      float eMean = (eDep + pulse.energyDep()) / 2.0;
      float sigmaR = 0.707 * sqrt(1.0 / eMean / nPEperMeV_ + noise2_ / eMean / eMean);

      if (abs(ratio) <= nSigmaNoise_ * sigmaR) {
        // combine the pulses
        pulse.setTime((pulse.time() + time) / 2.); // probably not necessary
        pulse.setEDep((pulse.energyDep() + eDep) / 2.);
        pulse.setNSiPMs(pulse.nSiPMs() + 1);
        addNewHit = false;
      } else if (eDep > pulse.energyDep()) {
        pulse.setTime(time); // probably not necessary
        pulse.setEDep(eDep);
        addNewHit = false;
      }

      // move the pulse in the final collection
      if (isCaphri) {
        hits_caphri->emplace_back(std::move(pulse));
      } else {
        hits_calo->emplace_back(std::move(pulse));
      }
      break;
    }
  }
  if (addNewHit) {
    pulseMap_[crystalID].emplace_back(mu2e::CaloHit(crystalID, 1, time, eDep));
  }
}

art::CaloHitsFromDataDTCEvents::CaloHitsFromDataDTCEvents(
    const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_(config().diagLevel()), digiSampling_(config().digiSampling()),
    deltaTPulses_(config().deltaTPulses()), hitEDepMax_(config().hitEDepMax()),
    hitEDepMin_(config().hitEDepMin()), caphriEDepMax_(config().caphriEDepMax()),
    caphriEDepMin_(config().caphriEDepMin()), nPEperMeV_(config().nPEperMeV()),
    noise2_(config().noiseLevelMeV() * config().noiseLevelMeV()),
    nSigmaNoise_(config().nSigmaNoise()), caloDAQUtil_("CaloHitsFromDataDTCEvents") {
  pulseMap_.reserve(4000);
  produces<mu2e::CaloHitCollection>("calo");
  produces<mu2e::CaloHitCollection>("caphri");
  produces<mu2e::IntensityInfoCalo>();
  total_events = 0;
  total_hits = 0;
  total_hits_good = 0;
  total_hits_bad = 0;
}

// ----------------------------------------------------------------------

void art::CaloHitsFromDataDTCEvents::produce(Event& event) {
  pulseMap_.clear();

  art::EventNumber_t eventNumber = event.event();

  mu2e::CaloDAQMap const& calodaqconds = _calodaqconds_h.get(event.id());

  // Collection of CaloHits for the event
  std::unique_ptr<mu2e::CaloHitCollection> calo_hits(new mu2e::CaloHitCollection);
  std::unique_ptr<mu2e::CaloHitCollection> caphri_hits(new mu2e::CaloHitCollection);

  // IntensityInfoCalo
  std::unique_ptr<mu2e::IntensityInfoCalo> int_info(new mu2e::IntensityInfoCalo);

  size_t totalSize = 0;
  size_t numCalDecoders = 0;
  unsigned short evtEnergy(0);

  std::unique_ptr<mu2e::CalorimeterDataDecoders> decoderColl(new mu2e::CalorimeterDataDecoders);
  artdaq::Fragments fragments = caloDAQUtil_.getFragments(event);
  for (const auto& frag : fragments) {
    mu2e::DTCEventFragment eventFragment(frag);
    auto caloSEvents =
        eventFragment.getSubsystemData(DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter);
    for (auto& subevent : caloSEvents) {
      decoderColl->emplace_back(subevent);
      auto& decoder = decoderColl->back();
      analyze_calorimeter_(calodaqconds, decoder, calo_hits, caphri_hits, evtEnergy);
      for (size_t i = 0; i < decoder.block_count(); ++i) {
        totalSize += decoder.blockSizeBytes(i);
      }
      numCalDecoders++;
    }
  }

  if (numCalDecoders == 0) {
    std::cout << "[CaloDigiFromDTCEvents::produce] found no Calorimeter decoders!" << std::endl;
    // Must put empty vectors anyway
    event.put(std::move(int_info));
    event.put(std::move(calo_hits), "calo");
    event.put(std::move(caphri_hits), "caphri");
    return;
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "[CaloHitsFromDTCEvents::produce] Run " << event.run() << ", subrun "
              << event.subRun() << ", event " << eventNumber << " has " << numCalDecoders
              << " CALO decoders." << std::endl;
    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  int_info->setCaloEnergy(evtEnergy);
  int_info->setNCaloHits(calo_hits->size());
  int_info->setNCaphriHits(caphri_hits->size());
  event.put(std::move(int_info));

  // Store the calo hits in the event
  event.put(std::move(calo_hits), "calo");
  event.put(std::move(caphri_hits), "caphri");

} // produce()

void art::CaloHitsFromDataDTCEvents::analyze_calorimeter_(
    mu2e::CaloDAQMap const& calodaqconds, const mu2e::CalorimeterDataDecoder& cc,
    std::unique_ptr<mu2e::CaloHitCollection> const& calo_hits,
    std::unique_ptr<mu2e::CaloHitCollection> const& caphri_hits, unsigned short& evtEnergy) {

  auto dtcID = cc.event_.GetDTCID();

  // Loop through the ROCs of this caloDecoder
  for (size_t iROC = 0; iROC < cc.block_count(); iROC++) {

    // only data_type 1 for now
    // if (data_type_ == 1)
    //

    auto hits = cc.GetCalorimeterHitTestForTrigger(iROC);
    if (hits == nullptr) {
      mf::LogError("CaloHitsFromDataDTCEvents") << "Error retrieving Calorimeter data from block "
                                                << iROC << "! Aborting processing of this block!";
      continue;
    }

    // Loop through the hits of this ROC
    total_hits += hits->size();
    for (size_t hitIdx = 0; hitIdx < hits->size(); hitIdx++) {

      mu2e::CalorimeterDataDecoder::CalorimeterHitTestDataPacket& thisHitPacket =
          hits->at(hitIdx).first;
      uint16_t thisHitPeak = hits->at(hitIdx).second;

      // Check if hit was decoded correctly
      auto errorCode = caloDAQUtil_.isHitGood(hits->at(hitIdx));
      if (errorCode) {
        failure_counter[errorCode]++;
        total_hits_bad++;
        if (diagLevel_ > 1) {
          std::cout << "[CaloDigisFromDataDecoders] BAD calo hit! DTC: " << dtcID
                    << ", ROC: " << iROC << ", hit number: " << hitIdx
                    << " [failure code: " << errorCode << "]" << std::endl;
          caloDAQUtil_.printCaloPulse(thisHitPacket);
        }
        continue;
      }
      total_hits_good++;

      if (diagLevel_ > 2) {
        std::cout << "[CaloHitsFromDataDTCEvents] calo hit " << hitIdx << std::endl;
        caloDAQUtil_.printCaloPulse(thisHitPacket);
      }

      // Fill the CaloHitCollection
      mu2e::CaloRawSiPMId rawId(thisHitPacket.BoardID, thisHitPacket.ChannelID);
      mu2e::CaloSiPMId offlineId = calodaqconds.offlineId(rawId);
      uint16_t crystalID = offlineId.crystal().id();
      uint16_t SiPMID = offlineId.id();

      size_t peakIndex = thisHitPacket.IndexOfMaxDigitizerSample;
      float eDep = thisHitPeak * peakADC2MeV_[SiPMID];
      float time = thisHitPacket.Time + peakIndex * digiSampling_ + timeCalib_[SiPMID];

      bool isCaphri = offlineId.crystal().isCaphri();

      // FIX ME! WE NEED TO CHECK IF TEH PULSE IS SATURATED HERE
      if (((eDep >= hitEDepMin_) || (isCaphri && (eDep >= caphriEDepMin_))) &&
          ((eDep < hitEDepMax_) || (isCaphri && (eDep < caphriEDepMax_)))) {
        addPulse(crystalID, time, eDep, calo_hits, caphri_hits);
        evtEnergy += eDep;
      }

    } // End loop over hits
  }   // End loop over ROCs
}

void art::CaloHitsFromDataDTCEvents::endJob() {

  if (diagLevel_ > 0) {
    std::cout << "\n ----- [CaloHitsFromDataDTCEvents] Decoding errors summary ----- " << std::endl;
    std::cout << "Total events: " << total_events << std::endl;
    std::cout << "Total hits: " << total_hits << std::endl;
    std::cout << "Total good hits: " << total_hits_good << std::endl;
    std::cout << "Total bad hits: " << total_hits_bad << std::endl;
    for (auto fail : failure_counter) {
      std::cout << "Failure mode " << fail.first << " [bad "
                << caloDAQUtil_.getCaloHitErrorName(fail.first) << "], count: " << fail.second
                << " (" << int(100. * fail.second / total_hits) << "%)\n";
    }
  }
}

// ======================================================================

DEFINE_ART_MODULE(art::CaloHitsFromDataDTCEvents)

// ======================================================================
