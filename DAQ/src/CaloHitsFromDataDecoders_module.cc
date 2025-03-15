// ======================================================================
//
// Make CaloHits directly from CaloDataDecoders
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include <artdaq-core/Data/Fragment.hh>

#include "Offline/DAQ/inc/CaloDAQUtilities.hh"

#include <iostream>

#include <string>

#include <array>
#include <list>
#include <memory>
#include <unordered_map>
#include <vector>

namespace art {
class CaloHitsFromDataDecoders;
}

// ======================================================================

class art::CaloHitsFromDataDecoders : public EDProducer {

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
    fhicl::Atom<art::InputTag> caloTag{fhicl::Name("caloTag"), fhicl::Comment("caloTag")};
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
  explicit CaloHitsFromDataDecoders(const art::EDProducer::Table<Config>& config);
  virtual ~CaloHitsFromDataDecoders() {}

  virtual void beginRun(art::Run&) override;

  // --- Production:
  virtual void produce(Event&);

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
      pulseMap_; // Temporary hack until the Calorimeter channel map is finialized
  mu2e::CaloDAQUtilities caloDAQUtil_;

  std::array<float, 674 * 4> peakADC2MeV_;
  std::array<float, 674 * 4> timeCalib_;

};

// ======================================================================

void art::CaloHitsFromDataDecoders::beginRun(art::Run& Run) {

  // FIX ME!
  // here we need to load the prodition-service with the calibrations
  // for converting teh peakADC into MeV. I decided to fill an array
  // with the conversion constants in order to speed up the access
  // NOW FILLING THE ARRAY WITH A DUMMY VALUE
  for (size_t i = 0; i < peakADC2MeV_.size(); ++i) {
    peakADC2MeV_[i] = 0.0461333;
    timeCalib_[i] = 0.;
  }

}

void art::CaloHitsFromDataDecoders::addPulse(
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

art::CaloHitsFromDataDecoders::CaloHitsFromDataDecoders(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, diagLevel_(config().diagLevel()),
    caloFragmentsTag_(config().caloTag()), digiSampling_(config().digiSampling()),
    deltaTPulses_(config().deltaTPulses()), hitEDepMax_(config().hitEDepMax()),
    hitEDepMin_(config().hitEDepMin()), caphriEDepMax_(config().caphriEDepMax()),
    caphriEDepMin_(config().caphriEDepMin()), nPEperMeV_(config().nPEperMeV()),
    noise2_(config().noiseLevelMeV() * config().noiseLevelMeV()),
    nSigmaNoise_(config().nSigmaNoise()), caloDAQUtil_("CaloHitsFromDataDecoders") {
  pulseMap_.reserve(4000);
  produces<mu2e::CaloHitCollection>("calo");
  produces<mu2e::CaloHitCollection>("caphri");
  produces<mu2e::IntensityInfoCalo>();
}

// ----------------------------------------------------------------------

void art::CaloHitsFromDataDecoders::produce(Event& event) {
  pulseMap_.clear();

  art::EventNumber_t eventNumber = event.event();

  mu2e::CaloDAQMap const& calodaqconds = _calodaqconds_h.get(event.id());

  // Collection of CaloHits for the event
  std::unique_ptr<mu2e::CaloHitCollection> calo_hits(new mu2e::CaloHitCollection);
  std::unique_ptr<mu2e::CaloHitCollection> caphri_hits(new mu2e::CaloHitCollection);

  // IntensityInfoCalo
  std::unique_ptr<mu2e::IntensityInfoCalo> int_info(new mu2e::IntensityInfoCalo);

  size_t totalSize(0), numCalFrags(0);
  unsigned short evtEnergy(0);

  auto fragmentHandle =
      event.getValidHandle<std::vector<mu2e::CalorimeterDataDecoder>>(caloFragmentsTag_);

  for (auto frag : *fragmentHandle) {
    analyze_calorimeter_(calodaqconds, frag, calo_hits, caphri_hits, evtEnergy);
    for (size_t i = 0; i < frag.block_count(); ++i) {
      totalSize += frag.blockSizeBytes(i);
    }
    numCalFrags++;
  }

  if (numCalFrags == 0) {
    std::cout << "[CaloHitsFromDataDecoders::produce] found no Calorimeter fragments!" << std::endl;
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::CaloHitsFromDataDecoders::produce exiting eventNumber=" << (int)(event.event())
              << " / timestamp=" << (int)eventNumber << std::endl;
  }

  int_info->setCaloEnergy(evtEnergy);
  int_info->setNCaloHits(calo_hits->size());
  event.put(std::move(int_info));

  // Store the calo hits in the event
  event.put(std::move(calo_hits), "calo");
  event.put(std::move(caphri_hits), "caphri");

} // produce()

void art::CaloHitsFromDataDecoders::analyze_calorimeter_(
    mu2e::CaloDAQMap const& calodaqconds,
    const mu2e::CalorimeterDataDecoder& cc, std::unique_ptr<mu2e::CaloHitCollection> const& calo_hits,
    std::unique_ptr<mu2e::CaloHitCollection> const& caphri_hits, unsigned short& evtEnergy) {

  if (diagLevel_ > 1) {
    caloDAQUtil_.printCaloFragmentInfo(cc);
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
      mf::LogError("CaloHitsFromDataDecoders")
          << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (diagLevel_ > 1) {
      caloDAQUtil_.printCaloFragmentHeader(hdr);
    }

    if (hdr->GetPacketCount() == 0)
      continue;

    auto calData = cc.GetCalorimeterHitData(curBlockIdx);
    if (calData == nullptr) {
      mf::LogError("CaloHitsFromDataDecoders")
          << "Error retrieving Calorimeter data from block " << curBlockIdx
          << "! Aborting processing of this block!";
      continue;
    }

    auto hits = cc.GetCalorimeterHitsForTrigger(curBlockIdx);
    bool err = false;
    for (size_t hitIdx = 0; hitIdx < calData->size(); hitIdx++) {

      // Fill the CaloDigiCollection
      if (hitIdx > hits.size()) {
        mf::LogError("CaloHitsFromDataDecoders")
            << "Error retrieving Calorimeter data from block " << curBlockIdx << " for hit "
            << hitIdx << "! Aborting processing of this block!";
        err = true;
        break;
      }

      if (diagLevel_ > 0) {
        std::cout << "[CaloHitsFromDataDecoders] calo hit " << hitIdx << std::endl;
        caloDAQUtil_.printCaloPulse(hits[hitIdx].first);
      }

      uint16_t packetid = hits[hitIdx].first.DIRACA;
      uint16_t dirac = packetid & 0xFF;
      uint16_t diracChannel = (packetid >>8) & 0x1F;
      mu2e::CaloRawSiPMId rawId(dirac,diracChannel);
      mu2e::CaloSiPMId offlineId = calodaqconds.offlineId(rawId);

      uint16_t crystalID = offlineId.crystal().id();
      uint16_t sipmID = offlineId.SiPMLocalId();

      size_t peakIndex = hits[hitIdx].first.IndexOfMaxDigitizerSample;
      // float  eDep(0);
      //      if (hits[hitIdx].first.IndexOfMaxDigitizerSample < hits[hitIdx].second.size()) {
      // eDep = hits[hitIdx].second.at(peakIndex) * peakADC2MeV_[sipmID];
      float eDep = hits[hitIdx].second * peakADC2MeV_[sipmID];
      //      }
      float time = hits[hitIdx].first.Time + peakIndex * digiSampling_ + timeCalib_[sipmID];

      bool  isCaphri = offlineId.crystal().isCaphri();
      // FIX ME! WE NEED TO CHECK IF TEH PULSE IS SATURATED HERE
      if (((eDep >= hitEDepMin_) || (isCaphri && (eDep >= caphriEDepMin_))) &&
          ((eDep < hitEDepMax_) || (isCaphri && (eDep < caphriEDepMax_)))) {
        addPulse(crystalID, time, eDep, calo_hits, caphri_hits);
        evtEnergy += eDep;
      }
      if (diagLevel_ > 1) {
        // Until we have the final mapping, the BoardID is just a placeholder
        // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);

        caloDAQUtil_.printAllHitInfo(crystalID, sipmID, hdr, hits[hitIdx].first,
                                     hits[hitIdx].second);
      } // End debug output

    } // End loop over readout channels in DataBlock

    if (err)
      continue;
  }
}
// ======================================================================

DEFINE_ART_MODULE(art::CaloHitsFromDataDecoders)

// ======================================================================
