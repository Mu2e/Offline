// ======================================================================
//
// CaloHitsFromFragments_plugin:  Add cal data products to the event
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

#include "RecoDataProducts/inc/CaloHit.hh"

#include <artdaq-core/Data/Fragment.hh>

#include "DAQ/inc/CaloDAQUtilities.hh"

#include <iostream>

#include <string>

#include <list>
#include <memory>
#include <unordered_map>

namespace art {
class CaloHitsFromFragments;
}

// ======================================================================

class art::CaloHitsFromFragments : public EDProducer {

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
    fhicl::Atom<float> pulseRatioMax{
        fhicl::Name("pulseRatioMax"),
        fhicl::Comment("Max of the ratio between two SiPM pulses coupled to the same crystal "
                       "within the time-gate")};
    fhicl::Atom<float> pulseRatioMin{
        fhicl::Name("pulseRatioMin"),
        fhicl::Comment("Min of the ratio between two SiPM pulses coupled to the same crystal "
                       "within the time-gate")};
  };

  // --- C'tor/d'tor:
  explicit CaloHitsFromFragments(const art::EDProducer::Table<Config>& config);
  virtual ~CaloHitsFromFragments() {}

  virtual void beginRun(art::Run&) override;

  // --- Production:
  virtual void produce(Event&);

private:
  void analyze_calorimeter_(const artdaq::Fragment& f,
                            std::unique_ptr<mu2e::CaloHitCollection> const& calo_hits,
                            std::unique_ptr<mu2e::CaloHitCollection> const& caphri_hits);

  void addPulse(uint16_t& crystalID, float& time, float& eDep);

  int diagLevel_;

  art::InputTag caloFragmentsTag_;
  float digiSampling_;
  float deltaTPulses_, pulseRatioMax_, pulseRatioMin_;

  const int hexShiftPrint = 7;

  std::unordered_map<uint16_t, std::list<CrystalInfo>>
      pulseMap_; // Temporary hack until the Calorimeter channel map is finialized
  mu2e::CaloDAQUtilities caloDAQUtil_;

  std::array<float, 674 * 4> peakADC2MeV_;
  std::array<int, 4> caphriCrystalID_;
};

// ======================================================================

void art::CaloHitsFromFragments::beginRun(art::Run& Run) {

  // FIX ME!
  // here we need to load the prodition-service with the calibrations
  // for converting teh peakADC into MeV. I decided to fill an array
  // with the conversion constants in order to speed up the access
  // NOW FILLING THE ARRAY WITH A DUMMY VALUE
  for (size_t i = 0; i < peakADC2MeV_.size(); ++i) {
    peakADC2MeV_[i] = 0.0461333;
  }

  // FIX ME!
  // the list of the SiPM-IDs of the channels that are used for the
  // calorimeter-lumi monitor should come from a pro-dition
  caphriCrystalID_ = {623, 624, 595, 596};
}

void art::CaloHitsFromFragments::addPulse(uint16_t& crystalID, float& time, float& eDep) {

  bool addNewHit(true);
  for (auto& pulse : pulseMap_[crystalID]) {
    // search if there is a hit matching in time and eDep
    if ((std::fabs(pulse._time - time) < deltaTPulses_) && (eDep / pulse._eDep <= pulseRatioMax_) &&
        (eDep / pulse._eDep >= pulseRatioMin_)) {
      // combine the pulses
      pulse._time = (pulse._time + time) / 2.;
      pulse._eDep = (pulse._eDep + eDep) / 2.;
      pulse._nSiPM += 1;
      addNewHit = false;
      break;
    }
  }
  if (addNewHit) {
    pulseMap_[crystalID].push_back(CrystalInfo(crystalID, 1, time, eDep));
  }
}

art::CaloHitsFromFragments::CaloHitsFromFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}, diagLevel_(config().diagLevel()),
    caloFragmentsTag_(config().caloTag()), digiSampling_(config().digiSampling()),
    deltaTPulses_(config().deltaTPulses()), pulseRatioMax_(config().pulseRatioMax()),
    pulseRatioMin_(config().pulseRatioMin()), caloDAQUtil_("CaloHitsFromFragments") {
  produces<mu2e::CaloHitCollection>();
  produces<mu2e::CaloHitCollection>("caphri");
}

// ----------------------------------------------------------------------

void art::CaloHitsFromFragments::produce(Event& event) {
  pulseMap_.clear();
  art::EventNumber_t eventNumber = event.event();

  // Collection of CaloHits for the event
  std::unique_ptr<mu2e::CaloHitCollection> calo_hits(new mu2e::CaloHitCollection);
  std::unique_ptr<mu2e::CaloHitCollection> caphri_hits(new mu2e::CaloHitCollection);

  art::Handle<artdaq::Fragments> calFragments;
  size_t numCalFrags(0);
  size_t totalSize = 0;
  event.getByLabel(caloFragmentsTag_, calFragments);
  if (!calFragments.isValid()) {
    std::cout << "[CaloHitsFromFragments::produce] found no Calorimeter fragments!" << std::endl;
    event.put(std::move(calo_hits));
    event.put(std::move(caphri_hits));
    return;
  }
  numCalFrags = calFragments->size();
  for (size_t idx = 0; idx < numCalFrags; ++idx) {
    auto size = ((*calFragments)[idx]).sizeBytes(); // * sizeof(artdaq::RawDataType);
    totalSize += size;
    analyze_calorimeter_((*calFragments)[idx], calo_hits, caphri_hits);
    //      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    std::cout << "Total Size: " << (int)totalSize << " bytes." << std::endl;
  }

  if (diagLevel_ > 0) {
    std::cout << "mu2e::CaloHitsFromFragments::produce exiting eventNumber=" << (int)(event.event())
              << " / timestamp=" << (int)eventNumber << std::endl;
  }

  // Store the calo hits in the event
  event.put(std::move(calo_hits));
  event.put(std::move(caphri_hits), "caphri");

} // produce()

void art::CaloHitsFromFragments::analyze_calorimeter_(
    const artdaq::Fragment& f, std::unique_ptr<mu2e::CaloHitCollection> const& calo_hits,
    std::unique_ptr<mu2e::CaloHitCollection> const& caphri_hits) {
  mu2e::CalorimeterFragment cc(f);

  if (diagLevel_ > 1) {
    caloDAQUtil_.printCaloFragmentInfo(f, cc);
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
      mf::LogError("CaloHitsFromFragments")
          << "Unable to retrieve block " << curBlockIdx << "!" << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (diagLevel_ > 1) {
      caloDAQUtil_.printCaloFragmentHeader(hdr);
    }

    if (hdr.GetPacketCount() == 0)
      continue;

    auto calData = cc.GetCalorimeterData(curBlockIdx);
    if (calData == nullptr) {
      mf::LogError("CaloHitsFromFragments")
          << "Error retrieving Calorimeter data from block " << curBlockIdx
          << "! Aborting processing of this block!";
      continue;
    }

    if (diagLevel_ > 0) {
      std::cout << "[CaloHitsFromFragments] NEW CALDATA: NumberOfHits " << calData->NumberOfHits
                << std::endl;
    }

    auto hits = cc.GetCalorimeterHits(curBlockIdx);
    bool err = false;
    for (size_t hitIdx = 0; hitIdx < calData->NumberOfHits; hitIdx++) {

      // Fill the CaloDigiCollection
      if (hitIdx > hits.size()) {
        mf::LogError("CaloHitsFromFragments")
            << "Error retrieving Calorimeter data from block " << curBlockIdx << " for hit "
            << hitIdx << "! Aborting processing of this block!";
        err = true;
        break;
      }

      if (diagLevel_ > 0) {
        std::cout << "[CaloHitsFromFragments] calo hit " << hitIdx << std::endl;
        caloDAQUtil_.printCaloPulse(hits[hitIdx].first);
      }

      // IMPORTANT NOTE: we don't have a final
      // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
      // placeholder. Because we still need to know which crystal a hit belongs to, we are
      // temporarily storing the 4-bit sipmID and 12-bit crystalID in the Reserved DIRAC A slot.
      // Also, note that until we have an actual map, channel index does not actually correspond
      // to the physical readout channel on a ROC.
      uint16_t crystalID =
          caloDAQUtil_.getCrystalID(hits[hitIdx].first); // hits[hitIdx].first.DIRACB & 0x0FFF;
      uint16_t sipmID =
          caloDAQUtil_.getSiPMID(hits[hitIdx].first); // hits[hitIdx].first.DIRACB >> 12;

      size_t peakIndex = hits[hitIdx].first.IndexOfMaxDigitizerSample;
      float eDep(0);
      if (hits[hitIdx].first.IndexOfMaxDigitizerSample < hits[hitIdx].second.size()) {
        eDep = hits[hitIdx].second.at(peakIndex) * peakADC2MeV_[sipmID];
      }
      float time = hits[hitIdx].first.Time + peakIndex * digiSampling_;

      addPulse(crystalID, time, eDep);

      if (diagLevel_ > 1) {
        // Until we have the final mapping, the BoardID is just a placeholder
        // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);

        caloDAQUtil_.printAllHitInfo(crystalID, sipmID, hdr, hits[hitIdx].first,
                                     hits[hitIdx].second);
      } // End debug output

    } // End loop over readout channels in DataBlock

    // now create the CaloHitCollection
    for (auto& crystal : pulseMap_) {
      bool isCaphri = std::find(caphriCrystalID_.begin(), caphriCrystalID_.end(), crystal.first) !=
                      caphriCrystalID_.end();
      for (auto& crystalInfo : crystal.second) {
        if (isCaphri) {
          caphri_hits->emplace_back(crystal.first, crystalInfo._nSiPM, crystalInfo._time,
                                    crystalInfo._eDep);
        } else {
          calo_hits->emplace_back(crystal.first, crystalInfo._nSiPM, crystalInfo._time,
                                  crystalInfo._eDep);
        }
      }
    }

    if (err)
      continue;
  }
}
// ======================================================================

DEFINE_ART_MODULE(art::CaloHitsFromFragments)

// ======================================================================
