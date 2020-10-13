// ======================================================================
//
// StrawAndCaloDigisFromFragments_plugin:  Add tracker/cal data products to the event
//
// ======================================================================

// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "mu2e-artdaq-core/Overlays/CalorimeterFragment.hh"
#include "mu2e-artdaq-core/Overlays/FragmentType.hh"
#include "mu2e-artdaq-core/Overlays/TrackerFragment.hh"

#include <artdaq-core/Data/Fragment.hh>

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/TrkTypes.hh"

#include <iostream>
#include <memory>
#include <string>

// ======================================================================
namespace mu2e {

class FragmentAna : public art::EDAnalyzer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    fhicl::Atom<int> parseCAL{fhicl::Name("parseCAL"), fhicl::Comment("parseCAL")};
    fhicl::Atom<int> parseTRK{fhicl::Name("parseTRK"), fhicl::Comment("parseTRK")};
    fhicl::Atom<art::InputTag> caloTag{fhicl::Name("caloTag"), fhicl::Comment("caloTag")};
    fhicl::Atom<art::InputTag> trkTag{fhicl::Name("trkTag"), fhicl::Comment("trkTag")};
  };

  explicit FragmentAna(const art::EDAnalyzer::Table<Config>& config);
  virtual ~FragmentAna() {}

  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void analyze(const art::Event& e) override;

private:
  void analyze_tracker_(const artdaq::Fragment& f);
  void analyze_calorimeter_(const artdaq::Fragment& f);

  int diagLevel_;

  int parseCAL_;
  int parseTRK_;

  art::InputTag trkFragmentsTag_;
  art::InputTag caloFragmentsTag_;

  TH1F* _hTrkNFragment;
  TH1F* _hTrkStrawId;
  TH1F* _hTrkTDC[4];
  TH1F* _hTrkTOT;
  TH1F* _hTrkPMP;
  TH1F* _hTrkMeanADC;
  TH1F* _hTrkMaxADC;
  TH1F* _hTrkWfSize;

  TH1F* _hCalNFragment;
  TH1F* _hCalROId;
  TH1F* _hCalT0;
  TH1F* _hCalPeakPos;
  TH1F* _hCalWfSize;
};

// ======================================================================

FragmentAna::FragmentAna(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}, diagLevel_(config().diagLevel()), parseCAL_(config().parseCAL()),
    parseTRK_(config().parseTRK()), trkFragmentsTag_(config().trkTag()),
    caloFragmentsTag_(config().caloTag()) {}

//--------------------------------------------------------------------------------
// create the histograms
//--------------------------------------------------------------------------------
void FragmentAna::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  art::TFileDirectory calDir = tfs->mkdir("calorimeter");
  art::TFileDirectory trkDir = tfs->mkdir("tracker");

  _hTrkNFragment = trkDir.make<TH1F>("hTrkNFragment", "n fragments from the tracker; nTrkFragments",
                                     1000, 0., 10000.);
  _hTrkStrawId =
      trkDir.make<TH1F>("hTrkStrawId", "trk fragment strawId; strawId", 20000, 0., 20000.);
  _hTrkTDC[0] = trkDir.make<TH1F>("hTrkTDC0", "trk fragment TDC0; TDC[0]", 264, 0., 264000);
  _hTrkTDC[1] = trkDir.make<TH1F>("hTrkTDC1", "trk fragment TDC1; TDC[1]", 264, 0., 264000);
  _hTrkTDC[2] = trkDir.make<TH1F>("hTrkTDCMean", "trk fragment average TDC; (TDC[0]+TDC[1])/2",
                                  2000, 0., 20000.);
  _hTrkTDC[3] = trkDir.make<TH1F>("hTrkTDCDelta", "trk fragment delta TDC; TDC[1]-TDC[0]", 220,
                                  -100., 10000.);
  _hTrkTOT =
      trkDir.make<TH1F>("hTrkTOT", "trk fragment average TOT; (TOT[0]+TOT[1])/2", 100, 0., 200.);
  _hTrkPMP = trkDir.make<TH1F>("hTrkPOP", "trk fragment average PMP; PMP", 100, 0., 200.);
  _hTrkMeanADC = trkDir.make<TH1F>("hTrkMeanADC", "trk fragment Mean ADC; <ADC>", 250, 0., 2500.);
  _hTrkMaxADC = trkDir.make<TH1F>("hTrkMaxADC", "trk fragment Max ADC; Max_ADC", 250, 0., 2500.);
  _hTrkWfSize = trkDir.make<TH1F>("hTrkWfSize", "trk fragment waveform size; trkFragment_wf_size",
                                  20, 0., 20.);

  _hCalNFragment = calDir.make<TH1F>(
      "hCalNFragment", "n fragments from the calorimeter; nCalFragments", 400, 0., 400.);
  _hCalROId =
      calDir.make<TH1F>("hCalROId", "calo fragment roId; calFragment_roId", 4000, 0., 4000.);
  _hCalT0 = calDir.make<TH1F>("hCalT0", "calo fragment t0; calFragment_t0 [ns]", 200, 0., 2000.);
  _hCalPeakPos =
      calDir.make<TH1F>("hCalPeakPos", "calo fragment peakPos; calFragment_peakPos", 100, 0., 100.);
  _hCalWfSize = calDir.make<TH1F>("hCalWfSize", "calo fragment waveform size; calFragment_wf_size",
                                  100, 0., 100.);
}

void FragmentAna::endJob() {}

//--------------------------------------------------------------------------------
void FragmentAna::analyze(const art::Event& event) {
  art::EventNumber_t eventNumber = event.event();

  art::Handle<artdaq::Fragments> trkFragments, calFragments;
  size_t numTrkFrags(0), numCalFrags(0);
  if (parseTRK_) {
    event.getByLabel(trkFragmentsTag_, trkFragments);
    if (!trkFragments.isValid()) {
      return;
    }
    numTrkFrags = trkFragments->size();
    _hTrkNFragment->Fill(numTrkFrags);
  }
  if (parseCAL_) {
    event.getByLabel(caloFragmentsTag_, calFragments);
    if (!calFragments.isValid()) {
      return;
    }
    numCalFrags = calFragments->size();
    _hCalNFragment->Fill(numCalFrags);
  }

  if (diagLevel_ > 1) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
              << ", event " << eventNumber << " has " << std::endl;
    std::cout << numTrkFrags << " TRK fragments, and ";
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    size_t totalSize = 0;
    for (size_t idx = 0; idx < numTrkFrags; ++idx) {
      auto size = ((*trkFragments)[idx]).size(); // * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tTRK Fragment " << idx << " has size " << size << std::endl;
    }
    for (size_t idx = 0; idx < numCalFrags; ++idx) {
      auto size = ((*calFragments)[idx]).size(); // * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
    }
    totalSize *= sizeof(artdaq::RawDataType);

    std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;
  }
  std::string curMode = "TRK";

  // Loop over the TRK and CAL fragments
  for (size_t idx = 0; idx < numTrkFrags + numCalFrags; ++idx) {

    auto curHandle = trkFragments;
    size_t curIdx = idx;
    if (idx >= numTrkFrags) {
      curIdx = idx - numTrkFrags;
      curHandle = calFragments;
    }
    const auto& fragment((*curHandle)[curIdx]);

    if (idx < numTrkFrags) {
      analyze_tracker_(fragment);
    } else {
      analyze_calorimeter_(fragment);
    }
  }
  if (diagLevel_ > 0) {
    std::cout << "mu2e::FragmentAna::produce exiting eventNumber=" << (int)(event.event())
              << " / timestamp=" << (int)eventNumber << std::endl;
  }
}

void FragmentAna::analyze_tracker_(const artdaq::Fragment& f) {
  mu2e::TrackerFragment cc(f);

  if (diagLevel_ > 1) {
    std::cout << std::endl;
    std::cout << "ArtFragment: ";
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

    auto block = cc.dataAtBlockIndex(curBlockIdx);
    if (block == nullptr) {
      mf::LogError("FragmentAna") << "Unable to retrieve header from block " << curBlockIdx << "!"
                                  << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    // Parse phyiscs information from TRK packets
    if (hdr.GetPacketCount() > 0 && parseTRK_ > 0) {

      // Create the StrawDigi data products
      auto trkDataVec = cc.GetTrackerData(curBlockIdx);
      if (trkDataVec.empty()) {
        mf::LogError("FragmentAna") << "Error retrieving Tracker data from DataBlock "
                                    << curBlockIdx << "! Aborting processing of this block!";
        continue;
      }

      for (auto& trkDataPair : trkDataVec) {
        mu2e::StrawId sid(trkDataPair.first->StrawIndex);
        mu2e::TrkTypes::TDCValues tdc = {trkDataPair.first->TDC0(), trkDataPair.first->TDC1()};
        mu2e::TrkTypes::TOTValues tot = {trkDataPair.first->TOT0, trkDataPair.first->TOT1};
        mu2e::TrkTypes::ADCValue pmp = trkDataPair.first->PMP;
        int sum{0};
        unsigned short maxadc{0};
        for (auto adc : trkDataPair.second) {
          sum += adc;
          maxadc = std::max(maxadc, adc);
        }
        // Fill the StrawDigiCollection
        _hTrkStrawId->Fill(sid.asUint16());
        _hTrkTDC[0]->Fill(tdc[0]);
        _hTrkTDC[1]->Fill(tdc[1]);
        _hTrkTDC[2]->Fill((tdc[0] + tdc[1]) / 2.);
        _hTrkTDC[3]->Fill(tdc[1] - tdc[0]);
        _hTrkTOT->Fill((tot[0] + tot[1]) / 2.);
        _hTrkPMP->Fill(pmp);
        int mean = (trkDataPair.second.size() != 0) ? sum / trkDataPair.second.size() : -1.;
        _hTrkMeanADC->Fill(mean);
        _hTrkMaxADC->Fill(maxadc);
        _hTrkWfSize->Fill(trkDataPair.second.size());
      }
    }
  }

  cc.ClearUpgradedPackets();
}

void FragmentAna::analyze_calorimeter_(const artdaq::Fragment& f) {
  mu2e::CalorimeterFragment cc(f);

  if (diagLevel_ > 1) {
    std::cout << std::endl;
    std::cout << "ArtFragment: ";
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

    auto block = cc.dataAtBlockIndex(curBlockIdx);
    if (block == nullptr) {
      mf::LogError("FragmentAna") << "Unable to retrieve header from block " << curBlockIdx << "!"
                                  << std::endl;
      continue;
    }
    auto hdr = block->GetHeader();

    if (hdr.GetPacketCount() > 0 &&
        parseCAL_ > 0) { // Parse phyiscs information from CAL packets

      auto calData = cc.GetCalorimeterData(curBlockIdx);
      if (calData == nullptr) {
        mf::LogError("FragmentAna") << "Error retrieving Calorimeter data from block "
                                    << curBlockIdx << "! Aborting processing of this block!";
        continue;
      }

      if (diagLevel_ > 0) {
        std::cout << "[StrawAndCaloDigiFromFragments] NEW CALDATA: NumberOfHits "
                  << calData->NumberOfHits << std::endl;
      }

      auto hits = cc.GetCalorimeterHits(curBlockIdx);

      bool err = false;
      for (auto& hit : hits) {

        // Fill the CaloDigiCollection
        
        // IMPORTANT NOTE: we don't have a final
        // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
        // placeholder. Because we still need to know which crystal a hit belongs to, we are
        // temporarily storing the 4-bit apdID and 12-bit crystalID in the Reserved DIRAC A slot.
        // Also, note that until we have an actual map, channel index does not actually correspond
        // to the physical readout channel on a ROC.
        uint16_t crystalID = hit.first.DIRACB & 0x0FFF;
        uint16_t roId = hit.first.DIRACB >> 12;
        _hCalROId->Fill(crystalID * 2 + roId);
        _hCalT0->Fill(hit.first.Time);
        _hCalPeakPos->Fill(hit.first.IndexOfMaxDigitizerSample);
        _hCalWfSize->Fill(hit.second.size());

      } // End loop over readout channels in DataBlock
      if (err)
        continue;
    }
  }
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::FragmentAna);
