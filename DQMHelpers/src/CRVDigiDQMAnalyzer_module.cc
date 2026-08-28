// Thin art analyzer that constructs, books, and fills CRVDigiDQM.
//
// Original Author: R. Mina

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/DQMHelpers/inc/CRVDigiDQM.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "TH1F.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

class CRVDigiDQMAnalyzer : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> crvDigiTag{
        Name("crvDigiTag"),
        Comment("CRV digi producer"),
        art::InputTag{"CrvDigi"}};
    fhicl::Atom<art::InputTag> crvStatusTag{
        Name("crvStatusTag"),
        Comment("CRV status producer"),
        art::InputTag{"CrvDigi"}};
    fhicl::Atom<std::string> outputTag{
        Name("outputTag"),
        Comment("TFileService subdirectory; empty books in the module directory"),
        ""};
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic level"), 0};

    fhicl::Atom<int> nBinsDigisPerEvt{
        Name("nBinsDigisPerEvt"), Comment("Bins for h1_digisPerEvt"), 200};
    fhicl::Atom<float> maxDigisPerEvt{
        Name("maxDigisPerEvt"), Comment("Upper edge for h1_digisPerEvt"), 4000.f};
    fhicl::Atom<int> nBinsPeakAdc{
        Name("nBinsPeakAdc"), Comment("Bins for h1_peakAdc"), 450};
    fhicl::Atom<float> maxPeakAdc{
        Name("maxPeakAdc"), Comment("Upper edge for h1_peakAdc"), 4500.f};
    fhicl::Atom<int> nBinsTdc{Name("nBinsTdc"), Comment("Bins for h1_tdc"), 400};
    fhicl::Atom<float> maxTdc{Name("maxTdc"), Comment("Upper edge for h1_tdc"), 40000.f};
    fhicl::Atom<double> cfFraction{
        Name("cfFraction"), Comment("Constant-fraction timing threshold"), 0.20};
    fhicl::Atom<float> dtBinSize{
        Name("dtBinSize"), Comment("CF dt histogram bin width [ns]"), 0.5f};
    fhicl::Atom<float> dtRange{
        Name("dtRange"), Comment("CF dt histogram +/- range [ns]"), 100.f};
    fhicl::Atom<float> dtVsFebBinSize{
        Name("dtVsFebBinSize"), Comment("dtVsFeb bin width [ns]"), 2.f};
    fhicl::Atom<float> dtVsFebRange{
        Name("dtVsFebRange"), Comment("dtVsFeb +/- range [ns]"), 500.f};
    fhicl::Atom<int> minAmplitude{
        Name("minAmplitude"), Comment("Minimum CF amplitude (peak-baseline)"), 10};
    fhicl::Atom<int> avgBlockSize{
        Name("avgBlockSize"), Comment("Events per g_digisAvgVsEwt point"), 30};
    fhicl::Atom<int> avgGraphPoints{
        Name("avgGraphPoints"), Comment("Max points in g_digisAvgVsEwt"), 1000};
    fhicl::Atom<int> channelsWindowEwts{
        Name("channelsWindowEwts"),
        Comment("EWT span for h1_channelsLastEwt"),
        50000};
    fhicl::Atom<bool> fillInclusive{
        Name("fillInclusive"),
        Comment("Also fill ValCrvDigi BarId/SiPM/ADC histograms"),
        true};
    fhicl::Atom<bool> fillCrvIdRates{
        Name("fillCrvIdRates"),
        Comment("Book CRVId rate maps and crvDigisPerChannel"),
        true};
    fhicl::Atom<bool> kppReadout{
        Name("kppReadout"),
        Comment("KPP FEB-axis sizing (ROC 1-2); ROC4->ROC2 is the unpacker's job"),
        true};
    fhicl::Atom<bool> fillSectorOccupancy{
        Name("fillSectorOccupancy"),
        Comment("Fill per-sector crvDigisPerChannelAndEvent_* using GeometryService"),
        false};
    fhicl::Atom<int> histDigisBins{
        Name("histDigisBins"),
        Comment("Bins for crvDigisPerChannelAndEvent_CRVsector*"),
        200};
    fhicl::Atom<double> histDigisStart{
        Name("histDigisStart"),
        Comment("Low edge for crvDigisPerChannelAndEvent_CRVsector*"),
        0.0};
    fhicl::Atom<double> histDigisEnd{
        Name("histDigisEnd"),
        Comment("High edge for crvDigisPerChannelAndEvent_CRVsector*"),
        0.1};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CRVDigiDQMAnalyzer(const Parameters& conf);

  void beginJob() override;
  void beginRun(const art::Run& run) override;
  void analyze(const art::Event& event) override;
  void endJob() override;

private:
  static CRVDigiDQM::Config makeHelperConfig(const Config& conf);

  art::InputTag crvDigiTag_;
  art::InputTag crvStatusTag_;
  std::string outputTag_;
  int diagLevel_;
  bool fillSectorOccupancy_;
  int histDigisBins_;
  double histDigisStart_;
  double histDigisEnd_;
  CRVDigiDQM dqm_;
  bool sectorMapSent_{false};
};

CRVDigiDQM::Config CRVDigiDQMAnalyzer::makeHelperConfig(const Config& conf)
{
  CRVDigiDQM::Config c;
  c.nBinsDigisPerEvt = conf.nBinsDigisPerEvt();
  c.maxDigisPerEvt = conf.maxDigisPerEvt();
  c.nBinsPeakAdc = conf.nBinsPeakAdc();
  c.maxPeakAdc = conf.maxPeakAdc();
  c.nBinsTdc = conf.nBinsTdc();
  c.maxTdc = conf.maxTdc();
  c.cfFraction = conf.cfFraction();
  c.dtBinSize = conf.dtBinSize();
  c.dtRange = conf.dtRange();
  c.dtVsFebBinSize = conf.dtVsFebBinSize();
  c.dtVsFebRange = conf.dtVsFebRange();
  c.minAmplitude = conf.minAmplitude();
  c.avgBlockSize = static_cast<std::size_t>(std::max(conf.avgBlockSize(), 1));
  c.avgGraphPoints = static_cast<std::size_t>(std::max(conf.avgGraphPoints(), 1));
  c.channelsWindowEwts =
      static_cast<std::size_t>(std::max(conf.channelsWindowEwts(), 1));
  c.fillInclusive = conf.fillInclusive();
  c.fillCrvIdRates = conf.fillCrvIdRates();
  c.kppReadout = conf.kppReadout();
  return c;
}

CRVDigiDQMAnalyzer::CRVDigiDQMAnalyzer(const Parameters& conf) :
    art::EDAnalyzer{conf},
    crvDigiTag_(conf().crvDigiTag()),
    crvStatusTag_(conf().crvStatusTag()),
    outputTag_(conf().outputTag()),
    diagLevel_(conf().diagLevel()),
    fillSectorOccupancy_(conf().fillSectorOccupancy()),
    histDigisBins_(std::max(conf().histDigisBins(), 1)),
    histDigisStart_(conf().histDigisStart()),
    histDigisEnd_(conf().histDigisEnd()),
    dqm_(makeHelperConfig(conf()))
{}

void CRVDigiDQMAnalyzer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  if (outputTag_.empty()) {
    // TFileService already gives each module its own directory; book into it.
    art::TFileDirectory dir = *tfs;
    dqm_.Book(dir);
  } else {
    dqm_.Book(tfs->mkdir(outputTag_));
  }
}

void CRVDigiDQMAnalyzer::beginRun(const art::Run&)
{
  if (!fillSectorOccupancy_ || sectorMapSent_) {
    return;
  }
  sectorMapSent_ = true;

  GeomHandle<CosmicRayShield> CRS;
  auto const& crvSectors = CRS->getCRSScintillatorShields();
  std::vector<std::string> sectorNames;
  sectorNames.reserve(crvSectors.size());
  for (std::size_t i = 0; i < crvSectors.size(); ++i) {
    sectorNames.emplace_back(crvSectors.at(i).name(""));
  }

  auto const& crvCounters = CRS->getAllCRSScintillatorBars();
  std::vector<int> channelToSector(crvCounters.size() * CRVId::nChanPerBar, -1);
  for (std::size_t channel = 0; channel < channelToSector.size(); ++channel) {
    CRSScintillatorBarIndex barIndex(
        static_cast<int>(channel / CRVId::nChanPerBar));
    channelToSector[channel] = CRS->getBar(barIndex).id().getShieldNumber();
  }

  dqm_.BookSectorOccupancy(sectorNames, channelToSector,
                           histDigisBins_, histDigisStart_, histDigisEnd_);
}

void CRVDigiDQMAnalyzer::analyze(const art::Event& event)
{
  art::Handle<CrvDigiCollection> digiHandle;
  event.getByLabel(crvDigiTag_, digiHandle);
  if (!digiHandle.isValid()) {
    if (diagLevel_ > 1) {
      std::cout << "[CRVDigiDQMAnalyzer] No CrvDigiCollection at "
                << crvDigiTag_ << std::endl;
    }
    return;
  }

  art::Handle<CrvStatusCollection> statusHandle;
  event.getByLabel(crvStatusTag_, statusHandle);
  const CrvStatusCollection emptyStatus;
  const CrvStatusCollection& status =
      (statusHandle.isValid() && statusHandle.product() != nullptr) ?
          *statusHandle :
          emptyStatus;

  dqm_.Fill(*digiHandle, status);
}

void CRVDigiDQMAnalyzer::endJob()
{
  dqm_.WriteGraphs();

  if (diagLevel_ > 0) {
    std::cout << "[CRVDigiDQMAnalyzer] Total events: " << dqm_.nEvents()
              << std::endl;
    std::cout << "[CRVDigiDQMAnalyzer] Total digis: " << dqm_.nDigis()
              << std::endl;
    std::cout << "[CRVDigiDQMAnalyzer] Active FEBs: " << dqm_.activeFEBs().size()
              << std::endl;
    for (const auto& [roc, febs] : dqm_.rocFEBMap()) {
      std::cout << "[CRVDigiDQMAnalyzer] ROC " << static_cast<int>(roc) << " has "
                << febs.size() << " FEBs:";
      for (auto feb : febs) {
        std::cout << " " << static_cast<int>(feb);
      }
      std::cout << std::endl;
    }
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CRVDigiDQMAnalyzer)
