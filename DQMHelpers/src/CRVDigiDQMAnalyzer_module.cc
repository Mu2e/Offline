//
// Thin art analyzer that constructs, books, and fills CRVDigiDQM.
//
// Original Author: R. Mina
//

#include "Offline/DQMHelpers/inc/CRVDigiDQM.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include <algorithm>
#include <iostream>
#include <string>

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
        Name("outputTag"), Comment("TFileService subdirectory"), "CRVDigiDQM"};
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
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CRVDigiDQMAnalyzer(const Parameters& conf);

  void beginJob() override;
  void analyze(const art::Event& event) override;
  void endJob() override;

private:
  static CRVDigiDQM::Config makeHelperConfig(const Config& conf);

  art::InputTag crvDigiTag_;
  art::InputTag crvStatusTag_;
  std::string outputTag_;
  int diagLevel_;
  CRVDigiDQM dqm_;
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
  c.minAmplitude = conf.minAmplitude();
  c.avgBlockSize = static_cast<std::size_t>(std::max(conf.avgBlockSize(), 1));
  c.avgGraphPoints = static_cast<std::size_t>(std::max(conf.avgGraphPoints(), 1));
  c.channelsWindowEwts =
      static_cast<std::size_t>(std::max(conf.channelsWindowEwts(), 1));
  c.fillInclusive = conf.fillInclusive();
  return c;
}

CRVDigiDQMAnalyzer::CRVDigiDQMAnalyzer(const Parameters& conf) :
    art::EDAnalyzer{conf},
    crvDigiTag_(conf().crvDigiTag()),
    crvStatusTag_(conf().crvStatusTag()),
    outputTag_(conf().outputTag()),
    diagLevel_(conf().diagLevel()),
    dqm_(makeHelperConfig(conf()))
{}

void CRVDigiDQMAnalyzer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  dqm_.Book(tfs->mkdir(outputTag_));
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
