// Thin art analyzer that constructs, books, and fills CRVStatusDQM.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/CRVStatusDQM.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include <algorithm>
#include <iostream>
#include <string>

namespace mu2e {

class CRVStatusDQMAnalyzer : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> crvStatusTag{
        Name("crvStatusTag"),
        Comment("CRV status producer"),
        art::InputTag{"CrvDigi"}};
    fhicl::Atom<art::InputTag> crvDaqErrorTag{
        Name("crvDaqErrorTag"),
        Comment("CRV DAQ-error producer"),
        art::InputTag{"CrvDigi"}};
    fhicl::Atom<std::string> outputTag{
        Name("outputTag"),
        Comment("TFileService subdirectory; empty books in the module directory"),
        ""};
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic level"), 0};

    fhicl::Atom<int> nBinsLatency{
        Name("nBinsLatency"), Comment("Bins for linkLatency"), 1024};
    fhicl::Atom<float> maxLinkLatency{
        Name("maxLinkLatency"), Comment("Upper edge for linkLatency"), 4096.f};
    fhicl::Atom<int> nBinsTriggerCount{
        Name("nBinsTriggerCount"), Comment("Bins for triggerCount"), 256};
    fhicl::Atom<float> maxTriggerCount{
        Name("maxTriggerCount"), Comment("Upper edge for triggerCount"), 65535.f};
    fhicl::Atom<int> nBinsWordCount{
        Name("nBinsWordCount"), Comment("Bins for wordCount"), 256};
    fhicl::Atom<float> maxWordCount{
        Name("maxWordCount"), Comment("Upper edge for wordCount"), 65535.f};
    fhicl::Atom<int> nBinsEwtMismatch{
        Name("nBinsEwtMismatch"), Comment("Bins for ewtMismatch"), 201};
    fhicl::Atom<float> maxEwtMismatch{
        Name("maxEwtMismatch"), Comment("Abs range for ewtMismatch"), 100.f};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CRVStatusDQMAnalyzer(const Parameters& conf);

  void beginJob() override;
  void analyze(const art::Event& event) override;
  void endSubRun(const art::SubRun& subrun) override;
  void endJob() override;

private:
  static CRVStatusDQM::Config makeHelperConfig(const Config& conf);

  art::InputTag crvStatusTag_;
  art::InputTag crvDaqErrorTag_;
  std::string outputTag_;
  int diagLevel_;
  CRVStatusDQM dqm_;
};

CRVStatusDQM::Config CRVStatusDQMAnalyzer::makeHelperConfig(const Config& conf)
{
  CRVStatusDQM::Config c;
  c.nBinsLatency = std::max(conf.nBinsLatency(), 1);
  c.maxLinkLatency = conf.maxLinkLatency();
  c.nBinsTriggerCount = std::max(conf.nBinsTriggerCount(), 1);
  c.maxTriggerCount = conf.maxTriggerCount();
  c.nBinsWordCount = std::max(conf.nBinsWordCount(), 1);
  c.maxWordCount = conf.maxWordCount();
  c.nBinsEwtMismatch = std::max(conf.nBinsEwtMismatch(), 1);
  c.maxEwtMismatch = conf.maxEwtMismatch();
  return c;
}

CRVStatusDQMAnalyzer::CRVStatusDQMAnalyzer(const Parameters& conf) :
    art::EDAnalyzer{conf},
    crvStatusTag_(conf().crvStatusTag()),
    crvDaqErrorTag_(conf().crvDaqErrorTag()),
    outputTag_(conf().outputTag()),
    diagLevel_(conf().diagLevel()),
    dqm_(makeHelperConfig(conf()))
{}

void CRVStatusDQMAnalyzer::beginJob()
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

void CRVStatusDQMAnalyzer::analyze(const art::Event& event)
{
  art::Handle<CrvStatusCollection> statusHandle;
  event.getByLabel(crvStatusTag_, statusHandle);
  const CrvStatusCollection emptyStatus;
  const CrvStatusCollection& status =
      (statusHandle.isValid() && statusHandle.product() != nullptr) ?
          *statusHandle :
          emptyStatus;

  art::Handle<CrvDAQerrorCollection> daqHandle;
  event.getByLabel(crvDaqErrorTag_, daqHandle);
  if (daqHandle.isValid() && daqHandle.product() != nullptr) {
    dqm_.Fill(status, *daqHandle);
  } else {
    dqm_.Fill(status);
  }

  if (diagLevel_ > 1) {
    std::cout << "[CRVStatusDQMAnalyzer] " << event.id()
              << " nStatus=" << status.size()
              << " nRocSnap=" << dqm_.lastEventRocs().size() << std::endl;
  }
}

void CRVStatusDQMAnalyzer::endSubRun(const art::SubRun& subrun)
{
  dqm_.EndSubRun(static_cast<int>(subrun.run()),
                 static_cast<int>(subrun.subRun()));
}

void CRVStatusDQMAnalyzer::endJob()
{
  dqm_.WriteGraphs();

  if (diagLevel_ > 0) {
    std::cout << "[CRVStatusDQMAnalyzer] Total events: " << dqm_.nEvents()
              << std::endl;
    std::cout << "[CRVStatusDQMAnalyzer] Events with ROC header: "
              << dqm_.nEventsWithRocHeader() << std::endl;
    std::cout << "[CRVStatusDQMAnalyzer] Events with firmware error bit: "
              << dqm_.nEventsWithAnyErrorBit() << std::endl;
    std::cout << "[CRVStatusDQMAnalyzer] Events with DAQ unpack error: "
              << dqm_.nEventsWithDaqErrors() << std::endl;
    std::cout << "[CRVStatusDQMAnalyzer] Distinct ROCs: "
              << dqm_.seenRocs().size() << std::endl;
    if (dqm_.nActiveFEBsMin() != 65535) {
      std::cout << "[CRVStatusDQMAnalyzer] Active FEBs min/mean/max: "
                << dqm_.nActiveFEBsMin() << " / " << dqm_.nActiveFEBsMean()
                << " / " << dqm_.nActiveFEBsMax() << std::endl;
    }
    for (int b = 0; b < CRVStatusDQM::kNErrorBits; ++b) {
      std::cout << "[CRVStatusDQMAnalyzer] " << CRVStatusDQM::errorBitLabel(b)
                << ": " << dqm_.errorBitCount(b) << std::endl;
    }
    for (const auto& roc : dqm_.seenRocs()) {
      std::cout << "[CRVStatusDQMAnalyzer] seen DTC " << static_cast<int>(roc.first)
                << " link " << static_cast<int>(roc.second) << std::endl;
    }
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CRVStatusDQMAnalyzer)
