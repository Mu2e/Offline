#ifndef DQMHelpers_inc_CRVStatusDQM_hh
#define DQMHelpers_inc_CRVStatusDQM_hh
// Standalone CRV ROC-status DQM helper: firmware-health histograms for the
// otsdaq online monitor and CrvDQMcollector. Status is per DTC link.
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/DQMSegmentation.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace mu2e {

class CRVStatusDQM {
public:
  struct Config {
    int nBinsLatency{1024};
    float maxLinkLatency{4096.f};
    int nBinsTriggerCount{256};
    float maxTriggerCount{65535.f};
    int nBinsWordCount{256};
    float maxWordCount{65535.f};
    int nBinsEwtMismatch{201};
    float maxEwtMismatch{100.f};
    int nBinsErrorsPerSubrun{100};
    float maxErrorsPerSubrun{10000.f};
    // TGraphs vs subrun. Online monitor only — they do not survive hadd.
    bool fillLivePlots{false};
    // Which histograms get per-subrun and last-N-events copies. Default is
    // job-only, which reproduces the pre-segmentation output exactly.
    DQMSegmentation::Config segmentation{};
  };

  static constexpr int kNErrorBits = 8;
  static constexpr int kErrorBitOffset = 24;
  static constexpr int kNPortFlags = static_cast<int>(CRVId::nFEBPerROC);
  static constexpr int kNLinksPerDTC = static_cast<int>(CRVId::nROCPerDTC);
  static constexpr int kNRocBins = static_cast<int>(CRVId::nROC);
  static constexpr int kNDaqErrorCodes = 7;

  static const char* errorBitLabel(int bitIndex);
  // Y-axis of errorBitsVsRoc: dtcId * nROCPerDTC + linkId. No roc==4 fold.
  static int rocBin(uint8_t dtcId, uint8_t linkId);
  static bool rocIndexed(uint8_t dtcId, uint8_t linkId);

  struct RocSnapshot {
    uint8_t dtcId{0};
    uint8_t linkId{0};
    uint64_t ewt{0};
    uint16_t triggerCount{0};
    uint16_t wordCount{0};
    uint16_t activeFebCount{0};
    uint16_t linkLatency{0};
    uint32_t microBunchStatus{0};
  };

  explicit CRVStatusDQM(const Config& config);

  void Book(art::TFileDirectory dir);
  void Fill(const CrvStatusCollection& crvStatus);
  void Fill(const CrvStatusCollection& crvStatus,
            const CrvDAQerrorCollection& crvDaqErrors);
  void BeginSubRun(int run, int subrun);
  void EndSubRun(int run, int subrun);
  void WriteGraphs();

  // Every histogram this helper books, and the segment copies of each.
  DQMSegmentation& segments() { return segments_; }
  const DQMSegmentation& segments() const { return segments_; }

  TH1F* nEventsHist() const { return h_nEvents_; }  //one count per event; hadd-safe
  TH1F* nRocHeaders() const { return h_nRocHeaders_; }  //ROC headers per event
  TH1F* activeFebCount() const { return h_activeFebCount_; }  //active FEBs per ROC header
  TH1F* triggerCount() const { return h_triggerCount_; }  //ROC TriggerCount word
  TH1F* wordCount() const { return h_wordCount_; }  //ROC ControllerEventWordCount word
  TH1F* linkLatency() const { return h_linkLatency_; }  //DTC link latency, all links
  TH1F* errorBits() const { return h_errorBits_; }  //firmware error bits 24-31, labelled
  TH2F* errorBitsVsRoc() const { return h_errorBitsVsRoc_; }  //those bits vs dtcId*6+linkId
  TH1F* portFlags() const { return h_portFlags_; }  //per-port problem flags, bits 0-23
  TH1F* rocCensus() const { return h_rocCensus_; }  //ROC headers seen per dtcId*6+linkId
  TH1F* eventHasError() const { return h_eventHasError_; }  //0/1 per event: any firmware bit
  TH1F* eventHasDaqError() const { return h_eventHasDaqError_; }  //0/1 per event: unpack error
  TH1F* daqErrorCode() const { return h_daqErrorCode_; }  //CrvDAQerror code, labelled
  TH1F* ewtMismatch() const { return h_ewtMismatch_; }  //ROC EWT minus DTC EWT
  TH1F* errorsPerSubrun() const { return h_errorsPerSubrun_; }  //error events in a subrun
  TH1F* meanLatencyPerSubrun() const { return h_meanLatencyPerSubrun_; }  //mean latency in a subrun
  TGraph* errorsVsSubrun() const { return g_errorsVsSubrun_; }  //error events vs subrun
  TGraph* meanLatencyVsSubrun() const { return g_meanLatencyVsSubrun_; }  //mean latency vs subrun

  //link latency, one hist per (dtcId, linkId)
  const std::map<std::pair<uint8_t, uint8_t>, DQMHist1<TH1F>>& linkLatencyByRoc() const
  {
    return h_linkLatencyByRoc_;
  }

  const std::vector<RocSnapshot>& lastEventRocs() const { return lastEventRocs_; }

  std::size_t nEvents() const { return nEvents_; }
  std::size_t nEventsWithRocHeader() const { return nEventsWithRocHeader_; }
  std::size_t nEventsWithAnyErrorBit() const { return nEventsWithAnyErrorBit_; }
  std::size_t nEventsWithDaqErrors() const { return nEventsWithDaqErrors_; }
  std::size_t nRocHeadersTotal() const { return nRocHeadersTotal_; }
  long long nUnindexedRocs() const { return nUnindexedRocs_; }
  uint16_t nActiveFEBsMin() const { return nActiveFEBsMin_; }
  uint16_t nActiveFEBsMax() const { return nActiveFEBsMax_; }
  double nActiveFEBsMean() const;
  const std::set<std::pair<uint8_t, uint8_t>>& seenRocs() const
  {
    return seenRocs_;
  }
  std::size_t errorBitCount(int bitIndex) const;
  bool booked() const { return booked_; }

private:
  void fillDaqErrors(const CrvDAQerrorCollection& crvDaqErrors);
  void persistGraph(TGraph* g);
  DQMHist1<TH1F> latencyHistFor(uint8_t dtcId, uint8_t linkId);
  void noteUnindexedRoc(uint8_t dtcId, uint8_t linkId);

  Config config_;
  bool booked_{false};
  std::optional<art::TFileDirectory> dir_;
  DQMSegmentation segments_;

  DQMHist1<TH1F> h_nEvents_;
  DQMHist1<TH1F> h_nRocHeaders_;
  DQMHist1<TH1F> h_activeFebCount_;
  DQMHist1<TH1F> h_triggerCount_;
  DQMHist1<TH1F> h_wordCount_;
  DQMHist1<TH1F> h_linkLatency_;
  DQMHist1<TH1F> h_errorBits_;
  DQMHist2<TH2F> h_errorBitsVsRoc_;
  DQMHist1<TH1F> h_portFlags_;
  DQMHist1<TH1F> h_rocCensus_;
  DQMHist1<TH1F> h_eventHasError_;
  DQMHist1<TH1F> h_eventHasDaqError_;
  DQMHist1<TH1F> h_daqErrorCode_;
  DQMHist1<TH1F> h_ewtMismatch_;
  DQMHist1<TH1F> h_errorsPerSubrun_;
  DQMHist1<TH1F> h_meanLatencyPerSubrun_;
  TGraph* g_errorsVsSubrun_{nullptr};
  TGraph* g_meanLatencyVsSubrun_{nullptr};
  std::map<std::pair<uint8_t, uint8_t>, DQMHist1<TH1F>> h_linkLatencyByRoc_;

  std::vector<RocSnapshot> lastEventRocs_;

  std::size_t nEvents_{0};
  std::size_t nEventsWithRocHeader_{0};
  std::size_t nEventsWithAnyErrorBit_{0};
  std::size_t nEventsWithDaqErrors_{0};
  std::size_t nRocHeadersTotal_{0};
  long long nUnindexedRocs_{0};
  bool warnedUnindexedRoc_{false};
  std::size_t errorBitCounts_[kNErrorBits]{};
  uint16_t nActiveFEBsMin_{std::numeric_limits<uint16_t>::max()};
  uint16_t nActiveFEBsMax_{0};
  std::uint64_t nActiveFEBsSum_{0};
  std::size_t nActiveFEBsSamples_{0};
  std::set<std::pair<uint8_t, uint8_t>> seenRocs_;

  std::size_t nEventsThisSubrun_{0};
  std::size_t nEventsWithAnyErrorBitThisSubrun_{0};
  std::uint64_t latencySumThisSubrun_{0};
  std::size_t latencyNThisSubrun_{0};
};

// At namespace scope rather than nested in the class: a nested sibling
// cannot default-construct Config, whose member initializers are only
// complete at the end of the enclosing class.
// Validated FHiCL for the Config above; see the note on CRVDigiDQMFhicl.
struct CRVStatusDQMFhicl {
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  // Not constexpr: Config carries the segmentation block, which holds
  // std::string and std::vector, so it is not a literal type.
  static inline const CRVStatusDQM::Config d{};

  fhicl::Atom<int> nBinsLatency{Name("nBinsLatency"), Comment("Bins for linkLatency"), d.nBinsLatency};
  fhicl::Atom<float> maxLinkLatency{Name("maxLinkLatency"), Comment("Upper edge for linkLatency"), d.maxLinkLatency};
  fhicl::Atom<int> nBinsTriggerCount{Name("nBinsTriggerCount"), Comment("Bins for triggerCount"), d.nBinsTriggerCount};
  fhicl::Atom<float> maxTriggerCount{Name("maxTriggerCount"), Comment("Upper edge for triggerCount"), d.maxTriggerCount};
  fhicl::Atom<int> nBinsWordCount{Name("nBinsWordCount"), Comment("Bins for wordCount"), d.nBinsWordCount};
  fhicl::Atom<float> maxWordCount{Name("maxWordCount"), Comment("Upper edge for wordCount"), d.maxWordCount};
  fhicl::Atom<int> nBinsEwtMismatch{Name("nBinsEwtMismatch"), Comment("Bins for ewtMismatch"), d.nBinsEwtMismatch};
  fhicl::Atom<float> maxEwtMismatch{Name("maxEwtMismatch"), Comment("+/- range for ewtMismatch"), d.maxEwtMismatch};
  fhicl::Atom<int> nBinsErrorsPerSubrun{Name("nBinsErrorsPerSubrun"), Comment("Bins for errorsPerSubrun"), d.nBinsErrorsPerSubrun};
  fhicl::Atom<float> maxErrorsPerSubrun{Name("maxErrorsPerSubrun"), Comment("Upper edge for errorsPerSubrun"), d.maxErrorsPerSubrun};
  fhicl::Atom<bool> fillLivePlots{
      Name("fillLivePlots"), Comment("Book TGraphs vs subrun (online only; not hadd-safe)"),
      d.fillLivePlots};
  fhicl::OptionalDelegatedParameter segmentation{
      Name("segmentation"),
      Comment("Per-subrun / last-N-events copies and publishing for this helper; "
              "see Offline/DQMHelpers/README.md")};
};

CRVStatusDQM::Config toConfig(const CRVStatusDQMFhicl& c);

} // namespace mu2e

#endif /* DQMHelpers_inc_CRVStatusDQM_hh */
