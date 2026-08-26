#ifndef DQMHelpers_inc_CRVStatusDQM_hh
#define DQMHelpers_inc_CRVStatusDQM_hh
//
// Standalone CRV ROC-status DQM helper. Books and fills firmware-health
// histograms used by both the otsdaq online monitor (CrvStatusMetrics) and
// the offline CRVStatusDQMAnalyzer. Does not fold roc==4 onto roc==2:
// status is per DTC link.
//
// Original Author: R. Mina
//

#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art_root_io/TFileDirectory.h"

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
  };

  static constexpr int kNErrorBits = 8;
  static constexpr int kErrorBitOffset = 24;
  static constexpr int kNPortFlags = 24;
  static constexpr int kNLinksPerDTC = 6;
  static constexpr int kNRocBins = 18; // CRVId::nROC; no GeometryService
  static constexpr int kNDaqErrorCodes = 7;

  static const char* errorBitLabel(int bitIndex);
  // Y-axis of errorBitsVsRoc: dtcId * 6 + linkId. No roc==4 fold.
  static int rocBin(uint8_t dtcId, uint8_t linkId);

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
  void EndSubRun(int run, int subrun);
  void WriteGraphs();

  TH1F* nRocHeaders() const { return h_nRocHeaders_; }
  TH1F* activeFebCount() const { return h_activeFebCount_; }
  TH1F* triggerCount() const { return h_triggerCount_; }
  TH1F* wordCount() const { return h_wordCount_; }
  TH1F* linkLatency() const { return h_linkLatency_; }
  TH1F* errorBits() const { return h_errorBits_; }
  TH2F* errorBitsVsRoc() const { return h_errorBitsVsRoc_; }
  TH1F* portFlags() const { return h_portFlags_; }
  TH1F* rocCensus() const { return h_rocCensus_; }
  TH1F* eventHasError() const { return h_eventHasError_; }
  TH1F* eventHasDaqError() const { return h_eventHasDaqError_; }
  TH1F* daqErrorCode() const { return h_daqErrorCode_; }
  TH1F* ewtMismatch() const { return h_ewtMismatch_; }
  TH1F* errorsPerSubrun() const { return h_errorsPerSubrun_; }
  TH1F* meanLatencyPerSubrun() const { return h_meanLatencyPerSubrun_; }
  TGraph* errorsVsSubrun() const { return g_errorsVsSubrun_; }
  TGraph* meanLatencyVsSubrun() const { return g_meanLatencyVsSubrun_; }

  const std::map<std::pair<uint8_t, uint8_t>, TH1F*>& linkLatencyByRoc() const
  {
    return h_linkLatencyByRoc_;
  }

  const std::vector<RocSnapshot>& lastEventRocs() const { return lastEventRocs_; }

  std::size_t nEvents() const { return nEvents_; }
  std::size_t nEventsWithRocHeader() const { return nEventsWithRocHeader_; }
  std::size_t nEventsWithAnyErrorBit() const { return nEventsWithAnyErrorBit_; }
  std::size_t nEventsWithDaqErrors() const { return nEventsWithDaqErrors_; }
  std::size_t nRocHeadersTotal() const { return nRocHeadersTotal_; }
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
  TH1F* latencyHistFor(uint8_t dtcId, uint8_t linkId);

  Config config_;
  bool booked_{false};
  std::optional<art::TFileDirectory> dir_;

  TH1F* h_nRocHeaders_{nullptr};
  TH1F* h_activeFebCount_{nullptr};
  TH1F* h_triggerCount_{nullptr};
  TH1F* h_wordCount_{nullptr};
  TH1F* h_linkLatency_{nullptr};
  TH1F* h_errorBits_{nullptr};
  TH2F* h_errorBitsVsRoc_{nullptr};
  TH1F* h_portFlags_{nullptr};
  TH1F* h_rocCensus_{nullptr};
  TH1F* h_eventHasError_{nullptr};
  TH1F* h_eventHasDaqError_{nullptr};
  TH1F* h_daqErrorCode_{nullptr};
  TH1F* h_ewtMismatch_{nullptr};
  TH1F* h_errorsPerSubrun_{nullptr};
  TH1F* h_meanLatencyPerSubrun_{nullptr};
  TGraph* g_errorsVsSubrun_{nullptr};
  TGraph* g_meanLatencyVsSubrun_{nullptr};
  std::map<std::pair<uint8_t, uint8_t>, TH1F*> h_linkLatencyByRoc_;

  std::vector<RocSnapshot> lastEventRocs_;

  std::size_t nEvents_{0};
  std::size_t nEventsWithRocHeader_{0};
  std::size_t nEventsWithAnyErrorBit_{0};
  std::size_t nEventsWithDaqErrors_{0};
  std::size_t nRocHeadersTotal_{0};
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

} // namespace mu2e

#endif /* DQMHelpers_inc_CRVStatusDQM_hh */
