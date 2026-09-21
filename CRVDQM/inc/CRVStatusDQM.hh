#ifndef CRVDQM_inc_CRVStatusDQM_hh
#define CRVDQM_inc_CRVStatusDQM_hh
// CRV ROC-status DQM client: firmware-health histograms for the otsdaq online
// monitor and the offline DQM modules. Status is per DTC link; every per-link
// object is booked up front for all CRVDQMRun1::kNLinks links, indexed
// dtcId * nROCPerDTC + linkId.
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVDQMRun1.hh"
#include "Offline/DQMHelpers/inc/DQMClient.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "TH1F.h"
#include "TH2F.h"

#include <cstdint>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace mu2e {

class CRVStatusDQM : public DQMClient {
public:
  static constexpr int kBinningVersion = 1;

  static constexpr int kNErrorBits = 8;       //MicroBunchStatus bits 24-31
  static constexpr int kErrorBitOffset = 24;
  static constexpr int kNPortFlags = 24;      //MicroBunchStatus bits 0-23
  static constexpr int kNDaqErrorCodes = 32;  //room for codes added upstream
  static constexpr uint16_t kLatencySentinel = 0xFFFF;
  static constexpr double kLatencyTickToUs = 0.064;

  static constexpr DQMAxis kRocHeaders = DQMAxis::Counts(0, 19);
  static constexpr DQMAxis kActiveFebs = DQMAxis::Counts(0, 24);
  static constexpr DQMAxis kTriggerCount{256, 0., 65535.};
  static constexpr DQMAxis kWordCount{512, 0., 32768.};
  static constexpr DQMAxis kLatency{1024, 0., 4096.};
  static constexpr DQMAxis kErrorBit = DQMAxis::Counts(0, kNErrorBits - 1);
  static constexpr DQMAxis kPortFlag = DQMAxis::Counts(0, kNPortFlags - 1);
  static constexpr DQMAxis kLink = DQMAxis::Counts(0, CRVDQMRun1::kNLinks - 1);
  static constexpr DQMAxis kBool = DQMAxis::Counts(0, 1);
  static constexpr DQMAxis kDaqErrorCode = DQMAxis::Counts(0, kNDaqErrorCodes - 1);
  static constexpr DQMAxis kEwtMismatch = DQMAxis::Counts(-100, 100);
  // Fraction, not count: a count's range depends on the subrun length.
  static constexpr DQMAxis kErrorFraction{101, -0.005, 1.005};
  // Per-link online views.
  static constexpr DQMAxis kLatencyUs{300, 0., 150.};
  static constexpr DQMAxis kLatencyUs2{2000, 0., 10000.};

  static const char* errorBitLabel(int bitIndex);

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

  explicit CRVStatusDQM(const DQMHistSet::Config& hists = {});

  void Fill(const CrvStatusCollection& crvStatus);
  void Fill(const CrvStatusCollection& crvStatus, const CrvDAQerrorCollection& crvDaqErrors);

  TH1F* nRocHeaders() const { return h_nRocHeaders_; }  //ROC headers per event
  TH1F* activeFebCount() const { return h_activeFebCount_; }  //active FEBs per ROC header
  TH1F* triggerCount() const { return h_triggerCount_; }  //ROC TriggerCount word
  TH1F* wordCount() const { return h_wordCount_; }  //ROC ControllerEventWordCount word
  TH1F* linkLatency() const { return h_linkLatency_; }  //DTC link latency, all links
  TH1F* linkLatencyInvalid() const { return h_linkLatencyInvalid_; }  //0xFFFF latency words per link
  TH1F* statusBlocksByLink() const { return h_statusBlocksByLink_; }  //its denominator
  TH1F* errorBits() const { return h_errorBits_; }  //firmware error bits 24-31, labelled
  TH2F* errorBitsVsRoc() const { return h_errorBitsVsRoc_; }  //those bits vs link
  TH1F* portFlags() const { return h_portFlags_; }  //per-port problem flags, bits 0-23
  TH1F* rocCensus() const { return h_rocCensus_; }  //ROC headers seen per link
  TH1F* eventHasError() const { return h_eventHasError_; }  //0/1 per event: any firmware bit
  TH1F* eventHasDaqError() const { return h_eventHasDaqError_; }  //0/1 per event: unpack error
  TH1F* daqErrorCode() const { return h_daqErrorCode_; }  //CrvDAQerror code, labelled
  TH1F* ewtMismatch() const { return h_ewtMismatch_; }  //ROC EWT minus DTC EWT
  TH1F* errorFractionPerSubrun() const { return h_errorFractionPerSubrun_; }
  TH1F* meanLatencyPerSubrun() const { return h_meanLatencyPerSubrun_; }
  //link latency, one hist per link
  TH1F* linkLatencyByLink(int link) const;

  const std::vector<RocSnapshot>& lastEventRocs() const { return lastEventRocs_; }

  std::size_t nEventsWithRocHeader() const { return nEventsWithRocHeader_; }
  std::size_t nEventsWithAnyErrorBit() const { return nEventsWithAnyErrorBit_; }
  std::size_t nEventsWithDaqErrors() const { return nEventsWithDaqErrors_; }
  std::size_t nRocHeadersTotal() const { return nRocHeadersTotal_; }
  uint16_t nActiveFEBsMin() const { return nActiveFEBsMin_; }
  uint16_t nActiveFEBsMax() const { return nActiveFEBsMax_; }
  double nActiveFEBsMean() const;
  const std::set<std::pair<uint8_t, uint8_t>>& seenRocs() const { return seenRocs_; }
  std::size_t errorBitCount(int bitIndex) const;

private:
  // Online per-link views, all booked in book().
  struct LinkHists {
    DQMH1<TH1F> latencyUs;
    DQMH1<TH1F> latencyUs2;
    DQMH1<TH1F> linkStatus;
    DQMH1<TH1F> rocStatus;
    DQMH1<TH1F> rocGroup;
    DQMH1<TH1F> portFlagBits;
  };

  void book() override;
  void endSubRun() override;
  void resetForNewRun() override;

  void fillDaqErrors(const CrvDAQerrorCollection& crvDaqErrors);
  void fillLink(const CrvStatus& status, int link);
  DQMSeries& linkSeries(const char* nameFmt, const char* titleFmt, int link);
  DQMSeries& bitSeries(const char* prefix, const char* bitName, int bit, int link);

  DQMH1<TH1F> h_nRocHeaders_;
  DQMH1<TH1F> h_activeFebCount_;
  DQMH1<TH1F> h_triggerCount_;
  DQMH1<TH1F> h_wordCount_;
  DQMH1<TH1F> h_linkLatency_;
  DQMH1<TH1F> h_linkLatencyInvalid_;
  DQMH1<TH1F> h_statusBlocksByLink_;
  DQMH1<TH1F> h_errorBits_;
  DQMH2<TH2F> h_errorBitsVsRoc_;
  DQMH1<TH1F> h_portFlags_;
  DQMH1<TH1F> h_rocCensus_;
  DQMH1<TH1F> h_eventHasError_;
  DQMH1<TH1F> h_eventHasDaqError_;
  DQMH1<TH1F> h_daqErrorCode_;
  DQMH1<TH1F> h_ewtMismatch_;
  DQMH1<TH1F> h_errorFractionPerSubrun_;
  DQMH1<TH1F> h_meanLatencyPerSubrun_;
  std::vector<DQMH1<TH1F>> h_linkLatencyByLink_;
  std::vector<LinkHists> linkHists_;
  DQMH2<TH2F> h2_rocStatusSummary_;
  DQMH2<TH2F> h2_rocGroupSummary_;
  DQMH2<TH2F> h2_portFlagBits_;
  DQMSeries* g_errorFractionVsSubrun_{nullptr};
  DQMSeries* g_meanLatencyVsSubrun_{nullptr};

  std::vector<RocSnapshot> lastEventRocs_;

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

#endif /* CRVDQM_inc_CRVStatusDQM_hh */
