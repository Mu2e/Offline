// Standalone CRV ROC-status DQM helper.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/CRVStatusDQM.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TString.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <optional>
#include <string>

namespace mu2e {

namespace {

const char* kErrorBitLabels[CRVStatusDQM::kNErrorBits] = {
    "FEBuBMismatch",
    "FEBBufferIssue",
    "FEBOverflow",
    "Group1Issue",
    "Group2Issue",
    "Group3Issue",
    "uBMatchError",
    "Truncation"};

} // namespace

const char* CRVStatusDQM::errorBitLabel(int bitIndex)
{
  if (bitIndex < 0 || bitIndex >= kNErrorBits) {
    return "";
  }
  return kErrorBitLabels[bitIndex];
}

int CRVStatusDQM::rocBin(uint8_t dtcId, uint8_t linkId)
{
  return static_cast<int>(dtcId) * kNLinksPerDTC + static_cast<int>(linkId);
}

bool CRVStatusDQM::rocIndexed(uint8_t dtcId, uint8_t linkId)
{
  if (static_cast<int>(linkId) >= kNLinksPerDTC) {
    return false;
  }
  const int y = rocBin(dtcId, linkId);
  return y >= 0 && y < kNRocBins;
}

void CRVStatusDQM::noteUnindexedRoc(uint8_t dtcId, uint8_t linkId)
{
  ++nUnindexedRocs_;
  if (warnedUnindexedRoc_) {
    return;
  }
  warnedUnindexedRoc_ = true;
  mf::LogWarning("CRVStatusDQM")
      << "status from DTC " << static_cast<int>(dtcId) << " link "
      << static_cast<int>(linkId) << " is outside the " << kNRocBins
      << "-bin ROC axis (link 0-" << (kNLinksPerDTC - 1) << ", dtcId*"
      << kNLinksPerDTC << "+linkId < " << kNRocBins
      << "). Per-link latency / rocCensus / errorBitsVsRoc skip it. "
      << "Reported once per job.";
}

CRVStatusDQM::CRVStatusDQM(const Config& config) : config_(config)
{
  segments_.SetConfig(config_.segmentation);
}

void CRVStatusDQM::Book(art::TFileDirectory dir)
{
  dir_ = dir;
  segments_.Book(dir);

  h_nEvents_ =
      segments_.book1<TH1F>("nEvents", "Events processed;;Events", 1, 0.5, 1.5);

  h_nRocHeaders_ = segments_.book1<TH1F>(
      "nRocHeaders", "ROC headers per event;N(ROC headers);Events", 20, -0.5, 19.5);

  h_activeFebCount_ = segments_.book1<TH1F>(
      "activeFebCount",
      "Active FEBs per ROC header;N(active FEBs);ROC headers",
      25,
      -0.5,
      24.5);

  h_triggerCount_ = segments_.book1<TH1F>(
      "triggerCount", "ROC TriggerCount;TriggerCount;ROC headers",
      config_.nBinsTriggerCount, 0, config_.maxTriggerCount);

  h_wordCount_ = segments_.book1<TH1F>(
      "wordCount", "ROC ControllerEventWordCount;WordCount;ROC headers",
      config_.nBinsWordCount, 0, config_.maxWordCount);

  h_linkLatency_ = segments_.book1<TH1F>(
      "linkLatency", "DTC link latency;Latency;Status blocks",
      config_.nBinsLatency, 0, config_.maxLinkLatency);

  h_errorBits_ = segments_.book1<TH1F>(
      "errorBits",
      "Firmware error-bit occupancy;Error bit;Counts",
      kNErrorBits,
      -0.5,
      kNErrorBits - 0.5);

  h_errorBitsVsRoc_ = segments_.book2<TH2F>(
      "errorBitsVsRoc",
      "Firmware error bits vs ROC;Error bit;DTC#times6 + link ID",
      kNErrorBits,
      -0.5,
      kNErrorBits - 0.5,
      kNRocBins,
      -0.5,
      kNRocBins - 0.5);

  h_errorBits_.ForEach([](TH1F* h) {
    for (int i = 0; i < kNErrorBits; ++i) {
      h->GetXaxis()->SetBinLabel(i + 1, kErrorBitLabels[i]);
    }
  });
  h_errorBitsVsRoc_.ForEach([](TH2F* h) {
    for (int i = 0; i < kNErrorBits; ++i) {
      h->GetXaxis()->SetBinLabel(i + 1, kErrorBitLabels[i]);
    }
  });

  h_portFlags_ = segments_.book1<TH1F>(
      "portFlags", "MicroBunch port flags (bits 0-23);Port;Counts",
      kNPortFlags, -0.5, kNPortFlags - 0.5);

  h_rocCensus_ = segments_.book1<TH1F>(
      "rocCensus",
      "ROC headers by DTC#times6 + link;DTC#times6 + link ID;ROC headers",
      kNRocBins,
      -0.5,
      kNRocBins - 0.5);

  h_eventHasError_ = segments_.book1<TH1F>(
      "eventHasError",
      "Event has any firmware error bit;0=ok 1=error;Events",
      2,
      -0.5,
      1.5);

  h_eventHasDaqError_ = segments_.book1<TH1F>(
      "eventHasDaqError",
      "Event has unpack DAQ error (excl. wrongSubsystemID);0=ok 1=error;Events",
      2,
      -0.5,
      1.5);

  h_daqErrorCode_ = segments_.book1<TH1F>(
      "daqErrorCode",
      "CrvDAQerror code;Error code;Counts",
      kNDaqErrorCodes,
      -0.5,
      kNDaqErrorCodes - 0.5);
  h_daqErrorCode_.ForEach([](TH1F* h) {
    h->GetXaxis()->SetBinLabel(1, "unknown");
    h->GetXaxis()->SetBinLabel(2, "unableToGetDataBlock");
    h->GetXaxis()->SetBinLabel(3, "invalidPacket");
    h->GetXaxis()->SetBinLabel(4, "wrongSubsystemID");
    h->GetXaxis()->SetBinLabel(5, "errorUnpackingStatusPacket");
    h->GetXaxis()->SetBinLabel(6, "errorUnpackingCrvHits");
    h->GetXaxis()->SetBinLabel(7, "byteCountMismatch");
  });

  h_ewtMismatch_ = segments_.book1<TH1F>(
      "ewtMismatch", "ROC EWT - DTC EWT;#Delta EWT;ROC headers",
      config_.nBinsEwtMismatch, -config_.maxEwtMismatch - 0.5f,
      config_.maxEwtMismatch + 0.5f);

  h_errorsPerSubrun_ = segments_.book1<TH1F>(
      "errorsPerSubrun",
      "Events with any firmware error bit, per subrun;Events with error;Subruns",
      config_.nBinsErrorsPerSubrun,
      -0.5f,
      config_.maxErrorsPerSubrun + 0.5f);

  h_meanLatencyPerSubrun_ = segments_.book1<TH1F>(
      "meanLatencyPerSubrun",
      "Mean link latency per subrun;Mean latency;Subruns",
      config_.nBinsLatency,
      0,
      config_.maxLinkLatency);

  if (config_.fillLivePlots) {
    g_errorsVsSubrun_ = dir.make<TGraph>();
    g_errorsVsSubrun_->SetName("g_errorsVsSubrun");
    g_errorsVsSubrun_->SetTitle(
        "Events with firmware error vs subrun;Subrun;Events with error");

    g_meanLatencyVsSubrun_ = dir.make<TGraph>();
    g_meanLatencyVsSubrun_->SetName("g_meanLatencyVsSubrun");
    g_meanLatencyVsSubrun_->SetTitle(
        "Mean link latency vs subrun;Subrun;Mean latency");
  }

  booked_ = true;
}

DQMHist1<TH1F> CRVStatusDQM::latencyHistFor(uint8_t dtcId, uint8_t linkId)
{
  if (!booked_ || !rocIndexed(dtcId, linkId)) {
    return DQMHist1<TH1F>();
  }
  const auto key = std::make_pair(dtcId, linkId);
  auto it = h_linkLatencyByRoc_.find(key);
  if (it != h_linkLatencyByRoc_.end()) {
    return it->second;
  }
  const std::string name =
      Form("linkLatency_dtc%u_roc%u", static_cast<unsigned>(dtcId),
           static_cast<unsigned>(linkId));
  const std::string title =
      Form("DTC %u ROC %u link latency;Latency;Status blocks",
           static_cast<unsigned>(dtcId),
           static_cast<unsigned>(linkId));
  DQMHist1<TH1F> h = segments_.book1<TH1F>(
      name, title, config_.nBinsLatency, 0, config_.maxLinkLatency);
  h_linkLatencyByRoc_[key] = h;
  return h;
}

void CRVStatusDQM::Fill(const CrvStatusCollection& crvStatus)
{
  ++nEvents_;
  ++nEventsThisSubrun_;

  // Rotate the segmentation ring between events, before anything is filled.
  const bool haveEwt = !crvStatus.empty();
  segments_.Advance(nEvents_,
                    haveEwt ? std::optional<uint64_t>(
                                  crvStatus.front().GetEventWindowTag())
                            : std::nullopt);

  h_nEvents_.Fill(1.f);
  lastEventRocs_.clear();

  int nHeadersThisEvent = 0;
  bool anyErrorThisEvent = false;

  for (const auto& status : crvStatus) {
    const uint8_t dtcId = status.GetDTCID();
    const uint8_t linkId = status.GetLinkID();
    const uint16_t latency = status.GetLinkLatency();
    const bool indexed = rocIndexed(dtcId, linkId);
    if (!indexed) {
      noteUnindexedRoc(dtcId, linkId);
    }

    latencySumThisSubrun_ += latency;
    ++latencyNThisSubrun_;

    if (booked_) {
      h_linkLatency_.Fill(latency);
      if (indexed) {
        latencyHistFor(dtcId, linkId).Fill(latency);
      }
    }

    if (!status.HasROCHeader()) {
      continue;
    }

    const auto& headers = status.GetROCHeader();
    const auto& roc = headers.front();

    ++nHeadersThisEvent;
    ++nRocHeadersTotal_;
    seenRocs_.insert({dtcId, linkId});
    const int ybin = rocBin(dtcId, linkId);
    if (indexed) {
      h_rocCensus_.Fill(ybin);
    }

    const std::bitset<kNPortFlags> activeFEBs = roc.GetActiveFEBFlags();
    const uint16_t nActive = static_cast<uint16_t>(activeFEBs.count());
    const uint16_t trigCount = roc.TriggerCount;
    const uint16_t wordCount = roc.ControllerEventWordCount;
    const uint32_t ubStatus = roc.GetMicroBunchStatus();
    const uint64_t rocEwt = roc.GetEventWindowTag();
    const uint64_t dtcEwt = status.GetEventWindowTag();

    RocSnapshot snap;
    snap.dtcId = dtcId;
    snap.linkId = linkId;
    snap.ewt = rocEwt;
    snap.triggerCount = trigCount;
    snap.wordCount = wordCount;
    snap.activeFebCount = nActive;
    snap.linkLatency = latency;
    snap.microBunchStatus = ubStatus;
    lastEventRocs_.push_back(snap);

    nActiveFEBsMin_ = std::min(nActiveFEBsMin_, nActive);
    nActiveFEBsMax_ = std::max(nActiveFEBsMax_, nActive);
    nActiveFEBsSum_ += nActive;
    ++nActiveFEBsSamples_;

    h_activeFebCount_.Fill(nActive);
    h_triggerCount_.Fill(trigCount);
    h_wordCount_.Fill(wordCount);
    h_ewtMismatch_.Fill(static_cast<float>(static_cast<int64_t>(rocEwt) -
                                           static_cast<int64_t>(dtcEwt)));

    const uint32_t portFlags = ubStatus & 0x00FFFFFFu;
    for (int p = 0; p < kNPortFlags; ++p) {
      if ((portFlags >> p) & 1u) {
        h_portFlags_.Fill(p);
      }
    }

    for (int b = 0; b < kNErrorBits; ++b) {
      if ((ubStatus >> (kErrorBitOffset + b)) & 1u) {
        ++errorBitCounts_[b];
        anyErrorThisEvent = true;
        h_errorBits_.Fill(b);
        if (indexed) {
          h_errorBitsVsRoc_.Fill(b, ybin);
        }
      }
    }
  }

  h_nRocHeaders_.Fill(nHeadersThisEvent);
  if (nHeadersThisEvent > 0) {
    ++nEventsWithRocHeader_;
  }
  if (anyErrorThisEvent) {
    ++nEventsWithAnyErrorBit_;
    ++nEventsWithAnyErrorBitThisSubrun_;
  }
  h_eventHasError_.Fill(anyErrorThisEvent ? 1.f : 0.f);
}

void CRVStatusDQM::Fill(const CrvStatusCollection& crvStatus,
                        const CrvDAQerrorCollection& crvDaqErrors)
{
  Fill(crvStatus);
  fillDaqErrors(crvDaqErrors);
}

void CRVStatusDQM::fillDaqErrors(const CrvDAQerrorCollection& crvDaqErrors)
{
  bool countedEvent = false;
  for (const auto& err : crvDaqErrors) {
    const int code = static_cast<int>(err.GetErrorCode());
    if (code >= 0 && code < kNDaqErrorCodes) {
      h_daqErrorCode_.Fill(code);
    }
    if (err.GetErrorCode() == CrvDAQerrorCode::wrongSubsystemID) {
      continue;
    }
    if (!countedEvent) {
      ++nEventsWithDaqErrors_;
      countedEvent = true;
    }
  }
  h_eventHasDaqError_.Fill(countedEvent ? 1.f : 0.f);
}

void CRVStatusDQM::BeginSubRun(int run, int subrun)
{
  segments_.BeginSubRun(run, subrun);
}

void CRVStatusDQM::EndSubRun(int run, int subrun)
{
  // Set the identity first: a caller that only calls EndSubRun still gets
  // correctly named per-subrun archive copies.
  segments_.BeginSubRun(run, subrun);
  segments_.EndSubRun();

  if (nEventsThisSubrun_ == 0) {
    return;
  }

  const double meanLat = (latencyNThisSubrun_ > 0) ?
                             static_cast<double>(latencySumThisSubrun_) /
                                 static_cast<double>(latencyNThisSubrun_) :
                             0.0;

  h_errorsPerSubrun_.Fill(static_cast<float>(nEventsWithAnyErrorBitThisSubrun_));
  h_meanLatencyPerSubrun_.Fill(static_cast<float>(meanLat));
  if (g_errorsVsSubrun_) {
    g_errorsVsSubrun_->SetPoint(g_errorsVsSubrun_->GetN(),
                                subrun,
                                static_cast<double>(nEventsWithAnyErrorBitThisSubrun_));
  }
  if (g_meanLatencyVsSubrun_) {
    g_meanLatencyVsSubrun_->SetPoint(
        g_meanLatencyVsSubrun_->GetN(), subrun, meanLat);
  }

  nEventsThisSubrun_ = 0;
  nEventsWithAnyErrorBitThisSubrun_ = 0;
  latencySumThisSubrun_ = 0;
  latencyNThisSubrun_ = 0;
}

void CRVStatusDQM::persistGraph(TGraph* g)
{
  if (!dir_ || g == nullptr) {
    return;
  }
  if (g->GetN() <= 0) {
    dir_->makeAndRegister<TGraph>(g->GetName(), g->GetTitle());
    return;
  }
  dir_->makeAndRegister<TGraph>(
      g->GetName(), g->GetTitle(), g->GetN(), g->GetX(), g->GetY());
}

void CRVStatusDQM::WriteGraphs()
{
  if (!booked_ || !dir_) {
    return;
  }
  segments_.Finalize();
  persistGraph(g_errorsVsSubrun_);
  persistGraph(g_meanLatencyVsSubrun_);
}

double CRVStatusDQM::nActiveFEBsMean() const
{
  if (nActiveFEBsSamples_ == 0) {
    return 0.0;
  }
  return static_cast<double>(nActiveFEBsSum_) /
         static_cast<double>(nActiveFEBsSamples_);
}

std::size_t CRVStatusDQM::errorBitCount(int bitIndex) const
{
  if (bitIndex < 0 || bitIndex >= kNErrorBits) {
    return 0;
  }
  return errorBitCounts_[bitIndex];
}

} // namespace mu2e
