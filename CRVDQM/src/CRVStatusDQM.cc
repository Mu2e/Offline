// CRV ROC-status DQM client.
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVStatusDQM.hh"

#include "TString.h"

#include <algorithm>
#include <bitset>
#include <optional>
#include <string>

namespace mu2e {

using namespace CRVDQMRun1;

namespace {

const char* kErrorBitLabels[CRVStatusDQM::kNErrorBits] = {
    "FEBuBMismatch", "FEBBufferIssue", "FEBOverflow", "Group1Issue",
    "Group2Issue",   "Group3Issue",    "uBMatchError", "Truncation"};

// Online link-status layout (otsdaq-mu2e-crv CrvDQM, mu2e/ots_ops).
constexpr int kNLinkStatusBits = 5;
constexpr uint8_t kLinkStatusBits[kNLinkStatusBits] = {0, 2, 3, 6, 7};
constexpr const char* kLinkStatusBitNames[kNLinkStatusBits] = {
    "ROCTimeout", "SeqNumErr", "CRCErr", "FatalErr", "Error"};

// MicroBunchStatus summary bits [24:31], split into error and group flags.
constexpr int kNRocStatusBits = 5;
constexpr uint8_t kRocStatusBits[kNRocStatusBits] = {24, 25, 26, 30, 31};
constexpr const char* kRocStatusBitNames[kNRocStatusBits] = {
    "FEBuBMismatch", "FEBBufferIssue", "FEBOverflow", "uBMatchErr", "Truncation"};

constexpr int kNRocGroupBits = 3;
constexpr uint8_t kRocGroupBits[kNRocGroupBits] = {27, 28, 29};
constexpr const char* kRocGroupBitNames[kNRocGroupBits] = {
    "Group1Issue", "Group2Issue", "Group3Issue"};

// Flag views count from 1 so the labels sit on bins 1..n.
constexpr DQMAxis flagAxis(int n) { return DQMAxis(n, 0.5, n + 0.5); }

template <class H>
void labelX(const H& h, int n, const char* const* labels)
{
  h.ForEach([n, labels](auto* p) {
    for (int i = 0; i < n; ++i) {
      p->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    }
  });
}

} // namespace

const char* CRVStatusDQM::errorBitLabel(int bitIndex)
{
  if (bitIndex < 0 || bitIndex >= kNErrorBits) {
    return "";
  }
  return kErrorBitLabels[bitIndex];
}

CRVStatusDQM::CRVStatusDQM(const DQMHistSet::Config& hists) :
    DQMClient("CRVStatusDQM", kBinningVersion, hists)
{}

TH1F* CRVStatusDQM::linkLatencyByLink(int link) const
{
  if (link < 0 || static_cast<std::size_t>(link) >= h_linkLatencyByLink_.size()) {
    return nullptr;
  }
  return h_linkLatencyByLink_[link];
}

void CRVStatusDQM::book()
{
  DQMHistSet& h = hists();

  h_nRocHeaders_ = h.book1<TH1F>("nRocHeaders", "ROC headers per event;N(ROC headers);Events",
                                 kRocHeaders);
  h_activeFebCount_ = h.book1<TH1F>(
      "activeFebCount", "Active FEBs per ROC header;N(active FEBs);ROC headers", kActiveFebs);
  h_triggerCount_ = h.book1<TH1F>("triggerCount", "ROC TriggerCount;TriggerCount;ROC headers",
                                  kTriggerCount);
  h_wordCount_ = h.book1<TH1F>(
      "wordCount", "ROC ControllerEventWordCount;WordCount;ROC headers", kWordCount);
  h_linkLatency_ = h.book1<TH1F>("linkLatency", "DTC link latency;Latency;Status blocks",
                                 kLatency);
  // Live check for the DTC's "no measurement" latency word, per link. With
  // statusBlocksByLink as denominator it is a rate that merges across files.
  h_linkLatencyInvalid_ = h.book1<TH1F>(
      "linkLatencyInvalid",
      "Link latency 0xFFFF (no measurement);DTC#times6 + link ID;Status blocks", kLink);
  h_statusBlocksByLink_ = h.book1<TH1F>(
      "statusBlocksByLink", "Status blocks by link;DTC#times6 + link ID;Status blocks", kLink);
  h_errorBits_ = h.book1<TH1F>("errorBits", "Firmware error-bit occupancy;Error bit;Counts",
                               kErrorBit);
  h_errorBitsVsRoc_ = h.book2<TH2F>(
      "errorBitsVsRoc", "Firmware error bits vs link;Error bit;DTC#times6 + link ID",
      kErrorBit, kLink);
  labelX(h_errorBits_, kNErrorBits, kErrorBitLabels);
  labelX(h_errorBitsVsRoc_, kNErrorBits, kErrorBitLabels);

  h_portFlags_ = h.book1<TH1F>("portFlags", "MicroBunch port flags (bits 0-23);Port;Counts",
                               kPortFlag);
  h_rocCensus_ = h.book1<TH1F>(
      "rocCensus", "ROC headers by DTC#times6 + link;DTC#times6 + link ID;ROC headers", kLink);
  h_eventHasError_ = h.book1<TH1F>(
      "eventHasError", "Event has any firmware error bit;0=ok 1=error;Events", kBool);
  h_eventHasDaqError_ = h.book1<TH1F>(
      "eventHasDaqError",
      "Event has unpack DAQ error;0=ok 1=error;Events", kBool);
  h_daqErrorCode_ = h.book1<TH1F>("daqErrorCode", "CrvDAQerror code;Error code;Counts",
                                  kDaqErrorCode);
  h_daqErrorCode_.ForEach([](TH1F* p) {
    for (const auto& [code, name] : CrvDAQerrorCodeDetail::names()) {
      const int c = static_cast<int>(code);
      if (c >= 0 && c < kNDaqErrorCodes) {
        p->GetXaxis()->SetBinLabel(c + 1, name.c_str());
      }
    }
  });
  h_ewtMismatch_ = h.book1<TH1F>("ewtMismatch", "ROC EWT - DTC EWT;#Delta EWT;ROC headers",
                                 kEwtMismatch);
  h_errorFractionPerSubrun_ = h.book1<TH1F>(
      "errorFractionPerSubrun",
      "Fraction of events with any firmware error bit, per subrun;Fraction;Subruns",
      kErrorFraction);
  h_meanLatencyPerSubrun_ = h.book1<TH1F>(
      "meanLatencyPerSubrun", "Mean link latency per subrun;Mean latency;Subruns", kLatency);

  h_linkLatencyByLink_.resize(kNLinks);
  linkHists_.resize(kNLinks);
  for (int link = 0; link < kNLinks; ++link) {
    const int dtc = link / kNLinksPerDTC;
    const int roc = link % kNLinksPerDTC;
    h_linkLatencyByLink_[link] = h.book1<TH1F>(
        Form("linkLatency_dtc%d_roc%d", dtc, roc),
        Form("DTC %d ROC %d link latency;Latency;Status blocks", dtc, roc), kLatency);

    LinkHists& l = linkHists_[link];
    l.latencyUs = h.book1<TH1F>(
        Form("h1_latency_link%d", link),
        Form("Link latency distribution (link %d);Latency [#mus];Entries", link), kLatencyUs);
    l.latencyUs2 = h.book1<TH1F>(
        Form("h1_latency2_link%d", link),
        Form("Link latency distribution (link %d);Latency [#mus];Entries", link), kLatencyUs2);
    l.linkStatus = h.book1<TH1F>(
        Form("h1_linkStatusSummary_link%d", link),
        Form("Link status bit occupancy (link %d);Error flag;Events with flag set", link),
        flagAxis(kNLinkStatusBits));
    l.rocStatus = h.book1<TH1F>(
        Form("h1_rocStatusSummary_link%d", link),
        Form("ROC status bit occupancy (link %d);Error flag;Events with flag set", link),
        flagAxis(kNRocStatusBits));
    l.rocGroup = h.book1<TH1F>(
        Form("h1_rocGroupSummary_link%d", link),
        Form("ROC group bit occupancy (link %d);Group flag;Events with flag set", link),
        flagAxis(kNRocGroupBits));
    l.portFlagBits = h.book1<TH1F>(
        Form("h1_portFlagBits_link%d", link),
        Form("Port flag bit occupancy (link %d);Port (bit);Events with bit set", link),
        flagAxis(kNPortFlags));
    labelX(l.linkStatus, kNLinkStatusBits, kLinkStatusBitNames);
    labelX(l.rocStatus, kNRocStatusBits, kRocStatusBitNames);
    labelX(l.rocGroup, kNRocGroupBits, kRocGroupBitNames);
  }

  h2_rocStatusSummary_ = h.book2<TH2F>(
      "h2_rocStatusSummary", "ROC status bit occupancy;Error flag;DTC#times6 + link ID",
      flagAxis(kNRocStatusBits), kLink);
  h2_rocGroupSummary_ = h.book2<TH2F>(
      "h2_rocGroupSummary", "ROC group bit occupancy;Group flag;DTC#times6 + link ID",
      flagAxis(kNRocGroupBits), kLink);
  h2_portFlagBits_ = h.book2<TH2F>(
      "h2_portFlagBits", "Port flag bit occupancy;Port (bit);DTC#times6 + link ID",
      flagAxis(kNPortFlags), kLink);
  labelX(h2_rocStatusSummary_, kNRocStatusBits, kRocStatusBitNames);
  labelX(h2_rocGroupSummary_, kNRocGroupBits, kRocGroupBitNames);

  g_errorFractionVsSubrun_ = &series().book(
      "g_errorFractionVsSubrun",
      "Fraction of events with a firmware error vs subrun;Subrun;Fraction");
  g_meanLatencyVsSubrun_ = &series().book(
      "g_meanLatencyVsSubrun", "Mean link latency vs subrun;Subrun;Mean latency");
}

// Per-link graphs are online-only and not part of the fixed set, so they are
// booked on a link's first status rather than for every possible link.
DQMSeries& CRVStatusDQM::linkSeries(const char* nameFmt, const char* titleFmt, int link)
{
  return series().book(Form(nameFmt, link), Form(titleFmt, link));
}

DQMSeries& CRVStatusDQM::bitSeries(const char* prefix, const char* bitName, int bit, int link)
{
  return series().book(
      Form("graphs/%s_link%d_%s", prefix, link, bitName),
      Form("Link %d %s (bit %d) vs EWT;Event window tag;Bit value", link, bitName, bit));
}

void CRVStatusDQM::fillLink(const CrvStatus& status, int link)
{
  LinkHists& l = linkHists_[link];
  const uint16_t latency = status.GetLinkLatency();
  const double ewt = static_cast<double>(status.GetEventWindowTag());
  const bool graphs = series().enabled();

  //0xFFFF would read as a 4.2 ms latency here; it has its own histogram
  const bool validLatency = latency != kLatencySentinel;
  const double latencyUs = latency * kLatencyTickToUs;
  if (validLatency) {
    l.latencyUs.Fill(latencyUs);
    l.latencyUs2.Fill(latencyUs);
  }
  if (graphs && validLatency) {
    linkSeries("graphs/g_linkLatencyVsEwt_link%d",
               "Link latency vs EWT (link %d);Event window tag;Latency [#mus]", link)
        .Step(ewt, latencyUs);
  }

  const uint8_t linkStatus = status.GetLinkStatus();
  for (int i = 0; i < kNLinkStatusBits; ++i) {
    if ((linkStatus >> kLinkStatusBits[i]) & 1u) {
      l.linkStatus.Fill(i + 1);
      if (graphs) {
        bitSeries("g_linkStatus", kLinkStatusBitNames[i], kLinkStatusBits[i], link)
            .Add(ewt, 1.);
      }
    }
  }

  if (!status.HasROCHeader()) {
    return;
  }
  const uint32_t ubStatus = status.GetROCHeader().front().GetMicroBunchStatus();
  if (graphs) {
    linkSeries("graphs/g_rocStatus_link%d",
               "ROC status vs EWT (link %d);Event window tag;ROC status", link)
        .Step(ewt, ubStatus);
  }
  for (int i = 0; i < kNRocStatusBits; ++i) {
    if ((ubStatus >> kRocStatusBits[i]) & 1u) {
      l.rocStatus.Fill(i + 1);
      h2_rocStatusSummary_.Fill(i + 1, link);
      if (graphs) {
        bitSeries("g_rocStatus", kRocStatusBitNames[i], kRocStatusBits[i], link).Add(ewt, 1.);
      }
    }
  }
  for (int i = 0; i < kNRocGroupBits; ++i) {
    if ((ubStatus >> kRocGroupBits[i]) & 1u) {
      l.rocGroup.Fill(i + 1);
      h2_rocGroupSummary_.Fill(i + 1, link);
      if (graphs) {
        bitSeries("g_rocStatus", kRocGroupBitNames[i], kRocGroupBits[i], link).Add(ewt, 1.);
      }
    }
  }

  const uint32_t portFlags = ubStatus & 0x00FFFFFFu;
  if (graphs) {
    linkSeries("graphs/g_portFlags_link%d",
               "Port flags vs EWT (link %d);Event window tag;Port flags", link)
        .Step(ewt, portFlags);
  }
  for (int bit = 0; bit < kNPortFlags; ++bit) {
    if ((portFlags >> bit) & 1u) {
      l.portFlagBits.Fill(bit + 1);
      h2_portFlagBits_.Fill(bit + 1, link);
    }
  }
}

void CRVStatusDQM::Fill(const CrvStatusCollection& crvStatus)
{
  if (!booked()) {
    return;
  }
  const bool haveEwt = !crvStatus.empty();
  beginEvent(haveEwt ? std::optional<uint64_t>(crvStatus.front().GetEventWindowTag())
                     : std::nullopt);
  ++nEventsThisSubrun_;
  lastEventRocs_.clear();

  int nHeadersThisEvent = 0;
  bool anyErrorThisEvent = false;

  for (const auto& status : crvStatus) {
    const uint8_t dtcId = status.GetDTCID();
    const uint8_t linkId = status.GetLinkID();
    const uint16_t latency = status.GetLinkLatency();
    const bool indexed = linkInRange(dtcId, linkId);
    const int link = indexed ? globalLink(dtcId, linkId) : -1;
    if (!indexed) {
      diag().Count("linkOutOfRange",
                   Form("status from DTC %d link %d is outside the %d-link axis; the "
                        "per-link histograms, rocCensus and errorBitsVsRoc skip it.",
                        dtcId, linkId, kNLinks));
    }

    //link -1 (off the link axis) lands in underflow, so nothing is silently dropped
    h_statusBlocksByLink_.Fill(link);
    h_linkLatency_.Fill(latency);
    if (latency == kLatencySentinel) {
      h_linkLatencyInvalid_.Fill(link);
      diag().Count("latencySentinel",
                   "link latency 0xFFFF (no measurement): counted per link in "
                   "linkLatencyInvalid, left out of the latency means and per-link views.");
    } else {
      latencySumThisSubrun_ += latency;
      ++latencyNThisSubrun_;
    }
    if (indexed) {
      h_linkLatencyByLink_[link].Fill(latency);
      fillLink(status, link);
    }

    if (!status.HasROCHeader()) {
      continue;
    }
    const auto& roc = status.GetROCHeader().front();

    ++nHeadersThisEvent;
    ++nRocHeadersTotal_;
    seenRocs_.insert({dtcId, linkId});
    if (indexed) {
      h_rocCensus_.Fill(link);
    }

    const std::bitset<kNPortFlags> activeFEBs = roc.GetActiveFEBFlags();
    const uint16_t nActive = static_cast<uint16_t>(activeFEBs.count());
    const uint32_t ubStatus = roc.GetMicroBunchStatus();
    const uint64_t rocEwt = roc.GetEventWindowTag();

    RocSnapshot snap;
    snap.dtcId = dtcId;
    snap.linkId = linkId;
    snap.ewt = rocEwt;
    snap.triggerCount = roc.TriggerCount;
    snap.wordCount = roc.ControllerEventWordCount;
    snap.activeFebCount = nActive;
    snap.linkLatency = latency;
    snap.microBunchStatus = ubStatus;
    lastEventRocs_.push_back(snap);

    nActiveFEBsMin_ = std::min(nActiveFEBsMin_, nActive);
    nActiveFEBsMax_ = std::max(nActiveFEBsMax_, nActive);
    nActiveFEBsSum_ += nActive;
    ++nActiveFEBsSamples_;

    h_activeFebCount_.Fill(nActive);
    h_triggerCount_.Fill(snap.triggerCount);
    h_wordCount_.Fill(snap.wordCount);
    h_ewtMismatch_.Fill(static_cast<double>(static_cast<int64_t>(rocEwt) -
                                            static_cast<int64_t>(status.GetEventWindowTag())));

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
          h_errorBitsVsRoc_.Fill(b, link);
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
  h_eventHasError_.Fill(anyErrorThisEvent ? 1. : 0.);
}

void CRVStatusDQM::Fill(const CrvStatusCollection& crvStatus,
                        const CrvDAQerrorCollection& crvDaqErrors)
{
  Fill(crvStatus);
  if (booked()) {
    fillDaqErrors(crvDaqErrors);
  }
}

void CRVStatusDQM::fillDaqErrors(const CrvDAQerrorCollection& crvDaqErrors)
{
  for (const auto& err : crvDaqErrors) {
    const int code = static_cast<int>(err.GetErrorCode());
    h_daqErrorCode_.Fill(code);
    if (code < 0 || code >= kNDaqErrorCodes) {
      diag().Count("daqErrorCodeOffAxis",
                   Form("CrvDAQerror code %d is past the %d-bin daqErrorCode axis "
                        "(it is in the overflow bin).", code, kNDaqErrorCodes));
    }
  }
  const bool hasError = !crvDaqErrors.empty();
  if (hasError) ++nEventsWithDaqErrors_;
  h_eventHasDaqError_.Fill(hasError ? 1. : 0.);
}

// Filled before the core archives the subrun, so a subrun copy holds its own value.
void CRVStatusDQM::endSubRun()
{
  if (nEventsThisSubrun_ == 0) {
    return;
  }
  const double fraction = static_cast<double>(nEventsWithAnyErrorBitThisSubrun_) /
                          static_cast<double>(nEventsThisSubrun_);
  const double meanLat = latencyNThisSubrun_ > 0
                             ? static_cast<double>(latencySumThisSubrun_) /
                                   static_cast<double>(latencyNThisSubrun_)
                             : 0.0;
  h_errorFractionPerSubrun_.Fill(fraction);
  h_meanLatencyPerSubrun_.Fill(meanLat);
  g_errorFractionVsSubrun_->Add(subrun(), fraction);
  g_meanLatencyVsSubrun_->Add(subrun(), meanLat);

  nEventsThisSubrun_ = 0;
  nEventsWithAnyErrorBitThisSubrun_ = 0;
  latencySumThisSubrun_ = 0;
  latencyNThisSubrun_ = 0;
}

void CRVStatusDQM::resetForNewRun()
{
  nEventsThisSubrun_ = 0;
  nEventsWithAnyErrorBitThisSubrun_ = 0;
  latencySumThisSubrun_ = 0;
  latencyNThisSubrun_ = 0;
}

double CRVStatusDQM::nActiveFEBsMean() const
{
  if (nActiveFEBsSamples_ == 0) {
    return 0.0;
  }
  return static_cast<double>(nActiveFEBsSum_) / static_cast<double>(nActiveFEBsSamples_);
}

std::size_t CRVStatusDQM::errorBitCount(int bitIndex) const
{
  if (bitIndex < 0 || bitIndex >= kNErrorBits) {
    return 0;
  }
  return errorBitCounts_[bitIndex];
}

} // namespace mu2e
