// Standalone CRV ROC-status DQM helper.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/CRVStatusDQM.hh"
#include "Offline/DQMHelpers/inc/DQMSegmentationConfig.hh"

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

// Online link-status layout (otsdaq-mu2e-crv CrvDQM, mu2e/ots_ops).
constexpr std::size_t kNLinkStatusBits = 5;
constexpr uint8_t kLinkStatusBits[kNLinkStatusBits] = {0, 2, 3, 6, 7};
constexpr const char* kLinkStatusBitNames[kNLinkStatusBits] = {
    "ROCTimeout", "SeqNumErr", "CRCErr", "FatalErr", "Error"};

// MicroBunchStatus summary bits [24:31], split into error and group flags.
constexpr std::size_t kNRocStatusBits = 5;
constexpr uint8_t kRocStatusBits[kNRocStatusBits] = {24, 25, 26, 30, 31};
constexpr const char* kRocStatusBitNames[kNRocStatusBits] = {
    "FEBuBMismatch", "FEBBufferIssue", "FEBOverflow", "uBMatchErr", "Truncation"};

constexpr std::size_t kNRocGroupBits = 3;
constexpr uint8_t kRocGroupBits[kNRocGroupBits] = {27, 28, 29};
constexpr const char* kRocGroupBitNames[kNRocGroupBits] = {
    "Group1Issue", "Group2Issue", "Group3Issue"};

constexpr int kNPortFlagBits = 24;

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

CRVStatusDQM::Config toConfig(const CRVStatusDQMFhicl& c)
{
  CRVStatusDQM::Config out;
  out.nBinsLatency = std::max(c.nBinsLatency(), 1);
  out.maxLinkLatency = c.maxLinkLatency();
  out.nBinsTriggerCount = std::max(c.nBinsTriggerCount(), 1);
  out.maxTriggerCount = c.maxTriggerCount();
  out.nBinsWordCount = std::max(c.nBinsWordCount(), 1);
  out.maxWordCount = c.maxWordCount();
  out.nBinsEwtMismatch = std::max(c.nBinsEwtMismatch(), 1);
  out.maxEwtMismatch = c.maxEwtMismatch();
  out.nBinsErrorsPerSubrun = std::max(c.nBinsErrorsPerSubrun(), 1);
  out.maxErrorsPerSubrun = c.maxErrorsPerSubrun();
  out.fillLivePlots = c.fillLivePlots();
  out.fillLinkPlots = c.fillLinkPlots();
  out.fillLinkGraphs = c.fillLinkGraphs();
  out.maxLinkGraphPoints = static_cast<std::size_t>(std::max(c.maxLinkGraphPoints(), 1));
  out.segmentation = parseSegmentation(c.segmentation);
  return out;
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
      nDaqErrorCodes(),
      -0.5,
      nDaqErrorCodes() - 0.5);
  h_daqErrorCode_.ForEach([](TH1F* h) {
    for (const auto& [code, name] : CrvDAQerrorCodeDetail::names()) {
      h->GetXaxis()->SetBinLabel(static_cast<int>(code) + 1, name.c_str());
    }
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

  if (config_.fillLinkPlots) {
    auto book2Labelled = [this](const char* name, const char* title, int nx,
                                const char* const* labels) {
      DQMHist2<TH2F> h = segments_.book2<TH2F>(name, title, nx, 0.5, nx + 0.5,
                                               kNLinkBins, -0.5, kNLinkBins - 0.5);
      if (labels != nullptr) {
        h.ForEach([nx, labels](TH2F* p) {
          for (int i = 0; i < nx; ++i) {
            p->GetXaxis()->SetBinLabel(i + 1, labels[i]);
          }
        });
      }
      return h;
    };
    h2_rocStatusSummary_ =
        book2Labelled("h2_rocStatusSummary", "ROC status bit occupancy;Error flag;Link",
                      kNRocStatusBits, kRocStatusBitNames);
    h2_rocGroupSummary_ =
        book2Labelled("h2_rocGroupSummary", "ROC group bit occupancy;Group flag;Link",
                      kNRocGroupBits, kRocGroupBitNames);
    h2_portFlagBits_ =
        book2Labelled("h2_portFlagBits", "Port flag bit occupancy;Port (bit);Link",
                      kNPortFlagBits, nullptr);
    if (config_.fillLinkGraphs) {
      graphDir_ = dir.mkdir("graphs");
    }
  }

  booked_ = true;
}

DQMHist1<TH1F> CRVStatusDQM::linkHist(std::map<uint8_t, DQMHist1<TH1F>>& hists,
                                      uint8_t linkId, const char* name,
                                      const char* title, int nBins, double lo,
                                      double hi, const char* const* labels)
{
  auto it = hists.find(linkId);
  if (it != hists.end()) {
    return it->second;
  }
  const std::string hname = Form(name, static_cast<int>(linkId));
  const std::string htitle = Form(title, static_cast<int>(linkId));
  DQMHist1<TH1F> h = segments_.book1<TH1F>(hname, htitle, nBins, lo, hi);
  if (labels != nullptr) {
    h.ForEach([nBins, labels](TH1F* p) {
      for (int i = 0; i < nBins; ++i) {
        p->GetXaxis()->SetBinLabel(i + 1, labels[i]);
      }
    });
  }
  hists[linkId] = h;
  return h;
}

TGraph* CRVStatusDQM::makeLinkGraph(const std::string& name, const std::string& title)
{
  TGraph* g = graphDir_->make<TGraph>();
  g->SetName(name.c_str());
  g->SetTitle(title.c_str());
  linkGraphs_.push_back(g);
  return g;
}

void CRVStatusDQM::stepPoint(std::map<uint8_t, TGraph*>& graphs,
                             std::map<uint8_t, uint32_t>& last, uint8_t linkId,
                             uint32_t value, uint64_t ewt, double scale,
                             const char* nameFmt, const char* titleFmt)
{
  const double x = static_cast<double>(ewt);
  const double y = value * scale;
  auto it = graphs.find(linkId);
  if (it == graphs.end()) {
    const std::string name = Form(nameFmt, static_cast<int>(linkId));
    const std::string title = Form(titleFmt, static_cast<int>(linkId));
    TGraph* g = makeLinkGraph(name, title);
    g->SetPoint(0, x, y);
    graphs[linkId] = g;
    last[linkId] = value;
    return;
  }
  TGraph* g = it->second;
  auto lastIt = last.find(linkId);
  if (lastIt == last.end()) {
    // First status after ResetForNewRun: re-seed the emptied graph.
    g->SetPoint(g->GetN(), x, y);
    last[linkId] = value;
    return;
  }
  if (value != lastIt->second) {
    g->SetPoint(g->GetN(), x, lastIt->second * scale);
    g->SetPoint(g->GetN(), x, y);
    lastIt->second = value;
    while (static_cast<std::size_t>(g->GetN()) > config_.maxLinkGraphPoints) {
      g->RemovePoint(0);
    }
    return;
  }
  // Unchanged: slide the last point forward so the step reaches the present.
  if (g->GetN() > 0) {
    g->SetPoint(g->GetN() - 1, x, y);
  }
}

void CRVStatusDQM::bookBitGraphs(std::map<std::pair<uint8_t, uint8_t>, TGraph*>& graphs,
                                 uint8_t linkId, const uint8_t* bits,
                                 const char* const* names, std::size_t n,
                                 const char* prefix)
{
  for (std::size_t i = 0; i < n; ++i) {
    const auto key = std::make_pair(linkId, bits[i]);
    if (graphs.find(key) != graphs.end()) {
      continue;
    }
    const std::string name =
        Form("%s_link%d_%s", prefix, static_cast<int>(linkId), names[i]);
    const std::string title = Form("Link %d %s (bit %d) vs EWT;Event window tag;Bit value",
                                   static_cast<int>(linkId), names[i],
                                   static_cast<int>(bits[i]));
    graphs[key] = makeLinkGraph(name, title);
  }
}

void CRVStatusDQM::addBitPoint(TGraph* g, uint64_t ewt) const
{
  g->SetPoint(g->GetN(), static_cast<double>(ewt), 1.0);
  while (static_cast<std::size_t>(g->GetN()) > config_.maxLinkGraphPoints) {
    g->RemovePoint(0);
  }
}

void CRVStatusDQM::fillLink(const CrvStatus& status)
{
  const uint8_t linkId = status.GetLinkID();
  const uint16_t latency = status.GetLinkLatency();
  const uint64_t ewt = status.GetEventWindowTag();
  const bool graphs = config_.fillLinkGraphs && graphDir_.has_value();

  const float latencyUs = static_cast<float>(latency * kLatencyTickToUs);
  linkHist(h_latencyUs_, linkId, "h1_latency_link%d",
           "Link latency distribution (link %d);Latency [#mus];Entries", 300, 0.0, 150.0)
      .Fill(latencyUs);
  linkHist(h_latencyUsWide_, linkId, "h1_latency2_link%d",
           "Link latency distribution (link %d);Latency [#mus];Entries", 2000, 0.0,
           10000.0)
      .Fill(latencyUs);
  if (graphs) {
    stepPoint(g_linkLatency_, lastLinkLatency_, linkId, latency, ewt, kLatencyTickToUs,
              "g_linkLatencyVsEwt_link%d",
              "Link latency vs EWT (link %d);Event window tag;Latency [#mus]");
  }

  const uint8_t linkStatus = status.GetLinkStatus();
  const DQMHist1<TH1F> hLinkStatus = linkHist(
      h_linkStatusSummary_, linkId, "h1_linkStatusSummary_link%d",
      "Link status bit occupancy (link %d);Error flag;Events with flag set",
      kNLinkStatusBits, 0.5, kNLinkStatusBits + 0.5, kLinkStatusBitNames);
  if (graphs) {
    bookBitGraphs(g_linkStatusBit_, linkId, kLinkStatusBits, kLinkStatusBitNames,
                  kNLinkStatusBits, "g_linkStatus");
  }
  for (std::size_t i = 0; i < kNLinkStatusBits; ++i) {
    if ((linkStatus >> kLinkStatusBits[i]) & 1u) {
      hLinkStatus.Fill(i + 1);
      if (graphs) {
        addBitPoint(g_linkStatusBit_[{linkId, kLinkStatusBits[i]}], ewt);
      }
    }
  }

  if (!status.HasROCHeader()) {
    return;
  }
  const uint32_t ubStatus = status.GetROCHeader().front().GetMicroBunchStatus();

  if (graphs) {
    stepPoint(g_rocStatus_, lastRocStatus_, linkId, ubStatus, ewt, 1.0,
              "g_rocStatus_link%d", "ROC status vs EWT (link %d);Event window tag;ROC status");
  }
  const DQMHist1<TH1F> hRocStatus = linkHist(
      h_rocStatusSummary_, linkId, "h1_rocStatusSummary_link%d",
      "ROC status bit occupancy (link %d);Error flag;Events with flag set",
      kNRocStatusBits, 0.5, kNRocStatusBits + 0.5, kRocStatusBitNames);
  const DQMHist1<TH1F> hRocGroup = linkHist(
      h_rocGroupSummary_, linkId, "h1_rocGroupSummary_link%d",
      "ROC group bit occupancy (link %d);Group flag;Events with flag set",
      kNRocGroupBits, 0.5, kNRocGroupBits + 0.5, kRocGroupBitNames);
  if (graphs) {
    bookBitGraphs(g_rocStatusBit_, linkId, kRocStatusBits, kRocStatusBitNames,
                  kNRocStatusBits, "g_rocStatus");
    bookBitGraphs(g_rocStatusBit_, linkId, kRocGroupBits, kRocGroupBitNames,
                  kNRocGroupBits, "g_rocStatus");
  }
  for (std::size_t i = 0; i < kNRocStatusBits; ++i) {
    if ((ubStatus >> kRocStatusBits[i]) & 1u) {
      hRocStatus.Fill(i + 1);
      h2_rocStatusSummary_.Fill(i + 1, linkId);
      if (graphs) {
        addBitPoint(g_rocStatusBit_[{linkId, kRocStatusBits[i]}], ewt);
      }
    }
  }
  for (std::size_t i = 0; i < kNRocGroupBits; ++i) {
    if ((ubStatus >> kRocGroupBits[i]) & 1u) {
      hRocGroup.Fill(i + 1);
      h2_rocGroupSummary_.Fill(i + 1, linkId);
      if (graphs) {
        addBitPoint(g_rocStatusBit_[{linkId, kRocGroupBits[i]}], ewt);
      }
    }
  }

  const uint32_t portFlags = ubStatus & 0x00FFFFFFu;
  const DQMHist1<TH1F> hPortFlags = linkHist(
      h_portFlagBits_, linkId, "h1_portFlagBits_link%d",
      "Port flag bit occupancy (link %d);Port (bit);Events with bit set",
      kNPortFlagBits, 0.5, kNPortFlagBits + 0.5);
  if (graphs) {
    stepPoint(g_portFlags_, lastPortFlags_, linkId, portFlags, ewt, 1.0,
              "g_portFlags_link%d", "Port flags vs EWT (link %d);Event window tag;Port flags");
  }
  for (int bit = 0; bit < kNPortFlagBits; ++bit) {
    if ((portFlags >> bit) & 1u) {
      hPortFlags.Fill(bit + 1);
      h2_portFlagBits_.Fill(bit + 1, linkId);
    }
  }
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
    if (config_.fillLinkPlots && booked_) {
      fillLink(status);
    }
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

int CRVStatusDQM::nDaqErrorCodes()
{
  const auto& names = CrvDAQerrorCodeDetail::names();
  return names.empty() ? 1 : static_cast<int>(names.rbegin()->first) + 1;
}

void CRVStatusDQM::fillDaqErrors(const CrvDAQerrorCollection& crvDaqErrors)
{
  bool countedEvent = false;
  for (const auto& err : crvDaqErrors) {
    const int code = static_cast<int>(err.GetErrorCode());
    if (code >= 0 && code < nDaqErrorCodes()) {
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

void CRVStatusDQM::persistGraph(TGraph* g, art::TFileDirectory* dir)
{
  art::TFileDirectory* d = dir != nullptr ? dir : (dir_ ? &*dir_ : nullptr);
  if (d == nullptr || g == nullptr) {
    return;
  }
  if (g->GetN() <= 0) {
    d->makeAndRegister<TGraph>(g->GetName(), g->GetTitle());
    return;
  }
  d->makeAndRegister<TGraph>(
      g->GetName(), g->GetTitle(), g->GetN(), g->GetX(), g->GetY());
}

void CRVStatusDQM::ResetForNewRun()
{
  segments_.ResetContents();
  for (TGraph* g : linkGraphs_) {
    g->Set(0);
  }
  for (TGraph* g : {g_errorsVsSubrun_, g_meanLatencyVsSubrun_}) {
    if (g != nullptr) {
      g->Set(0);
    }
  }
  lastLinkLatency_.clear();
  lastRocStatus_.clear();
  lastPortFlags_.clear();
  nEventsThisSubrun_ = 0;
  nEventsWithAnyErrorBitThisSubrun_ = 0;
  latencySumThisSubrun_ = 0;
  latencyNThisSubrun_ = 0;
}

void CRVStatusDQM::WriteGraphs()
{
  if (!booked_ || !dir_) {
    return;
  }
  segments_.Finalize();
  persistGraph(g_errorsVsSubrun_);
  persistGraph(g_meanLatencyVsSubrun_);
  if (graphDir_) {
    for (TGraph* g : linkGraphs_) {
      persistGraph(g, &*graphDir_);
    }
  }
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
