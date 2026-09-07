// Standalone CRV digi DQM helper.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/CRVDigiDQM.hh"
#include "Offline/DQMHelpers/inc/DQMSegmentationConfig.hh"
#include "Offline/DQMHelpers/inc/CRVCFTime.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "TString.h"
#include "TH1.h"

namespace mu2e {

int CRVDigiDQM::globalFebId(uint8_t roc, uint8_t feb)
{
  return (static_cast<int>(roc) - 1) * kNFebSlotsPerROC + feb;
}

int CRVDigiDQM::globalChannelId(uint8_t roc, uint8_t feb, uint8_t febChannel)
{
  return globalFebId(roc, feb) * static_cast<int>(CRVId::nChanPerFEB) +
         febChannel;
}

CRVDigiDQM::Config toConfig(const CRVDigiDQMFhicl& c)
{
  CRVDigiDQM::Config out;
  out.nBinsDigisPerEvt = c.nBinsDigisPerEvt();
  out.maxDigisPerEvt = c.maxDigisPerEvt();
  out.nBinsPeakAdc = c.nBinsPeakAdc();
  out.maxPeakAdc = c.maxPeakAdc();
  out.nBinsTdc = c.nBinsTdc();
  out.maxTdc = c.maxTdc();
  out.cfFraction = c.cfFraction();
  out.dtBinSize = c.dtBinSize();
  out.dtRange = c.dtRange();
  out.dtVsFebBinSize = c.dtVsFebBinSize();
  out.dtVsFebRange = c.dtVsFebRange();
  out.minAmplitude = c.minAmplitude();
  out.avgBlockSize = static_cast<std::size_t>(std::max(c.avgBlockSize(), 1));
  out.avgGraphPoints = static_cast<std::size_t>(std::max(c.avgGraphPoints(), 1));
  out.channelsWindowEwts =
      static_cast<std::size_t>(std::max(c.channelsWindowEwts(), 1));
  out.fillInclusive = c.fillInclusive();
  out.fillCrvIdRates = c.fillCrvIdRates();
  out.kppReadout = c.kppReadout();
  out.fillLivePlots = c.fillLivePlots();
  out.segmentation = parseSegmentation(c.segmentation);
  return out;
}

CRVDigiDQM::CRVDigiDQM(const Config& config) :
    config_(config), nDigisOffline_(CRVId::nChannels, 0)
{
  if (config_.channelsWindowEwts == 0) {
    config_.channelsWindowEwts = 1;
  }
  if (config_.avgBlockSize == 0) {
    config_.avgBlockSize = 1;
  }
  if (config_.dtBinSize > 0) {
    nBinsDt_ = static_cast<int>(2.0 * config_.dtRange / config_.dtBinSize);
  }
  if (nBinsDt_ < 1) {
    nBinsDt_ = 1;
  }
  if (config_.dtVsFebBinSize > 0) {
    nBinsDtVsFeb_ =
        static_cast<int>(2.0 * config_.dtVsFebRange / config_.dtVsFebBinSize);
  }
  if (nBinsDtVsFeb_ < 1) {
    nBinsDtVsFeb_ = 1;
  }

  // channelsWindowEwts is the helper's own window default; a rule that names a
  // span of its own wins, so the two knobs cannot disagree silently.
  DQMSegmentation::Config segmentation = config_.segmentation;
  for (auto& rule : segmentation.rules) {
    if (rule.window.enabled && !rule.window.spanSet) {
      rule.window.span = config_.channelsWindowEwts;
    }
  }
  segments_.SetConfig(segmentation);
}

void CRVDigiDQM::Book(art::TFileDirectory dir)
{
  dir_ = dir;
  segments_.Book(dir);

  h2_dtVsFeb_ = segments_.book2<TH2F>(
      "dtVsFeb",
      "First-hit time vs median of other FEBs;Global FEB ID;#Deltat [ns]",
      nFebIdBins(),
      -0.5,
      nFebIdBins() - 0.5,
      nBinsDtVsFeb_,
      -config_.dtVsFebRange,
      config_.dtVsFebRange);

  // Per-FEB desync counters: a slipped FEB stands above its neighbours.
  h_dtOutOfRangePerFeb_ = segments_.book1<TH1F>(
      "dtOutOfRangePerFeb",
      Form("Events with |#Deltat| > %.0f ns;Global FEB ID;Events",
           config_.dtVsFebRange),
      nFebIdBins(), -0.5, nFebIdBins() - 0.5);

  h1_digisPerEvt_ = segments_.book1<TH1F>("h1_digisPerEvt",
                                          "Hits / event;Hits / event;Events",
                                          config_.nBinsDigisPerEvt,
                                          -0.5,
                                          config_.maxDigisPerEvt + 0.5);
  h1_digisPerEvt_.ForEach([](TH1F* h) { h->SetMinimum(0.5); });

  h_nEvents_ =
      segments_.book1<TH1F>("nEvents", "Events processed;;Events", 1, 0.5, 1.5);

  h1_peakAdc_ = segments_.book1<TH1F>("h1_peakAdc",
                                      "Max sample ADC;Max sample ADC;Hits",
                                      config_.nBinsPeakAdc,
                                      0,
                                      config_.maxPeakAdc);

  h1_tdc_ = segments_.book1<TH1F>(
      "h1_tdc",
      "Start timestamp of digi in units of 12.5ns;Start timestamp of digi;Digis",
      config_.nBinsTdc,
      0,
      config_.maxTdc);

  // KPP-only: these axes cover 33 FEB slots, so on full CRV they would be all
  // overflow. crvDigisPerChannel / crvDigiRates cover the full detector instead.
  if (config_.kppReadout) {
    h1_channels_ = segments_.book1<TH1F>("h1_channels",
                                         "Channel occupancy;Global channel ID;Hits",
                                         kNGlobalChannelBins,
                                         -0.5,
                                         kNGlobalChannelBins - 0.5);
    h1_channels_.ForEach([](TH1F* h) { h->SetMinimum(0.5); });

    h2_channels_ = segments_.book2<TH2F>(
        "h2_channels",
        "FEB vs channel hit map;Channel;FEB",
        static_cast<int>(CRVId::nChanPerFEB),
        0.5,
        static_cast<double>(CRVId::nChanPerFEB) + 0.5,
        kNGlobalFebBins,
        0.5,
        kNGlobalFebBins + 0.5);
  }

  if (config_.fillLivePlots) {
    g_digisVsEwt_ = dir.make<TGraph>();
    g_digisVsEwt_->SetName("g_digisVsEwt");
    g_digisVsEwt_->SetTitle(Form("Hits in last %zu EWTs (%zu points);"
                                 "Event window tag;Hits",
                                 kEwtWindow,
                                 kGraphPoints));
    g_digisVsEwt_->SetPoint(0, 0, 0);
    g_digisVsEwt_->SetPoint(1, 1, 1);

    g_digisAvgVsEwt_ = dir.make<TGraph>();
    g_digisAvgVsEwt_->SetName("g_digisAvgVsEwt");
    g_digisAvgVsEwt_->SetTitle(Form("Mean hits per event (averaged over %zu events);"
                                    "Event window tag; ",
                                    config_.avgBlockSize));
    g_digisAvgVsEwt_->SetPoint(0, 0, 0);
    g_digisAvgVsEwt_->SetPoint(1, 1, 1);
  }

  if (config_.fillInclusive) {
    hBarId_ = segments_.book1<TH1D>("BarId", "Bar ID", 200, -0.5,
                                    static_cast<double>(CRVId::nBars) - 0.5);
    hSiPM_ = segments_.book1<TH1D>("SiPM", "SiPM", 4, -0.5, 3.5);
    hADC_ = segments_.book1<TH1D>("ADC", "ADC in waveform", 100, 0.0, 3000.0);
  }

  if (config_.fillCrvIdRates) {
    h_crvDigiRatesROC_.assign(CRVId::nROC, DQMHist1<TH1F>());
    for (std::size_t roc = 1; roc <= CRVId::nROC; ++roc) {
      h_crvDigiRatesROC_[roc - 1] = segments_.book1<TH1F>(
          Form("crvDigiRates_ROC%zu", roc),
          Form("crvDigiRates_ROC%zu;Online channel in ROC;Digis", roc),
          static_cast<int>(CRVId::nFEBPerROC * CRVId::nChanPerFEB),
          0,
          static_cast<double>(CRVId::nFEBPerROC * CRVId::nChanPerFEB));
    }
    h_crvDigiRates_ = segments_.book2<TH2F>(
        "crvDigiRates",
        "crvDigiRates:FEBchannel:FEB;FEB channel;FEB port",
        static_cast<int>(CRVId::nChanPerFEB),
        0,
        static_cast<double>(CRVId::nChanPerFEB),
        static_cast<int>(CRVId::nROC * CRVId::nFEBPerROC),
        0,
        static_cast<double>(CRVId::nROC * CRVId::nFEBPerROC));
    h_crvDigisPerChannel_ = segments_.book1<TH1F>(
        "crvDigisPerChannel",
        "Digis vs offline channel;Offline channel (bar#times4+SiPM);Digis",
        static_cast<int>(CRVId::nChannels),
        -0.5,
        static_cast<double>(CRVId::nChannels) - 0.5);
  }

  booked_ = true;
}

void CRVDigiDQM::Fill(const CrvDigiCollection& crvDigis,
                      const CrvStatusCollection& crvStatus)
{
  if (!booked_) {
    return;
  }

  ++nEvents_;

  // Rotate the segmentation ring between events, before anything is filled, so
  // no entry can straddle two segments.
  const bool haveEwt = !crvStatus.empty();
  const uint64_t ewt = haveEwt ? crvStatus.front().GetEventWindowTag() : 0;
  segments_.Advance(nEvents_,
                    haveEwt ? std::optional<uint64_t>(ewt) : std::nullopt);

  h_nEvents_.Fill(1.f);

  const int nDigis = static_cast<int>(crvDigis.size());
  std::map<int, std::map<uint8_t, std::vector<FpgaHit>>> hitTimes;

  for (const auto& digi : crvDigis) {
    const uint8_t roc = digi.GetROC();
    const uint8_t feb = digi.GetFEB();
    const uint8_t febChannel = digi.GetFEBchannel();

    const int febId = globalFebId(roc, feb);
    if (febId > maxFebIdSeen_) {
      maxFebIdSeen_ = febId;
    }
    if (febId < 0 || febId >= nFebIdBins()) {
      ++nFebIdOutOfAxis_;
      // Once per job: the axis size comes from Config::kppReadout, so this is
      // a fixed configuration error, not a per-event condition.
      if (!warnedOffAxisFeb_) {
        warnedOffAxisFeb_ = true;
        mf::LogWarning("CRVDigiDQM")
            << "digi from ROC " << static_cast<int>(roc)
            << " FEB " << static_cast<int>(feb) << " gives globalFebId " << febId
            << ", outside the " << nFebIdBins() << "-bin FEB axis. It and any "
            << "others like it are hidden in the overflow bin of dtVsFeb"
            << (config_.kppReadout ? " and h2_channels" : "") << ". "
            << (config_.kppReadout
                    ? "KPP is ROC 1-2 -- check CrvDigisFromArtdaqFragmentsFEBII "
                      "useROC4asROC2, or set kppReadout=false for a larger CRV."
                    : "Raise kNFebIdBins.")
            << " Reported once per job.";
      }
    }

    if (config_.kppReadout) {
      h1_channels_.Fill(globalChannelId(roc, feb, febChannel));
      h2_channels_.Fill(febChannel + 1, febId);
    }

    const auto& adcs = digi.GetADCs();
    if (!adcs.empty()) {
      const int16_t maxSample = *std::max_element(adcs.begin(), adcs.end());
      h1_peakAdc_.Fill(maxSample);
    }

    h1_tdc_.Fill(digi.GetStartTDC());

    const int barIndex = digi.GetScintillatorBarIndex().asUint();
    const int sipm = digi.GetSiPMNumber();
    if (sipm >= 0) {
      const std::size_t offlineChannel =
          static_cast<std::size_t>(barIndex) * CRVId::nChanPerBar +
          static_cast<std::size_t>(sipm);
      if (offlineChannel < nDigisOffline_.size()) {
        ++nDigisOffline_[offlineChannel];
        h_crvDigisPerChannel_.Fill(static_cast<float>(offlineChannel));
      }
    }

    if (config_.fillCrvIdRates) {
      const int rocId = static_cast<int>(roc);
      const int febIdRaw = static_cast<int>(feb);
      const int febCh = static_cast<int>(febChannel);
      if (rocId >= 1 && static_cast<std::size_t>(rocId) <= CRVId::nROC &&
          febIdRaw >= 1 && static_cast<std::size_t>(febIdRaw) <= CRVId::nFEBPerROC &&
          febCh >= 0 && static_cast<std::size_t>(febCh) < CRVId::nChanPerFEB) {
        const int rocChannel =
            (febIdRaw - 1) * static_cast<int>(CRVId::nChanPerFEB) + febCh;
        if (!h_crvDigiRatesROC_.empty()) {
          h_crvDigiRatesROC_[static_cast<std::size_t>(rocId) - 1].Fill(rocChannel);
        }
        const int portIndex =
            (rocId - 1) * static_cast<int>(CRVId::nFEBPerROC) + febIdRaw - 1;
        h_crvDigiRates_.Fill(febCh, portIndex);
      } else {
        ++nCrvIdOutOfRange_;
        if (!warnedCrvIdOutOfRange_) {
          warnedCrvIdOutOfRange_ = true;
          mf::LogWarning("CRVDigiDQM")
              << "digi from ROC " << rocId << " FEB " << febIdRaw
              << " channel " << febCh << " is outside CRVId range (ROC 1-"
              << CRVId::nROC << ", FEB 1-" << CRVId::nFEBPerROC
              << ", channel 0-" << (CRVId::nChanPerFEB - 1)
              << "). It is omitted from crvDigiRates*; occupancy and timing "
              << "are still filled. Reported once per job.";
        }
      }
    }

    if (config_.fillInclusive) {
      hBarId_.Fill(digi.GetScintillatorBarIndex().asInt());
      hSiPM_.Fill(digi.GetSiPMNumber());
      for (auto a : adcs) {
        hADC_.Fill(a);
      }
    }

    const CFResult cf = cfTime(adcs, config_.cfFraction, config_.minAmplitude);
    if (cf.valid) {
      const double absTime_ns =
          cf.time_ns + digi.GetStartTDC() * kDigitizationPeriodNs;
      const uint8_t fpga =
          static_cast<uint8_t>(febChannel / CRVId::nChanPerFPGA);
      hitTimes[febId][fpga].push_back({absTime_ns, febChannel});
    }

    activeROCs_.insert(roc);
    activeFEBs_.insert(febId);
    rocFEBMap_[roc].insert(feb);
  }

  nDigis_ += static_cast<std::size_t>(nDigis);
  h1_digisPerEvt_.Fill(nDigis);

  fillTiming(hitTimes);

  if (haveEwt && config_.fillLivePlots) {
    fillEwtSeries(ewt, nDigis);
    fillMicroBunchStatus(crvStatus);
  }
}

void CRVDigiDQM::BeginSubRun(int run, int subrun)
{
  segments_.BeginSubRun(run, subrun);
}

void CRVDigiDQM::EndSubRun()
{
  segments_.EndSubRun();
}

void CRVDigiDQM::fillEwtSeries(uint64_t ewt, int nDigis)
{
  if (g_digisVsEwt_ == nullptr || g_digisAvgVsEwt_ == nullptr) {
    return;
  }
  if (avgBlockCount_ == 0) {
    avgBlockFirstEwt_ = ewt;
  }
  avgBlockSum_ += nDigis;
  ++avgBlockCount_;
  if (avgBlockCount_ >= config_.avgBlockSize) {
    const double midEwt =
        0.5 * (static_cast<double>(avgBlockFirstEwt_) + static_cast<double>(ewt));
    const double mean = static_cast<double>(avgBlockSum_) /
                        static_cast<double>(avgBlockCount_);

    if (!avgSeedsCleared_) {
      g_digisAvgVsEwt_->Set(0);
      avgSeedsCleared_ = true;
    }
    g_digisAvgVsEwt_->SetPoint(g_digisAvgVsEwt_->GetN(), midEwt, mean);

    while (static_cast<std::size_t>(g_digisAvgVsEwt_->GetN()) >
           config_.avgGraphPoints) {
      g_digisAvgVsEwt_->RemovePoint(0);
    }

    avgBlockSum_ = 0;
    avgBlockCount_ = 0;
  }

  ewtWindow_.emplace_back(ewt, nDigis);
  ewtWindowSum_ += nDigis;
  while (ewtWindow_.size() > kEwtWindow) {
    ewtWindowSum_ -= ewtWindow_.front().second;
    ewtWindow_.pop_front();
  }
  if (ewtWindow_.size() == 1 && g_digisVsEwt_->GetN() == 2) {
    g_digisVsEwt_->Set(0);
  }
  g_digisVsEwt_->SetPoint(g_digisVsEwt_->GetN(),
                          static_cast<double>(ewt),
                          static_cast<double>(ewtWindowSum_));

  while (static_cast<std::size_t>(g_digisVsEwt_->GetN()) > kGraphPoints) {
    g_digisVsEwt_->RemovePoint(0);
  }

  const double currentEwt = static_cast<double>(ewt);
  double xLo = std::max(0.0, currentEwt - kEwtXRange);
  double xHi = currentEwt;
  if (xHi <= xLo) {
    xHi = xLo + 1.0;
  }
  if (TH1* frame = g_digisVsEwt_->GetHistogram()) {
    frame->GetXaxis()->SetLimits(xLo, xHi);
  }

  if (g_digisAvgVsEwt_->GetN() > 0) {
    double* ax = g_digisAvgVsEwt_->GetX();
    const int nAvg = g_digisAvgVsEwt_->GetN();
    double aLo = *std::min_element(ax, ax + nAvg);
    double aHi = *std::max_element(ax, ax + nAvg);
    if (aHi <= aLo) {
      aHi = aLo + 1.0;
    }
    if (TH1* frame = g_digisAvgVsEwt_->GetHistogram()) {
      frame->GetXaxis()->SetLimits(aLo, aHi);
    }
  }
}

// Median of sorted values with the entry at sorted position p removed.
namespace {
double medianExcluding(const std::vector<double>& sorted, std::size_t p)
{
  const std::size_t m = sorted.size() - 1;
  auto at = [&](std::size_t i) { return sorted[i < p ? i : i + 1]; };
  if (m % 2 == 1) {
    return at(m / 2);
  }
  return 0.5 * (at(m / 2 - 1) + at(m / 2));
}
} // namespace

void CRVDigiDQM::fillTiming(
    const std::map<int, std::map<uint8_t, std::vector<FpgaHit>>>& hitTimes)
{
  std::vector<std::pair<int, double>> febFirstHit;
  for (const auto& [febId, fpgaMap] : hitTimes) {
    // TODO: the median hit time may be a more noise-robust FEB reference than
    // the first hit, which an early dark-noise pulse can hijack.
    double earliest = std::numeric_limits<double>::max();
    for (const auto& [fpga, hits] : fpgaMap) {
      for (const auto& hit : hits) {
        if (hit.time_ns < earliest) {
          earliest = hit.time_ns;
        }
      }
    }
    febFirstHit.push_back({febId, earliest});
  }

  // Leave-one-out reference isolates the slipped FEB instead of smearing its
  // partners, which a common reference would do at low FEB multiplicity.
  if (febFirstHit.size() >= 2) {
    std::vector<double> sorted;
    sorted.reserve(febFirstHit.size());
    for (const auto& entry : febFirstHit) {
      sorted.push_back(entry.second);
    }
    std::sort(sorted.begin(), sorted.end());

    for (const auto& entry : febFirstHit) {
      const std::size_t p = static_cast<std::size_t>(
          std::lower_bound(sorted.begin(), sorted.end(), entry.second) -
          sorted.begin());
      const double dt = entry.second - medianExcluding(sorted, p);
      h2_dtVsFeb_.Fill(entry.first, dt);

      // Measured on the raw dt, so an arbitrarily large slip is still counted
      // even when the y-axis would bury it in the overflow bin.
      const double absDt = std::abs(dt);
      if (absDt > maxAbsDtSeen_) {
        maxAbsDtSeen_ = absDt;
      }
      if (absDt > config_.dtVsFebRange) {
        ++nDtOutOfRange_;
        const int febId = entry.first;
        if (febId >= 0 && febId < nFebIdBins()) {
          h_dtOutOfRangePerFeb_.Fill(febId);
        }
      }
    }
  }

  if (!segments_.booked()) {
    return;
  }

  for (const auto& [febId, fpgaMap] : hitTimes) {
    for (auto itA = fpgaMap.begin(); itA != fpgaMap.end(); ++itA) {
      for (auto itB = itA; itB != fpgaMap.end(); ++itB) {
        const uint8_t fpgaA = itA->first;
        const uint8_t fpgaB = itB->first;
        const auto& hitsA = itA->second;
        const auto& hitsB = itB->second;

        const uint8_t pairCode =
            static_cast<uint8_t>(fpgaA * CRVId::nFPGAPerFEB + fpgaB);
        const auto key = std::make_pair(febId, pairCode);

        if (h1_dtFpgaPairs_.find(key) == h1_dtFpgaPairs_.end()) {
          // Booked on first use, so a KPP-sized readout pays only for the FEB
          // pairs it actually has.
          const std::string name =
              Form("timing_fpga/dt_feb%02d_fpga%d_fpga%d", febId, fpgaA, fpgaB);
          const std::string title = Form(
              "#Deltat FEB %02d FPGA %d - FPGA %d;#Deltat [ns];Entries",
              febId,
              fpgaA,
              fpgaB);
          h1_dtFpgaPairs_[key] = segments_.book1<TH1F>(name,
                                                       title,
                                                       nBinsDt_,
                                                       -config_.dtRange,
                                                       config_.dtRange);
        }

        const DQMHist1<TH1F>& h = h1_dtFpgaPairs_[key];
        if (fpgaA == fpgaB) {
          for (std::size_t ia = 0; ia < hitsA.size(); ++ia) {
            for (std::size_t ib = ia + 1; ib < hitsA.size(); ++ib) {
              if (hitsA[ia].channel != hitsA[ib].channel) {
                h.Fill(hitsA[ib].time_ns - hitsA[ia].time_ns);
              }
            }
          }
        } else {
          for (const auto& hA : hitsA) {
            for (const auto& hB : hitsB) {
              h.Fill(hB.time_ns - hA.time_ns);
            }
          }
        }
      }
    }
  }
}

void CRVDigiDQM::fillMicroBunchStatus(const CrvStatusCollection& crvStatus)
{
  if (!dir_ || !config_.fillLivePlots) {
    return;
  }

  for (const auto& status : crvStatus) {
    if (!status.HasROCHeader()) {
      continue;
    }

    const uint8_t linkID = status.GetLinkID();
    const uint32_t ubStatus = status.GetMicroBunchStatus();
    const uint64_t ewt = status.GetEventWindowTag();

    if (g_ubStatusVsEwt_.find(linkID) == g_ubStatusVsEwt_.end()) {
      const std::string name = Form("g_ubStatusVsEwt_link%d", linkID);
      const std::string title = Form("MicroBunchStatus vs EWT (link %d);"
                                     "Event window tag;MicroBunchStatus",
                                     linkID);
      TGraph* g = dir_->make<TGraph>();
      g->SetName(name.c_str());
      g->SetTitle(title.c_str());
      g->SetPoint(0, static_cast<double>(ewt), static_cast<double>(ubStatus));
      g_ubStatusVsEwt_[linkID] = g;
      lastMicroBunchStatus_[linkID] = ubStatus;
      continue;
    }

    if (ubStatus != lastMicroBunchStatus_[linkID]) {
      TGraph* g = g_ubStatusVsEwt_[linkID];
      g->SetPoint(g->GetN(),
                  static_cast<double>(ewt),
                  static_cast<double>(ubStatus));
      lastMicroBunchStatus_[linkID] = ubStatus;
    }
  }
}

void CRVDigiDQM::persistGraph(TGraph* g)
{
  if (!dir_ || g == nullptr) {
    return;
  }
  if (g->GetN() <= 0) {
    dir_->makeAndRegister<TGraph>(g->GetName(), g->GetTitle());
    return;
  }
  dir_->makeAndRegister<TGraph>(g->GetName(),
                                g->GetTitle(),
                                g->GetN(),
                                g->GetX(),
                                g->GetY());
}

int CRVDigiDQM::nDigisOffline(std::size_t channel) const
{
  if (channel >= nDigisOffline_.size()) {
    return 0;
  }
  return nDigisOffline_[channel];
}

void CRVDigiDQM::BookSectorOccupancy(const std::vector<std::string>& sectorNames,
                                     const std::vector<int>& channelToSector,
                                     int nBins, double lo, double hi)
{
  if (!segments_.booked() || !h_sectorOccupancy_.empty()) {
    return;
  }
  channelToSector_ = channelToSector;
  h_sectorOccupancy_.reserve(sectorNames.size());
  for (const auto& sector : sectorNames) {
    const std::string name = "crvDigisPerChannelAndEvent_CRVsector" + sector;
    h_sectorOccupancy_.push_back(segments_.book1<TH1F>(name, name, nBins, lo, hi));
  }
}

void CRVDigiDQM::fillSectorOccupancy()
{
  // Per-file rate distribution. After hadd, rebuild from crvDigisPerChannel / nEvents.
  if (h_sectorOccupancy_.empty() || nEvents_ == 0) {
    return;
  }
  const float invN = 1.0f / static_cast<float>(nEvents_);
  for (std::size_t channel = 0; channel < channelToSector_.size(); ++channel) {
    const int sector = channelToSector_[channel];
    if (sector < 0 ||
        static_cast<std::size_t>(sector) >= h_sectorOccupancy_.size()) {
      continue;
    }
    h_sectorOccupancy_[sector].Fill(nDigisOffline(channel) * invN);
  }
}

void CRVDigiDQM::WriteGraphs()
{
  if (!booked_) {
    return;
  }
  // Overflow bins record that something was off-axis / off-scale, not these
  // extrema, so report them only when they themselves sit past the axis.
  if (maxFebIdSeen_ >= nFebIdBins()) {
    mf::LogWarning("CRVDigiDQM")
        << "max globalFebId " << maxFebIdSeen_ << " is outside the "
        << nFebIdBins() << "-bin FEB axis. Occupancy and dtVsFeb overflow "
        << "bins do not store this value. Reported once at end of job.";
  }
  if (maxAbsDtSeen_ > config_.dtVsFebRange) {
    mf::LogWarning("CRVDigiDQM")
        << "max |dt| " << maxAbsDtSeen_ << " ns exceeds dtVsFebRange "
        << config_.dtVsFebRange << " ns. dtVsFeb overflow and "
        << "dtOutOfRangePerFeb do not store this magnitude. "
        << "Reported once at end of job.";
  }
  if (!dir_) {
    return;
  }
  fillSectorOccupancy();
  // Stamp every copy with the range it actually covers before anything is
  // written; a partially filled live copy must not advertise a full span.
  segments_.Finalize();
  if (!config_.fillLivePlots) {
    return;
  }
  persistGraph(g_digisVsEwt_);
  persistGraph(g_digisAvgVsEwt_);
  for (auto& entry : g_ubStatusVsEwt_) {
    persistGraph(entry.second);
  }
}

} // namespace mu2e
