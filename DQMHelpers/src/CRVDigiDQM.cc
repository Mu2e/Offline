//
// Standalone CRV digi DQM helper.
//
// Original Author: R. Mina
//

#include "Offline/DQMHelpers/inc/CRVDigiDQM.hh"
#include "Offline/DQMHelpers/inc/CRVCFTime.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "TString.h"
#include "TH1.h"

namespace mu2e {

uint8_t CRVDigiDQM::foldedROC(uint8_t roc)
{
  return (roc == kFoldFromROC) ? kFoldToROC : roc;
}

int CRVDigiDQM::globalFebId(uint8_t roc, uint8_t feb)
{
  return (static_cast<int>(foldedROC(roc)) - 1) * kNFebSlotsPerROC + feb;
}

int CRVDigiDQM::globalChannelId(uint8_t roc, uint8_t feb, uint8_t febChannel)
{
  return globalFebId(roc, feb) * kNChanPerFEB + febChannel;
}

CRVDigiDQM::CRVDigiDQM(const Config& config) : config_(config)
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
}

void CRVDigiDQM::Book(art::TFileDirectory dir)
{
  dir_ = dir;
  timingFebDir_ = dir.mkdir("timing_feb");
  timingFpgaDir_ = dir.mkdir("timing_fpga");

  h1_digisPerEvt_ = dir.make<TH1F>("h1_digisPerEvt",
                                   "Hits / event;Hits / event;Events",
                                   config_.nBinsDigisPerEvt,
                                   0.5,
                                   config_.maxDigisPerEvt + 0.5);
  h1_digisPerEvt_->SetMinimum(0.5);

  h1_peakAdc_ = dir.make<TH1F>("h1_peakAdc",
                               "Max sample ADC;Max sample ADC;Hits",
                               config_.nBinsPeakAdc,
                               0,
                               config_.maxPeakAdc);

  h1_tdc_ = dir.make<TH1F>(
      "h1_tdc",
      "Start timestamp of digi in units of 12.5ns;Start timestamp of digi;Digis",
      config_.nBinsTdc,
      0,
      config_.maxTdc);

  h1_channels_ = dir.make<TH1F>("h1_channels",
                                "Channel occupancy;Global channel ID;Hits",
                                kNGlobalChannelBins,
                                -0.5,
                                kNGlobalChannelBins - 0.5);
  h1_channels_->SetMinimum(0.5);

  h1_channelsLastEwt_ = dir.make<TH1F>(
      "h1_channelsLastEwt",
      Form("Channel occupancy (EWT span %zu);Global channel ID;Hits",
           config_.channelsWindowEwts),
      kNGlobalChannelBins,
      -0.5,
      kNGlobalChannelBins - 0.5);
  h1_channelsLastEwt_->SetMinimum(0.5);

  h2_channels_ = dir.make<TH2F>("h2_channels",
                                "FEB vs channel hit map;Channel;FEB",
                                kNChanPerFEB,
                                0.5,
                                kNChanPerFEB + 0.5,
                                kNGlobalFebBins,
                                0.5,
                                kNGlobalFebBins + 0.5);

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

  if (config_.fillInclusive) {
    hBarId_ = dir.make<TH1D>("BarId", "Bar ID", 200, -0.5, 5503.5);
    hSiPM_ = dir.make<TH1D>("SiPM", "SiPM", 4, -0.5, 3.5);
    hADC_ = dir.make<TH1D>("ADC", "ADC in waveform", 100, 0.0, 3000.0);
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

  const bool haveEwt = !crvStatus.empty();
  const uint64_t ewt = haveEwt ? crvStatus.front().GetEventWindowTag() : 0;

  const int nDigis = static_cast<int>(crvDigis.size());
  std::vector<uint16_t> eventChannelHits;
  std::map<uint8_t, std::map<uint8_t, std::vector<FpgaHit>>> hitTimes;

  for (const auto& digi : crvDigis) {
    const uint8_t roc = digi.GetROC();
    const uint8_t feb = digi.GetFEB();
    const uint8_t febChannel = digi.GetFEBchannel();

    const int febId = globalFebId(roc, feb);
    const int channelId = globalChannelId(roc, feb, febChannel);

    h1_channels_->Fill(channelId);
    eventChannelHits.push_back(static_cast<uint16_t>(channelId));
    h2_channels_->Fill(febChannel + 1, febId);

    const auto& adcs = digi.GetADCs();
    if (!adcs.empty()) {
      const int16_t maxSample = *std::max_element(adcs.begin(), adcs.end());
      h1_peakAdc_->Fill(maxSample);
    }

    h1_tdc_->Fill(digi.GetStartTDC());

    if (config_.fillInclusive) {
      if (hBarId_) {
        hBarId_->Fill(digi.GetScintillatorBarIndex().asInt());
      }
      if (hSiPM_) {
        hSiPM_->Fill(digi.GetSiPMNumber());
      }
      if (hADC_) {
        for (auto a : adcs) {
          hADC_->Fill(a);
        }
      }
    }

    const CFResult cf = cfTime(adcs, config_.cfFraction, config_.minAmplitude);
    if (cf.valid) {
      const double absTime_ns =
          cf.time_ns + digi.GetStartTDC() * kDigitizationPeriodNs;
      const uint8_t fpga = febChannel / 16;
      hitTimes[static_cast<uint8_t>(febId)][fpga].push_back(
          {absTime_ns, febChannel});
    }

    const uint8_t rocFolded = foldedROC(roc);
    activeROCs_.insert(rocFolded);
    activeFEBs_.insert(static_cast<uint8_t>(febId));
    rocFEBMap_[rocFolded].insert(feb);
  }

  nDigis_ += static_cast<std::size_t>(nDigis);
  h1_digisPerEvt_->Fill(nDigis);

  fillTiming(hitTimes);

  if (haveEwt) {
    fillEwtSeries(ewt, nDigis);
    fillRollingOccupancy(ewt, eventChannelHits);
    fillMicroBunchStatus(crvStatus);
  }
}

void CRVDigiDQM::fillEwtSeries(uint64_t ewt, int nDigis)
{
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

void CRVDigiDQM::fillRollingOccupancy(
    uint64_t ewt, const std::vector<uint16_t>& eventChannelHits)
{
  recentChannelHitsByEwt_.emplace_back(ewt, eventChannelHits);
  for (const auto channelId : recentChannelHitsByEwt_.back().second) {
    if (channelId < kNGlobalChannelBins) {
      h1_channelsLastEwt_->AddBinContent(static_cast<Int_t>(channelId) + 1, 1.0);
    }
  }

  const uint64_t minKeepEwt =
      (ewt > config_.channelsWindowEwts) ? ewt - config_.channelsWindowEwts : 0;
  while (!recentChannelHitsByEwt_.empty() &&
         recentChannelHitsByEwt_.front().first < minKeepEwt) {
    for (const auto channelId : recentChannelHitsByEwt_.front().second) {
      if (channelId < kNGlobalChannelBins) {
        h1_channelsLastEwt_->AddBinContent(static_cast<Int_t>(channelId) + 1,
                                           -1.0);
      }
    }
    recentChannelHitsByEwt_.pop_front();
  }
}

void CRVDigiDQM::fillTiming(
    const std::map<uint8_t, std::map<uint8_t, std::vector<FpgaHit>>>& hitTimes)
{
  if (!timingFebDir_ || !timingFpgaDir_) {
    return;
  }

  std::vector<std::pair<uint8_t, double>> febFirstHit;
  for (const auto& [febId, fpgaMap] : hitTimes) {
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

  for (std::size_t i = 0; i < febFirstHit.size(); ++i) {
    for (std::size_t j = i + 1; j < febFirstHit.size(); ++j) {
      const uint8_t lo = febFirstHit[i].first;
      const uint8_t hi = febFirstHit[j].first;
      const double dt = febFirstHit[j].second - febFirstHit[i].second;
      const auto key = std::make_pair(lo, hi);

      if (h1_dtFebPairs_.find(key) == h1_dtFebPairs_.end()) {
        const std::string name = Form("dt_feb%02d_feb%02d", lo, hi);
        const std::string title =
            Form("#Deltat FEB %02d - FEB %02d;#Deltat [ns];Entries", lo, hi);
        h1_dtFebPairs_[key] = timingFebDir_->make<TH1F>(name.c_str(),
                                                        title.c_str(),
                                                        nBinsDt_,
                                                        -config_.dtRange,
                                                        config_.dtRange);
      }
      h1_dtFebPairs_[key]->Fill(dt);
    }
  }

  for (const auto& [febId, fpgaMap] : hitTimes) {
    for (auto itA = fpgaMap.begin(); itA != fpgaMap.end(); ++itA) {
      for (auto itB = itA; itB != fpgaMap.end(); ++itB) {
        const uint8_t fpgaA = itA->first;
        const uint8_t fpgaB = itB->first;
        const auto& hitsA = itA->second;
        const auto& hitsB = itB->second;

        const uint8_t pairCode = static_cast<uint8_t>(fpgaA * 4 + fpgaB);
        const auto key = std::make_pair(febId, pairCode);

        if (h1_dtFpgaPairs_.find(key) == h1_dtFpgaPairs_.end()) {
          const std::string name =
              Form("dt_feb%02d_fpga%d_fpga%d", febId, fpgaA, fpgaB);
          const std::string title = Form(
              "#Deltat FEB %02d FPGA %d - FPGA %d;#Deltat [ns];Entries",
              febId,
              fpgaA,
              fpgaB);
          h1_dtFpgaPairs_[key] = timingFpgaDir_->make<TH1F>(name.c_str(),
                                                            title.c_str(),
                                                            nBinsDt_,
                                                            -config_.dtRange,
                                                            config_.dtRange);
        }

        TH1F* h = h1_dtFpgaPairs_[key];
        if (fpgaA == fpgaB) {
          for (std::size_t ia = 0; ia < hitsA.size(); ++ia) {
            for (std::size_t ib = ia + 1; ib < hitsA.size(); ++ib) {
              if (hitsA[ia].channel != hitsA[ib].channel) {
                h->Fill(hitsA[ib].time_ns - hitsA[ia].time_ns);
              }
            }
          }
        } else {
          for (const auto& hA : hitsA) {
            for (const auto& hB : hitsB) {
              h->Fill(hB.time_ns - hA.time_ns);
            }
          }
        }
      }
    }
  }
}

void CRVDigiDQM::fillMicroBunchStatus(const CrvStatusCollection& crvStatus)
{
  if (!dir_) {
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

void CRVDigiDQM::WriteGraphs()
{
  if (!booked_ || !dir_) {
    return;
  }
  persistGraph(g_digisVsEwt_);
  persistGraph(g_digisAvgVsEwt_);
  for (auto& entry : g_ubStatusVsEwt_) {
    persistGraph(entry.second);
  }
}

} // namespace mu2e
