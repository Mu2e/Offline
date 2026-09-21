// CRV digi DQM client.
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVDigiDQM.hh"
#include "Offline/CRVDQM/inc/CRVCFTime.hh"

#include "TString.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>

namespace mu2e {

using namespace CRVDQMRun1;

const char* CRVDigiDQM::dtClassName(int c)
{
  switch (c) {
    case kDtSameModuleSameSide: return "SameModuleSameSide";
    case kDtSameModuleOtherSide: return "SameModuleOtherSide";
    case kDtAdjModuleSameSide: return "AdjModuleSameSide";
    case kDtAdjModuleOtherSide: return "AdjModuleOtherSide";
    default: return "Unknown";
  }
}

CRVDigiDQM::CRVDigiDQM(const DQMHistSet::Config& hists) :
    DQMClient("CRVDigiDQM", kBinningVersion, hists),
    nDigisOffline_(kNOfflineChannels, 0)
{}

TH2F* CRVDigiDQM::dtPartner(int dtClass) const
{
  if (dtClass < 0 || static_cast<std::size_t>(dtClass) >= h2_dtPartner_.size()) {
    return nullptr;
  }
  return h2_dtPartner_[dtClass];
}

void CRVDigiDQM::book()
{
  DQMHistSet& h = hists();

  h1_digisPerEvt_ = h.book1<TH1F>("h1_digisPerEvt", "Hits / event;Hits / event;Events",
                                  kDigisPerEvt);
  h1_digisPerEvt2_ = h.book1<TH1F>("h1_digisPerEvt2", "Hits / event;Hits / event;Events",
                                   kDigisPerEvt2);
  h1_peakAdc_ = h.book1<TH1F>("h1_peakAdc", "Max sample ADC;Max sample ADC;Hits", kPeakAdc);
  h1_tdc_ = h.book1<TH1F>(
      "h1_tdc", "Start timestamp of digi in units of 12.5ns;Start timestamp of digi;Digis",
      kTdc);
  h1_tdc2_ = h.book1<TH1F>(
      "h1_tdc2", "Start timestamp of digi in units of 12.5ns;Start timestamp of digi;Digis",
      kTdc2);

  h1_channels_ = h.book1<TH1F>(
      "h1_channels", "Channel occupancy;Online channel (FEB port#times64+FEB channel);Hits",
      kOnlineChannel);
  h2_channels_ = h.book2<TH2F>(
      "h2_channels", "FEB vs channel hit map;FEB channel;FEB port (ROC-1)#times24+(FEB-1)",
      kFebChannel, kFebPort);

  hBarId_ = h.book1<TH1D>("BarId", "Bar ID", kBarId);
  hSiPM_ = h.book1<TH1D>("SiPM", "SiPM", kSiPM);
  hADC_ = h.book1<TH1D>("ADC", "ADC in waveform", kAdcSample);

  h_crvDigiRatesROC_.resize(kNROC);
  for (int roc = 1; roc <= kNROC; ++roc) {
    h_crvDigiRatesROC_[roc - 1] =
        h.book1<TH1F>(Form("crvDigiRates_ROC%d", roc),
                      Form("crvDigiRates_ROC%d;Online channel in ROC;Digis", roc),
                      kRocChannelEdges);
  }
  h_crvDigiRates_ = h.book2<TH2F>("crvDigiRates",
                                  "crvDigiRates:FEBchannel:FEB;FEB channel;FEB port",
                                  kFebChannelEdges, kFebPortEdges);
  h_crvDigisPerChannel_ = h.book1<TH1F>(
      "crvDigisPerChannel", "Digis vs offline channel;Offline channel (bar#times4+SiPM);Digis",
      kOfflineChannel);

  h2_dtFpgaPairs_ = h.book2<TH2F>(
      "dtFpgaPairs",
      "#Deltat between hits on one FEB;FEB port#times10 + FPGA pair;#Deltat [ns]",
      kFpgaPair, kDtFpga);

  h2_dtPartner_.resize(kNDtClasses);
  for (int c = 0; c < kNDtClasses; ++c) {
    h2_dtPartner_[c] = h.book2<TH2F>(
        Form("dtPartner_%s", dtClassName(c)),
        Form("#Deltat to partner FEB, %s;FEB port;#Deltat [ns]", dtClassName(c)),
        kFebPort, kDtPartner[c]);
  }
  h_layersPerGroup_ = h.book1<TH1F>(
      "layersPerGroup", "Layers in a coincidence group;Layers;Groups", kLayersPerGroup);
  h_groupsPerEvent_ = h.book1<TH1F>(
      "groupsPerEvent", "Coincidence groups per event;Groups;Events", kGroupsPerEvent);
  h_groupsPerEvent2_ = h.book1<TH1F>(
      "groupsPerEvent2", "Coincidence groups per event (air showers);Groups;Events",
      kGroupsPerEvent2);
  h_sectorsPerEvent_ = h.book1<TH1F>(
      "sectorsPerEvent", "CRV sectors with a group;Sectors;Events", kSectorsPerEvent);
  h_febNoGroup_ = h.book1<TH1F>(
      "febNoGroup", "Events where this FEB had hits but no group formed;FEB port;Events",
      kFebPort);

  // Every configuration's family is booked; a job fills only its own.
  h_sectorOccupancy_.resize(kNConfigurations);
  h_sectorOccupancy2_.resize(kNConfigurations);
  for (int c = 0; c < kNConfigurations; ++c) {
    for (int s = 0; s < kLayouts[c].nSectors; ++s) {
      const std::string tag = "_CRVsector" + sectorTag(c, s);
      const std::string title = ";Digis per channel and event;Channels";
      h_sectorOccupancy_[c].push_back(h.bookSummary1<TH1F>(
          "crvDigisPerChannelAndEvent" + tag, "crvDigisPerChannelAndEvent" + tag + title,
          kDigisPerChannelAndEvent));
      h_sectorOccupancy2_[c].push_back(h.bookSummary1<TH1F>(
          "crvDigisPerChannelAndEvent2" + tag, "crvDigisPerChannelAndEvent2" + tag + title,
          kDigisPerChannelAndEvent2));
    }
  }

  g_digisVsEwt_ = &series().book(
      "g_digisVsEwt",
      Form("Hits in last %zu events;Event window tag;Hits", kEwtWindow));
  g_digisAvgVsEwt_ = &series().book(
      "g_digisAvgVsEwt",
      Form("Mean hits per event (averaged over %zu events);Event window tag;Hits",
           kAvgBlockSize),
      kAvgGraphPoints);
}

void CRVDigiDQM::SetConfiguration(int configuration, const std::vector<int>& channelToSector)
{
  configuration_ = configuration;
  channelToSector_ = channelToSector;
}

void CRVDigiDQM::SetFebTopology(const std::vector<FebTopology>& febTopology,
                                const std::vector<int>& channelToLayer)
{
  febTopology_ = febTopology;
  channelToLayer_ = channelToLayer;
}

void CRVDigiDQM::Fill(const CrvDigiCollection& crvDigis, const CrvStatusCollection& crvStatus)
{
  if (!booked()) {
    return;
  }
  const bool haveEwt = !crvStatus.empty();
  const uint64_t ewt = haveEwt ? crvStatus.front().GetEventWindowTag() : 0;
  beginEvent(haveEwt ? std::optional<uint64_t>(ewt) : std::nullopt);

  const int nDigis = static_cast<int>(crvDigis.size());
  std::map<int, std::map<uint8_t, std::vector<FpgaHit>>> hitTimes;  //by FEB port
  std::vector<PartnerHit> partnerHits;

  for (const auto& digi : crvDigis) {
    const int roc = static_cast<int>(digi.GetROC());
    const int feb = static_cast<int>(digi.GetFEB());
    const int febChannel = static_cast<int>(digi.GetFEBchannel());
    const bool online = onlineIdInRange(roc, feb, febChannel);
    const int port = online ? febPort(roc, feb) : -1;

    if (online) {
      h1_channels_.Fill(onlineChannel(roc, feb, febChannel));
      h2_channels_.Fill(febChannel, port);
      h_crvDigiRatesROC_[roc - 1].Fill(rocChannel(feb, febChannel));
      h_crvDigiRates_.Fill(febChannel, port);
      activeFebPorts_.insert(port);
    } else {
      diag().Count("onlineIdOutOfRange",
                   "a digi's ROC/FEB/channel is outside CRVId (ROC 1-18, FEB 1-24, "
                   "channel 0-63); it is left out of the online-indexed histograms "
                   "and the timing.");
    }

    const auto& adcs = digi.GetADCs();
    if (!adcs.empty()) {
      h1_peakAdc_.Fill(*std::max_element(adcs.begin(), adcs.end()));
    }
    h1_tdc_.Fill(digi.GetStartTDC());
    h1_tdc2_.Fill(digi.GetStartTDC());

    const int barIndex = digi.GetScintillatorBarIndex().asInt();
    const int sipm = digi.GetSiPMNumber();
    hBarId_.Fill(barIndex);
    hSiPM_.Fill(sipm);
    for (auto a : adcs) {
      hADC_.Fill(a);
    }
    int offline = -1;
    if (sipm >= 0 && barIndex >= 0) {
      offline = barIndex * static_cast<int>(CRVId::nChanPerBar) + sipm;
      if (offline < kNOfflineChannels) {
        ++nDigisOffline_[offline];
        h_crvDigisPerChannel_.Fill(offline);
      } else {
        offline = -1;
      }
    }

    const CFResult cf = cfTime(adcs, kCFFraction, kCFMinAmplitude);
    if (cf.valid && online) {
      const double absTime_ns = cf.time_ns + digi.GetStartTDC() * CRVDigitizationPeriod;
      const uint8_t fpga = static_cast<uint8_t>(febChannel / kNChanPerFPGA);
      hitTimes[port][fpga].push_back({absTime_ns, static_cast<uint8_t>(febChannel)});

      // Partner timing takes a harder amplitude cut: below it the sample is
      // dominated by hits uncorrelated with the traversal.
      const int amplitude = static_cast<int>(cf.peak) - cf.baseline;
      if (amplitude >= kDtMinAmplitude && offline >= 0 &&
          static_cast<std::size_t>(port) < febTopology_.size() &&
          febTopology_[port].valid &&
          static_cast<std::size_t>(offline) < channelToLayer_.size() &&
          channelToLayer_[offline] >= 0) {
        const FebTopology& t = febTopology_[port];
        partnerHits.push_back(
            {absTime_ns, port, t.sector, t.module, t.side, channelToLayer_[offline]});
      }
    }

    activeROCs_.insert(static_cast<uint8_t>(roc));
    rocFEBMap_[static_cast<uint8_t>(roc)].insert(static_cast<uint8_t>(feb));
  }

  nDigis_ += static_cast<std::size_t>(nDigis);
  h1_digisPerEvt_.Fill(nDigis);
  h1_digisPerEvt2_.Fill(nDigis);

  fillFpgaTiming(hitTimes);
  fillPartnerTiming(partnerHits);

  if (haveEwt) {
    fillEwtSeries(ewt, nDigis);
  }
}

// Every hit pair on one FEB: different channels on one FPGA, or any two FPGAs.
// A clock or FPGA offset moves the correlated core of its column.
void CRVDigiDQM::fillFpgaTiming(
    const std::map<int, std::map<uint8_t, std::vector<FpgaHit>>>& hitTimes)
{
  for (const auto& [port, fpgaMap] : hitTimes) {
    for (auto itA = fpgaMap.begin(); itA != fpgaMap.end(); ++itA) {
      for (auto itB = itA; itB != fpgaMap.end(); ++itB) {
        const double x = port * kNFpgaPairs + fpgaPairIndex(itA->first, itB->first);
        const auto& hitsA = itA->second;
        const auto& hitsB = itB->second;
        if (itA == itB) {
          for (std::size_t ia = 0; ia < hitsA.size(); ++ia) {
            for (std::size_t ib = ia + 1; ib < hitsA.size(); ++ib) {
              if (hitsA[ia].channel != hitsA[ib].channel) {
                h2_dtFpgaPairs_.Fill(x, hitsA[ib].time_ns - hitsA[ia].time_ns);
              }
            }
          }
        } else {
          for (const auto& hA : hitsA) {
            for (const auto& hB : hitsB) {
              h2_dtFpgaPairs_.Fill(x, hB.time_ns - hA.time_ns);
            }
          }
        }
      }
    }
  }
}

int CRVDigiDQM::dtClassFor(int febA, int febB) const
{
  const FebTopology& a = febTopology_[febA];
  const FebTopology& b = febTopology_[febB];
  if (a.sector != b.sector) {
    return -1;
  }
  const int dModule = std::abs(a.module - b.module);
  const bool sameSide = a.side == b.side;
  if (dModule == 0) {
    //same module, same side means the other layer pair: one FEB covers two layers
    return sameSide ? kDtSameModuleSameSide : kDtSameModuleOtherSide;
  }
  if (dModule == 1) {
    return sameSide ? kDtAdjModuleSameSide : kDtAdjModuleOtherSide;
  }
  return -1;  //not a partner: no timing relation is expected
}

// One qualifying group: every ordered FEB pair that is a partner contributes,
// so a slipped FEB appears as a displaced column in its class.
void CRVDigiDQM::fillGroup(const std::vector<const PartnerHit*>& group)
{
  std::map<int, double> earliest;  //per FEB port
  std::set<int> layers;
  for (const PartnerHit* hit : group) {
    layers.insert(hit->layer);
    auto it = earliest.find(hit->febPort);
    if (it == earliest.end() || hit->time_ns < it->second) {
      earliest[hit->febPort] = hit->time_ns;
    }
  }

  h_layersPerGroup_.Fill(static_cast<double>(layers.size()));
  ++nGroups_;

  for (const auto& [febA, tA] : earliest) {
    for (const auto& [febB, tB] : earliest) {
      if (febA == febB) {
        continue;
      }
      const int c = dtClassFor(febA, febB);
      if (c >= 0) {
        h2_dtPartner_[c].Fill(febA, tA - tB);
      }
    }
  }
}

void CRVDigiDQM::fillPartnerTiming(std::vector<PartnerHit>& hits)
{
  if (hits.empty()) {
    h_groupsPerEvent_.Fill(0.);
    h_groupsPerEvent2_.Fill(0.);
    h_sectorsPerEvent_.Fill(0.);
    return;
  }

  // Sector first, then time: a muon in one sector never times against another,
  // and the two sides of the tracker are separate traversals.
  std::sort(hits.begin(), hits.end(), [](const PartnerHit& a, const PartnerHit& b) {
    return a.sector != b.sector ? a.sector < b.sector : a.time_ns < b.time_ns;
  });

  std::set<int> febsInGroup;
  std::set<int> sectorsWithGroup;
  std::set<int> febsWithHits;
  std::size_t nGroupsThisEvent = 0;
  for (const PartnerHit& hit : hits) {
    febsWithHits.insert(hit.febPort);
  }

  std::size_t i = 0;
  while (i < hits.size()) {
    //one time cluster: consecutive hits of one sector no further apart than the window
    std::size_t j = i + 1;
    while (j < hits.size() && hits[j].sector == hits[i].sector &&
           hits[j].time_ns - hits[j - 1].time_ns <= kDtCoincWindow) {
      ++j;
    }

    // Within the cluster, a group is one module, or two adjacent modules when
    // the muon crossed between them. Layers are counted over the group.
    std::map<int, std::vector<const PartnerHit*>> byModule;
    for (std::size_t k = i; k < j; ++k) {
      byModule[hits[k].module].push_back(&hits[k]);
    }

    std::set<int> used;
    for (const auto& [module, moduleHits] : byModule) {
      if (used.count(module) > 0) {
        continue;
      }
      std::set<int> layers;
      for (const PartnerHit* hit : moduleHits) {
        layers.insert(hit->layer);
      }
      std::vector<const PartnerHit*> group = moduleHits;
      std::set<int> merged{module};
      if (static_cast<int>(layers.size()) < kDtMinLayers) {
        //a muon crossing between adjacent modules leaves its layers split
        auto next = byModule.find(module + 1);
        if (next != byModule.end() && used.count(module + 1) == 0) {
          for (const PartnerHit* hit : next->second) {
            layers.insert(hit->layer);
            group.push_back(hit);
          }
          merged.insert(module + 1);
        }
      }
      if (static_cast<int>(layers.size()) < kDtMinLayers) {
        continue;
      }
      used.insert(merged.begin(), merged.end());
      fillGroup(group);
      ++nGroupsThisEvent;
      sectorsWithGroup.insert(hits[i].sector);
      for (const PartnerHit* hit : group) {
        febsInGroup.insert(hit->febPort);
      }
    }
    i = j;
  }

  h_groupsPerEvent_.Fill(static_cast<double>(nGroupsThisEvent));
  h_groupsPerEvent2_.Fill(static_cast<double>(nGroupsThisEvent));
  h_sectorsPerEvent_.Fill(static_cast<double>(sectorsWithGroup.size()));
  for (int feb : febsWithHits) {
    if (febsInGroup.count(feb) == 0) {
      h_febNoGroup_.Fill(feb);
    }
  }
}

void CRVDigiDQM::fillEwtSeries(uint64_t ewt, int nDigis)
{
  if (!series().enabled()) {
    return;
  }
  if (avgBlockCount_ == 0) {
    avgBlockFirstEwt_ = ewt;
  }
  avgBlockSum_ += nDigis;
  ++avgBlockCount_;
  if (avgBlockCount_ >= kAvgBlockSize) {
    const double midEwt =
        0.5 * (static_cast<double>(avgBlockFirstEwt_) + static_cast<double>(ewt));
    g_digisAvgVsEwt_->Add(midEwt, static_cast<double>(avgBlockSum_) /
                                      static_cast<double>(avgBlockCount_));
    avgBlockSum_ = 0;
    avgBlockCount_ = 0;
  }

  ewtWindow_.emplace_back(ewt, nDigis);
  ewtWindowSum_ += nDigis;
  while (ewtWindow_.size() > kEwtWindow) {
    ewtWindowSum_ -= ewtWindow_.front().second;
    ewtWindow_.pop_front();
  }
  g_digisVsEwt_->Add(static_cast<double>(ewt), static_cast<double>(ewtWindowSum_));
}

// Per-file rate distribution; after hadd, rebuild it from crvDigisPerChannel / nEvents.
void CRVDigiDQM::fillSectorOccupancy()
{
  if (configuration_ < 0 || nEvents() == 0) {
    if (configuration_ < 0) {
      diag().Count("noConfiguration",
                   "no CRV configuration was set; the per-sector occupancy stays empty.");
    }
    return;
  }
  const double invN = 1.0 / static_cast<double>(nEvents());
  const auto& family = h_sectorOccupancy_[configuration_];
  const auto& family2 = h_sectorOccupancy2_[configuration_];
  for (std::size_t channel = 0; channel < channelToSector_.size(); ++channel) {
    const int sector = channelToSector_[channel];
    if (sector < 0 || static_cast<std::size_t>(sector) >= family.size() ||
        channel >= nDigisOffline_.size()) {
      continue;
    }
    family[sector].Fill(nDigisOffline_[channel] * invN);
    family2[sector].Fill(nDigisOffline_[channel] * invN);
  }
}

void CRVDigiDQM::endJob()
{
  fillSectorOccupancy();
}

void CRVDigiDQM::resetForNewRun()
{
  ewtWindow_.clear();
  ewtWindowSum_ = 0;
  avgBlockSum_ = 0;
  avgBlockCount_ = 0;
  avgBlockFirstEwt_ = 0;
  nGroups_ = 0;
}

} // namespace mu2e
