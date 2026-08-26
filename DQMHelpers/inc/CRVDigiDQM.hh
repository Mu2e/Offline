#ifndef DQMHelpers_inc_CRVDigiDQM_hh
#define DQMHelpers_inc_CRVDigiDQM_hh
//
// Standalone CRV digi DQM helper. Books and fills the histograms used by both
// the otsdaq online monitor and offline DQM art modules.
//
// Channel-ID convention (VST/KPP slot map, not CRVId.hh):
//   if (roc == 4) roc = 2;                 // DTC link 3 folded onto ROC 2
//   globalFebId     = (roc-1)*25 + feb;    // 25 FEB slots per ROC
//   globalChannelId = globalFebId*64 + febChannel;  // 2112 occupancy bins
//
// Original Author: R. Mina
//

#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art_root_io/TFileDirectory.h"

#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

#include <cstdint>
#include <deque>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace mu2e {

class CRVDigiDQM {
public:
  // Binning and CF-timing defaults match otsdaq-mu2e-crv CrvDQM_module.cc
  // (branch mu2e/ots_ops).
  struct Config {
    int nBinsDigisPerEvt{200};
    float maxDigisPerEvt{4000.f};
    int nBinsPeakAdc{450};
    float maxPeakAdc{4500.f};
    int nBinsTdc{400};
    float maxTdc{40000.f};
    double cfFraction{0.20};
    float dtBinSize{0.5f};
    float dtRange{100.f};
    int minAmplitude{10};
    std::size_t avgBlockSize{30};
    std::size_t avgGraphPoints{1000};
    std::size_t channelsWindowEwts{50000};
    bool fillInclusive{true};
  };

  // VST/KPP readout geography used by the online occupancy plots.
  static constexpr int kNFebSlotsPerROC = 25;
  static constexpr int kNChanPerFEB = 64;
  static constexpr int kNGlobalChannelBins = 2112;
  static constexpr int kNGlobalFebBins = 32;
  static constexpr uint8_t kFoldFromROC = 4;
  static constexpr uint8_t kFoldToROC = 2;

  static constexpr std::size_t kEwtWindow = 1000;
  static constexpr std::size_t kGraphPoints = 10000;
  static constexpr double kEwtXRange = 1000000;

  static uint8_t foldedROC(uint8_t roc);
  static int globalFebId(uint8_t roc, uint8_t feb);
  static int globalChannelId(uint8_t roc, uint8_t feb, uint8_t febChannel);

  explicit CRVDigiDQM(const Config& config);

  void Book(art::TFileDirectory dir);
  void Fill(const CrvDigiCollection& crvDigis,
            const CrvStatusCollection& crvStatus);
  void WriteGraphs();

  TH1F* h1_digisPerEvt() const { return h1_digisPerEvt_; }
  TH1F* h1_peakAdc() const { return h1_peakAdc_; }
  TH1F* h1_tdc() const { return h1_tdc_; }
  TH1F* h1_channels() const { return h1_channels_; }
  TH1F* h1_channelsLastEwt() const { return h1_channelsLastEwt_; }
  TH2F* h2_channels() const { return h2_channels_; }
  TGraph* g_digisVsEwt() const { return g_digisVsEwt_; }
  TGraph* g_digisAvgVsEwt() const { return g_digisAvgVsEwt_; }

  TH1D* BarId() const { return hBarId_; }
  TH1D* SiPM() const { return hSiPM_; }
  TH1D* ADC() const { return hADC_; }

  const std::map<std::pair<uint8_t, uint8_t>, TH1F*>& dtFebPairs() const
  {
    return h1_dtFebPairs_;
  }
  const std::map<std::pair<uint8_t, uint8_t>, TH1F*>& dtFpgaPairs() const
  {
    return h1_dtFpgaPairs_;
  }
  const std::map<uint8_t, TGraph*>& ubStatusVsEwt() const
  {
    return g_ubStatusVsEwt_;
  }

  std::size_t nEvents() const { return nEvents_; }
  std::size_t nDigis() const { return nDigis_; }
  bool hasEwtWindow() const { return !ewtWindow_.empty(); }
  uint64_t lastEwt() const
  {
    return ewtWindow_.empty() ? 0 : ewtWindow_.back().first;
  }
  const std::set<uint8_t>& activeFEBs() const { return activeFEBs_; }
  const std::set<uint8_t>& activeROCs() const { return activeROCs_; }
  const std::map<uint8_t, std::set<uint8_t>>& rocFEBMap() const
  {
    return rocFEBMap_;
  }

private:
  struct FpgaHit {
    double time_ns;
    uint8_t channel;
  };

  void fillEwtSeries(uint64_t ewt, int nDigis);
  void fillRollingOccupancy(uint64_t ewt,
                            const std::vector<uint16_t>& eventChannelHits);
  void fillTiming(const std::map<uint8_t, std::map<uint8_t, std::vector<FpgaHit>>>& hitTimes);
  void fillMicroBunchStatus(const CrvStatusCollection& crvStatus);
  void persistGraph(TGraph* g);

  Config config_;
  int nBinsDt_{400};
  bool booked_{false};

  std::optional<art::TFileDirectory> dir_;
  std::optional<art::TFileDirectory> timingFebDir_;
  std::optional<art::TFileDirectory> timingFpgaDir_;

  TH1F* h1_digisPerEvt_{nullptr};
  TH1F* h1_peakAdc_{nullptr};
  TH1F* h1_tdc_{nullptr};
  TH1F* h1_channels_{nullptr};
  TH1F* h1_channelsLastEwt_{nullptr};
  TH2F* h2_channels_{nullptr};
  TGraph* g_digisVsEwt_{nullptr};
  TGraph* g_digisAvgVsEwt_{nullptr};

  TH1D* hBarId_{nullptr};
  TH1D* hSiPM_{nullptr};
  TH1D* hADC_{nullptr};

  std::map<std::pair<uint8_t, uint8_t>, TH1F*> h1_dtFebPairs_;
  std::map<std::pair<uint8_t, uint8_t>, TH1F*> h1_dtFpgaPairs_;
  std::map<uint8_t, TGraph*> g_ubStatusVsEwt_;
  std::map<uint8_t, uint32_t> lastMicroBunchStatus_;

  std::size_t nEvents_{0};
  std::size_t nDigis_{0};
  std::set<uint8_t> activeFEBs_;
  std::set<uint8_t> activeROCs_;
  std::map<uint8_t, std::set<uint8_t>> rocFEBMap_;

  std::deque<std::pair<uint64_t, int>> ewtWindow_;
  long long ewtWindowSum_{0};
  std::deque<std::pair<uint64_t, std::vector<uint16_t>>> recentChannelHitsByEwt_;

  long long avgBlockSum_{0};
  std::size_t avgBlockCount_{0};
  uint64_t avgBlockFirstEwt_{0};
  bool avgSeedsCleared_{false};
};

} // namespace mu2e

#endif /* DQMHelpers_inc_CRVDigiDQM_hh */
