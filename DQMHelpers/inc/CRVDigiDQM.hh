#ifndef DQMHelpers_inc_CRVDigiDQM_hh
#define DQMHelpers_inc_CRVDigiDQM_hh
//
// Standalone CRV digi DQM helper. Books and fills the histograms used by both
// the otsdaq online monitor and offline DQM art modules. No GeometryService.
//
// Two channel-ID conventions, selected by Config::kppReadout (never both):
//   kppReadout=true (default, and the only mode any existing data needs) -
//     KPP occupancy (h1_channels / h2_channels):
//     if (roc == 4) roc = 2;                 // DTC link 3 folded onto ROC 2
//     globalFebId     = (roc-1)*25 + feb;    // 25 FEB slots per ROC
//     globalChannelId = globalFebId*64 + febChannel;  // 2112 occupancy bins
//   kppReadout=false - full CRV, which does not exist yet. No fold, and the
//     occupancy pair is not booked: with the fold on, ROC 2 and ROC 4 would
//     merge, and 394 of 432 FEBs fall past the 2112-bin axis. crvDigisPerChannel
//     and crvDigiRates below already cover the full detector correctly binned.
//   CRVId rate maps (crvDigiRates_ROC*, crvDigiRates, crvDigisPerChannel):
//     raw GetROC()/GetFEB()/GetFEBchannel() (no fold; CRVId 24 FEBs/ROC).
//     Offline channel = barIndex*4 + SiPM. Scaled by 1/nEvents in WriteGraphs.
//     Per-sector crvDigisPerChannelAndEvent_* is filled by the art module.
//
// Inter-FEB sync (dtVsFeb): each FEB's first CF hit time is compared against the
// median of the other FEBs' first hit times in the same event, so a FEB whose
// clock has slipped shows up as a displaced vertical stripe at its globalFebId.
// This replaces the former per-FEB-pair dt_febXX_febYY histograms, which encoded
// the same N offsets in O(N^2) histograms. Intra-FEB FPGA timing is unchanged.
//
// Original Author: R. Mina
//

#include "Offline/DataProducts/inc/CRVId.hh"
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
    // dtVsFeb needs a range wide enough to show a firmware clock slip, which is
    // much larger than the few-ns spread dtRange is binned for.
    float dtVsFebBinSize{2.f};
    float dtVsFebRange{500.f};
    int minAmplitude{10};
    std::size_t avgBlockSize{30};
    std::size_t avgGraphPoints{1000};
    std::size_t channelsWindowEwts{50000};
    bool fillInclusive{true};
    // CRVId-indexed rate maps and offline-channel occupancy (no GeometryService).
    bool fillCrvIdRates{true};
    // Single-DTC KPP cabling: fold ROC 4 onto ROC 2 and book the 33-slot
    // occupancy pair. False for full CRV; see the header comment above.
    bool kppReadout{true};
  };

  // KPP readout geography used by the online occupancy plots.
  static constexpr int kNFebSlotsPerROC = 25;
  static constexpr int kNChanPerFEB = 64;
  static constexpr int kNGlobalChannelBins = 2112;
  static constexpr int kNGlobalFebBins = 32;
  static constexpr uint8_t kFoldFromROC = 4;
  static constexpr uint8_t kFoldToROC = 2;

  // Full-CRV globalFebId range, used only when kppReadout is false.
  static constexpr int kNFebIdBins = kNFebSlotsPerROC * static_cast<int>(CRVId::nROC);

  // dtVsFeb x-axis: same geography as the occupancy axes, so both modes stay
  // readable and the two stay consistent by construction.
  int nFebIdBins() const
  {
    return config_.kppReadout ? kNGlobalFebBins + 1 : kNFebIdBins;
  }

  static constexpr std::size_t kEwtWindow = 1000;
  static constexpr std::size_t kGraphPoints = 10000;
  static constexpr double kEwtXRange = 1000000;

  // Not static: the fold and the FEB stride depend on Config::kppReadout.
  uint8_t foldedROC(uint8_t roc) const;
  int globalFebId(uint8_t roc, uint8_t feb) const;
  int globalChannelId(uint8_t roc, uint8_t feb, uint8_t febChannel) const;

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

  TH1F* crvDigisPerChannel() const { return h_crvDigisPerChannel_; }
  TH2F* crvDigiRates() const { return h_crvDigiRates_; }
  const std::vector<TH1F*>& crvDigiRatesROC() const { return h_crvDigiRatesROC_; }
  int nDigisOffline(std::size_t channel) const;
  std::size_t nOfflineChannels() const { return nDigisOffline_.size(); }
  bool ratesScaled() const { return ratesScaled_; }

  TH2F* dtVsFeb() const { return h2_dtVsFeb_; }

  // Per-FEB desync counters: a FEB whose clock has slipped stands above the
  // others. Cumulative for offline, rolling EWT window for the online monitor.
  TH1F* dtOutOfRangePerFeb() const { return h_dtOutOfRangePerFeb_; }
  TH1F* dtOutOfRangePerFebLastEwt() const
  {
    return h_dtOutOfRangePerFebLastEwt_;
  }
  const std::map<std::pair<int, uint8_t>, TH1F*>& dtFpgaPairs() const
  {
    return h1_dtFpgaPairs_;
  }
  const std::map<uint8_t, TGraph*>& ubStatusVsEwt() const
  {
    return g_ubStatusVsEwt_;
  }

  // Axis-coverage diagnostics. Only the off-axis FEB case is logged (once per
  // job, from Fill); the dt counts are reported through the histograms above.
  int maxFebIdSeen() const { return maxFebIdSeen_; }
  long long nFebIdOutOfAxis() const { return nFebIdOutOfAxis_; }
  long long nDtOutOfRange() const { return nDtOutOfRange_; }
  double maxAbsDtSeen() const { return maxAbsDtSeen_; }

  std::size_t nEvents() const { return nEvents_; }
  std::size_t nDigis() const { return nDigis_; }
  bool hasEwtWindow() const { return !ewtWindow_.empty(); }
  uint64_t lastEwt() const
  {
    return ewtWindow_.empty() ? 0 : ewtWindow_.back().first;
  }
  const std::set<int>& activeFEBs() const { return activeFEBs_; }
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
  void fillTiming(const std::map<int, std::map<uint8_t, std::vector<FpgaHit>>>& hitTimes);
  void fillMicroBunchStatus(const CrvStatusCollection& crvStatus);
  void persistGraph(TGraph* g);
  void scaleRateHists();
  void fillRollingDtOutOfRange(uint64_t ewt);

  Config config_;
  int nBinsDt_{400};
  int nBinsDtVsFeb_{500};
  bool booked_{false};

  std::optional<art::TFileDirectory> dir_;
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

  TH1F* h_crvDigisPerChannel_{nullptr};
  TH2F* h_crvDigiRates_{nullptr};
  std::vector<TH1F*> h_crvDigiRatesROC_;
  std::vector<int> nDigisOffline_;
  bool ratesScaled_{false};

  TH2F* h2_dtVsFeb_{nullptr};
  TH1F* h_dtOutOfRangePerFeb_{nullptr};
  TH1F* h_dtOutOfRangePerFebLastEwt_{nullptr};
  std::map<std::pair<int, uint8_t>, TH1F*> h1_dtFpgaPairs_;
  std::map<uint8_t, TGraph*> g_ubStatusVsEwt_;
  std::map<uint8_t, uint32_t> lastMicroBunchStatus_;

  bool warnedOffAxisFeb_{false};
  int maxFebIdSeen_{-1};
  long long nFebIdOutOfAxis_{0};
  long long nDtOutOfRange_{0};
  double maxAbsDtSeen_{0.0};
  // FEBs whose dt was out of range in the event being filled.
  std::vector<int> dtOutOfRangeThisEvent_;

  std::size_t nEvents_{0};
  std::size_t nDigis_{0};
  std::set<int> activeFEBs_;
  std::set<uint8_t> activeROCs_;
  std::map<uint8_t, std::set<uint8_t>> rocFEBMap_;

  std::deque<std::pair<uint64_t, int>> ewtWindow_;
  long long ewtWindowSum_{0};
  std::deque<std::pair<uint64_t, std::vector<uint16_t>>> recentChannelHitsByEwt_;
  std::deque<std::pair<uint64_t, std::vector<int>>> recentDtByEwt_;

  long long avgBlockSum_{0};
  std::size_t avgBlockCount_{0};
  uint64_t avgBlockFirstEwt_{0};
  bool avgSeedsCleared_{false};
};

} // namespace mu2e

#endif /* DQMHelpers_inc_CRVDigiDQM_hh */
