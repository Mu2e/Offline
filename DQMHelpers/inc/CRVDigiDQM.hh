#ifndef DQMHelpers_inc_CRVDigiDQM_hh
#define DQMHelpers_inc_CRVDigiDQM_hh
// Standalone CRV digi DQM helper: books and fills the histograms used by the
// otsdaq online monitor and the offline DQM modules. No GeometryService.
// Original Author: R. Mina

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "art_root_io/TFileDirectory.h"

#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

#include <cstdint>
#include <string>
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
    // CRVId-indexed occupancy maps (no GeometryService). Stored as raw
    // counts; divide by nEvents after hadd.
    bool fillCrvIdRates{true};
    // KPP FEB-axis sizing: book the 33-slot occupancy pair and size the FEB
    // axes for ROC 1-2. False for full CRV; see the header comment above.
    bool kppReadout{true};
    // TGraphs vs EWT and *LastEwt snapshots. Online monitor only — they do
    // not survive hadd. Offline analyzers leave this false.
    bool fillLivePlots{false};
  };

  // KPP readout geography used by the online occupancy plots.
  // 25-slot stride is the online convention (not CRVId::nFEBPerROC = 24).
  static constexpr int kNFebSlotsPerROC = 25;
  static constexpr int kNGlobalFebBins = 32;
  static constexpr int kNGlobalChannelBins =
      (kNGlobalFebBins + 1) * static_cast<int>(CRVId::nChanPerFEB);

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

  static int globalFebId(uint8_t roc, uint8_t feb);
  static int globalChannelId(uint8_t roc, uint8_t feb, uint8_t febChannel);

  explicit CRVDigiDQM(const Config& config);

  void Book(art::TFileDirectory dir);
  // Caller injects the geometry-derived sector map so the helper needs no
  // GeometryService. Negative sector skips the channel. Call once after Book().
  void BookSectorOccupancy(const std::vector<std::string>& sectorNames,
                           const std::vector<int>& channelToSector,
                           int nBins, double lo, double hi);
  //digis/event per channel, one hist per CRV sector
  const std::vector<TH1F*>& sectorOccupancy() const { return h_sectorOccupancy_; }
  void Fill(const CrvDigiCollection& crvDigis,
            const CrvStatusCollection& crvStatus);
  void WriteGraphs();

  TH1F* h1_digisPerEvt() const { return h1_digisPerEvt_; }  //digis per event
  TH1F* nEventsHist() const { return h_nEvents_; }  //one count per event; hadd-safe
  TH1F* h1_peakAdc() const { return h1_peakAdc_; }  //largest ADC sample of a digi
  TH1F* h1_tdc() const { return h1_tdc_; }  //digi start time in 12.5ns ticks
  TH1F* h1_channels() const { return h1_channels_; }  //occupancy vs global channel ID
  TH1F* h1_channelsLastEwt() const { return h1_channelsLastEwt_; }  //same, rolling EWT window
  TH2F* h2_channels() const { return h2_channels_; }  //FEB vs FEB channel hit map
  TGraph* g_digisVsEwt() const { return g_digisVsEwt_; }  //digis in the last EWTs vs EWT
  TGraph* g_digisAvgVsEwt() const { return g_digisAvgVsEwt_; }  //mean digis/event vs EWT

  TH1D* BarId() const { return hBarId_; }  //ValCrvDigi: scintillator bar index
  TH1D* SiPM() const { return hSiPM_; }  //ValCrvDigi: SiPM number within the bar
  TH1D* ADC() const { return hADC_; }  //ValCrvDigi: every ADC sample

  //CRVId occupancy maps (raw counts; divide by nEvents after hadd)
  TH1F* crvDigisPerChannel() const { return h_crvDigisPerChannel_; }  //vs offline channel bar*4+SiPM
  TH2F* crvDigiRates() const { return h_crvDigiRates_; }  //FEB channel vs FEB port
  const std::vector<TH1F*>& crvDigiRatesROC() const { return h_crvDigiRatesROC_; }  //one per ROC
  int nDigisOffline(std::size_t channel) const;
  std::size_t nOfflineChannels() const { return nDigisOffline_.size(); }

  //inter-FEB sync: first-hit time minus the median of the other FEBs, per FEB
  TH2F* dtVsFeb() const { return h2_dtVsFeb_; }

  //a FEB whose clock has slipped stands above the others in these
  TH1F* dtOutOfRangePerFeb() const { return h_dtOutOfRangePerFeb_; }  //events with |dt| off scale
  TH1F* dtOutOfRangePerFebLastEwt() const  //same, rolling EWT window
  {
    return h_dtOutOfRangePerFebLastEwt_;
  }
  //intra-FEB timing, keyed (globalFebId, fpgaA*nFPGAPerFEB+fpgaB)
  const std::map<std::pair<int, uint8_t>, TH1F*>& dtFpgaPairs() const
  {
    return h1_dtFpgaPairs_;
  }
  //ROC MicroBunchStatus vs EWT, one graph per DTC link
  const std::map<uint8_t, TGraph*>& ubStatusVsEwt() const
  {
    return g_ubStatusVsEwt_;
  }

  // Axis-coverage diagnostics. Off-axis FEB and CRVId-range: one LogWarning
  // from Fill. maxFebIdSeen / maxAbsDtSeen: one LogWarning from WriteGraphs
  // if they sit past the hist axes (overflow bins do not store those values).
  int maxFebIdSeen() const { return maxFebIdSeen_; }
  long long nFebIdOutOfAxis() const { return nFebIdOutOfAxis_; }
  long long nCrvIdOutOfRange() const { return nCrvIdOutOfRange_; }
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
  void fillRollingDtOutOfRange(uint64_t ewt);
  void fillSectorOccupancy();

  Config config_;
  int nBinsDt_{400};
  int nBinsDtVsFeb_{500};
  bool booked_{false};

  std::optional<art::TFileDirectory> dir_;
  std::optional<art::TFileDirectory> timingFpgaDir_;

  TH1F* h1_digisPerEvt_{nullptr};
  TH1F* h_nEvents_{nullptr};
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

  TH2F* h2_dtVsFeb_{nullptr};
  std::vector<TH1F*> h_sectorOccupancy_;
  std::vector<int> channelToSector_;
  TH1F* h_dtOutOfRangePerFeb_{nullptr};
  TH1F* h_dtOutOfRangePerFebLastEwt_{nullptr};
  std::map<std::pair<int, uint8_t>, TH1F*> h1_dtFpgaPairs_;
  std::map<uint8_t, TGraph*> g_ubStatusVsEwt_;
  std::map<uint8_t, uint32_t> lastMicroBunchStatus_;

  bool warnedOffAxisFeb_{false};
  bool warnedCrvIdOutOfRange_{false};
  int maxFebIdSeen_{-1};
  long long nFebIdOutOfAxis_{0};
  long long nCrvIdOutOfRange_{0};
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
