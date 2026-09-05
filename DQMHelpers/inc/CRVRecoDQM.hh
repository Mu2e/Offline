#ifndef DQMHelpers_inc_CRVRecoDQM_hh
#define DQMHelpers_inc_CRVRecoDQM_hh
// Standalone CRV reco DQM helper: the coincidence-cluster sector-type spectrum
// and the per-channel PE spectra of the pulses in those clusters, reduced at
// end of job to Landau(x)Gauss MPV maps. No GeometryService, no Proditions.
// Original Author: R. Mina

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"

#include "art_root_io/TFileDirectory.h"

#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace mu2e {

class CRVRecoDQM {
public:
  // Binning defaults match the CrvDQMcollector FHiCL defaults, so a cutover on
  // the same file is comparable.
  struct Config {
    // per-channel PE spectra: the input to the Landau(x)Gauss fit
    int nBinsPEs{75};
    double minPEs{0.0};
    double maxPEs{150.0};
    // fit window as a fraction of the spectrum peak, and its hard low edge
    double PEfitRangeStart{0.7};
    double PEfitRangeEnd{2.0};
    double PEstart{15.0};
    // CrvSectorType axis of crvCoincidencesClusters
    int nSectorTypeBins{10};
    // ~50k per-channel spectra. Expert output; off by default.
    bool writePerChannelPE{false};
    // The per-event DqmCrv / ValCrv* plots.
    bool fillInclusive{true};

    // Axes of the inclusive plots that depend on where the detector is and how
    // long its readout window is. Defaults are the DqmCrv / ValCrv* values, cut
    // for the full CRV in Mu2e coordinates; on the extracted CRV the y and z
    // defaults are entirely overflow (it sits at y ~ 4237-4649, z ~ 21440-23065),
    // which is what these exist to fix. Everything else in the inclusive block
    // is detector-independent and stays fixed, so those series cannot drift.
    //
    // Each time quantity is booked twice, as a short view and a full-window
    // view. The three short axes share one definition and the three long axes
    // another, because they are the same quantity in the same window and the
    // point of them is that pulse, leading-edge and cluster times overlay.
    int nBinsTime{100};      // PulseTime, LeadingTime, tc
    double minTime{0.0};
    double maxTime{2000.0};
    int nBinsTime2{100};     // PulseTime2, LeadingTime2, t2c
    double minTime2{0.0};
    double maxTime2{100.0e3};
    int nBinsPos{100};       // X, Y, Z share a granularity, not a range
    double minX{-3904 - 3000};
    double maxX{-3904 + 3000};
    double minY{0.0};
    double maxY{3000.0};
    double minZ{-3500.0};
    double maxZ{20000.0};
  };

  // Online geography. The MPV index and the MPV axes are both derived from
  // these, so the two cannot drift apart.
  static constexpr int kNChanPerFEB = static_cast<int>(CRVId::nChanPerFEB);
  static constexpr int kNFebPerROC = static_cast<int>(CRVId::nFEBPerROC);
  static constexpr int kNROC = static_cast<int>(CRVId::nROC);
  static constexpr int kNChannelsPerROC = kNFebPerROC * kNChanPerFEB;
  static constexpr int kNFebPorts = kNROC * kNFebPerROC;
  static constexpr int kNOnlineChannels = kNROC * kNChannelsPerROC;

  // ROC and FEB are 1-based on the wire; FEB channel is 0-based.
  static bool onlineIdInRange(int roc, int feb, int febChannel);
  static int rocChannel(int feb, int febChannel);
  static int febPort(int roc, int feb);
  static int onlineChannelIndex(int roc, int feb, int febChannel);

  struct LandauGaussResult {
    float mpv{0.f};
    float chi2PerNdf{std::numeric_limits<float>::quiet_NaN()};
    bool fitted{false};  //false when the spectrum was too sparse to fit
  };

  explicit CRVRecoDQM(const Config& config);

  void Book(art::TFileDirectory dir);
  // Caller injects the geometry-derived sector map so the helper needs no
  // GeometryService. Negative sector skips the channel, which is how a caller
  // drops its Proditions notConnected channels. Call once after Book().
  void BookSectorMPV(const std::vector<std::string>& sectorNames,
                     const std::vector<int>& channelToSector);
  // X/Y/Z are booked apart from Book() because their ranges depend on where the
  // detector is, which Book() (beginJob) cannot know. A caller holding geometry
  // injects the CRV envelope as a min/max corner pair in Mu2e coordinates --
  // built from the counters, so it works on a countersOnly geometry -- and the
  // helper keeps no GeometryService dependency of its own. Optional: a caller
  // that never calls either overload gets the Config ranges booked on first
  // fill. An empty or degenerate envelope falls back to Config too.
  void BookPositionAxes(const std::vector<double>& minPoint,
                        const std::vector<double>& maxPoint,
                        double margin = 0.05);
  void BookPositionAxes();  //from Config
  void Fill(const CrvCoincidenceClusterCollection& clusters);
  // The inclusive reco-pulse plots need every pulse, not just the ones a
  // coincidence kept, so they are only filled through this overload.
  void Fill(const CrvCoincidenceClusterCollection& clusters,
            const CrvRecoPulseCollection& recoPulses);
  // endJob hook, named for symmetry with CRVDigiDQM / CRVStatusDQM. Runs the
  // per-channel fits and fills the MPV maps; there is nothing to fill before it.
  void WriteGraphs();

  TH1F* nEventsHist() const { return h_nEvents_; }  //one count per event; hadd-safe
  TH1F* nEventsWithClustersHist() const { return h_nEventsWithClusters_; }  //denominator for the cluster rate
  TH1I* coincidenceClusters() const { return h_coincidenceClusters_; }  //clusters by CrvSectorType

  //Landau(x)Gauss MPV of each channel's PE spectrum, one hist per CRV sector
  const std::vector<TH1F*>& PEsMPVSector() const { return h_PEsMPVSector_; }
  //the same MPVs indexed by online channel within a ROC, one hist per ROC
  const std::vector<TH1F*>& PEsMPVROC() const { return h_PEsMPVROC_; }
  //and as a whole-detector FEB-channel vs FEB-port map
  TH2F* PEsMPV() const { return h_PEsMPV_; }

  //per-channel PE spectra, indexed bar*nChanPerBar+SiPM and by online channel
  TH1F* PEsOffline(std::size_t channel) const;
  TH1F* PEsOnline(std::size_t onlineChannel) const;

  //DqmCrv reco-pulse block, filled from the whole CrvRecoPulseCollection
  TH1D* NPulses() const { return h_NPulses_; }  //reco pulses per event, 0-100
  TH1D* NPulse2() const { return h_NPulse2_; }  //same, 0-2000
  TH1D* BarIdr() const { return h_BarIdr_; }  //reco-pulse occupancy vs bar index
  TH1D* SiPMr() const { return h_SiPMr_; }  //which of the 4 SiPM positions
  TH1D* PEr() const { return h_PEr_; }  //PE from the pulse-area fit
  TH1D* PEHeight() const { return h_PEHeight_; }  //PE from the pulse height
  TH1D* PulseTime() const { return h_PulseTime_; }  //fitted peak time, 0-2000 ns
  TH1D* PulseTime2() const { return h_PulseTime2_; }  //same, 0-100 us
  TH1D* chi2() const { return h_chi2_; }  //pulse-fit chi2
  TH1D* logchi2() const { return h_logchi2_; }  //log10 of the same
  TH1D* LeadingTime() const { return h_LeadingTime_; }  //LE time, 0-2000 ns
  TH1D* LeadingTime2() const { return h_LeadingTime2_; }  //same, 0-100 us

  //DqmCrv coincidence-cluster block
  TH1D* NClus() const { return h_NClus_; }  //clusters per event
  TH1D* NPc() const { return h_NPc_; }  //reco pulses per cluster
  TH1D* PEc() const { return h_PEc_; }  //total PE per cluster
  TH1D* tc() const { return h_tc_; }  //cluster start time, 0-2000 ns
  TH1D* t2c() const { return h_t2c_; }  //same, 0-100 us
  TH1D* X() const { return h_X_; }  //cluster average hit position
  TH1D* Y() const { return h_Y_; }
  TH1D* Z() const { return h_Z_; }

  std::size_t nEvents() const { return nEvents_; }
  std::size_t nEventsWithClusters() const { return nEventsWithClusters_; }
  std::size_t nClusters() const { return nClusters_; }
  std::size_t nRecoPulses() const { return nRecoPulses_; }  //cluster members only
  std::size_t nInclusiveRecoPulses() const { return nInclusiveRecoPulses_; }  //whole collection

  // Coverage diagnostics, read rather than logged per event. Out-of-range
  // ROC/FEB and null Ptrs raise one LogWarning each from Fill.
  long long nOnlineIdOutOfRange() const { return nOnlineIdOutOfRange_; }
  long long nOfflineChannelOutOfRange() const { return nOfflineChannelOutOfRange_; }
  long long nNullRecoPulsePtrs() const { return nNullRecoPulsePtrs_; }

  // Fit quality. An MPV with no goodness of fit next to it is not actionable,
  // and a channel that never fitted is not the same as one whose MPV is zero.
  std::size_t nFits() const { return nFits_; }
  std::size_t nFitsSucceeded() const { return nFitsSucceeded_; }
  double meanFitChi2PerNdf() const;

  bool booked() const { return booked_; }

private:
  LandauGaussResult fitChannel(TH1F* h);
  TH1F* spectrum(std::vector<TH1F*>& hists, std::size_t index,
                 const char* namePrefix);
  void bookInclusive(art::TFileDirectory& dir);
  void bookPositionHists(const double lo[3], const double hi[3]);
  void ensurePositionAxes();
  void fillClusters(const CrvCoincidenceClusterCollection& clusters);
  void fillInclusiveRecoPulses(const CrvRecoPulseCollection& recoPulses);
  void fitSectorMPV();
  void fitOnlineMPV();

  Config config_;
  bool booked_{false};

  std::optional<art::TFileDirectory> dir_;
  std::optional<art::TFileDirectory> spectraDir_;

  TH1F* h_nEvents_{nullptr};
  TH1F* h_nEventsWithClusters_{nullptr};
  TH1I* h_coincidenceClusters_{nullptr};

  std::vector<TH1F*> h_PEsMPVSector_;
  std::vector<TH1F*> h_PEsMPVROC_;
  TH2F* h_PEsMPV_{nullptr};

  // Booked on first fill, so a KPP-sized geometry pays for its own channels
  // rather than for all of CRVId. Views into either ownedSpectra_ (default) or
  // the output file's perChannelPE directory (Config::writePerChannelPE).
  std::vector<TH1F*> h_PEsOffline_;
  std::vector<TH1F*> h_PEsOnline_;
  std::vector<std::unique_ptr<TH1F>> ownedSpectra_;

  std::vector<int> channelToSector_;

  TH1D* h_NPulses_{nullptr};
  TH1D* h_NPulse2_{nullptr};
  TH1D* h_BarIdr_{nullptr};
  TH1D* h_SiPMr_{nullptr};
  TH1D* h_PEr_{nullptr};
  TH1D* h_PEHeight_{nullptr};
  TH1D* h_PulseTime_{nullptr};
  TH1D* h_PulseTime2_{nullptr};
  TH1D* h_chi2_{nullptr};
  TH1D* h_logchi2_{nullptr};
  TH1D* h_LeadingTime_{nullptr};
  TH1D* h_LeadingTime2_{nullptr};

  TH1D* h_NClus_{nullptr};
  TH1D* h_NPc_{nullptr};
  TH1D* h_PEc_{nullptr};
  TH1D* h_tc_{nullptr};
  TH1D* h_t2c_{nullptr};
  TH1D* h_X_{nullptr};
  TH1D* h_Y_{nullptr};
  TH1D* h_Z_{nullptr};

  std::size_t nEvents_{0};
  std::size_t nEventsWithClusters_{0};
  std::size_t nClusters_{0};
  std::size_t nRecoPulses_{0};
  std::size_t nInclusiveRecoPulses_{0};

  long long nOnlineIdOutOfRange_{0};
  long long nOfflineChannelOutOfRange_{0};
  long long nNullRecoPulsePtrs_{0};
  bool warnedOnlineIdOutOfRange_{false};
  bool warnedOfflineChannelOutOfRange_{false};
  bool warnedNullRecoPulsePtr_{false};

  std::size_t nFits_{0};
  std::size_t nFitsSucceeded_{0};
  double fitChi2Sum_{0.0};
  std::size_t fitChi2N_{0};
};

} // namespace mu2e

#endif /* DQMHelpers_inc_CRVRecoDQM_hh */
