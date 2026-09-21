#ifndef CRVDQM_inc_CRVRecoDQM_hh
#define CRVDQM_inc_CRVRecoDQM_hh
// CRV reco DQM client: coincidence-cluster and reco-pulse histograms, and the
// PE spectrum of every channel's cluster pulses as one mergeable TH2, reduced
// at end of job to Landau(x)Gauss MPV maps. No GeometryService, no Proditions:
// the caller injects the configuration and sector map.
//
// The MPV maps are fit results and do not merge. After hadd, refit the merged
// crvPEsVsChannel with FitPEsVsChannel().
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVDQMRun1.hh"
#include "Offline/DQMHelpers/inc/DQMClient.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"

#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace mu2e {

class CRVRecoDQM : public DQMClient {
public:
  static constexpr int kBinningVersion = 1;

  static constexpr DQMAxis kOne{1, 0.5, 1.5};
  static constexpr DQMAxis kSectorType{10, 0., 10.};
  static constexpr DQMAxis kPE{75, 0., 150.};
  static constexpr DQMAxis kOfflineChannel =
      DQMAxis::Counts(0, CRVDQMRun1::kNOfflineChannels - 1);
  static constexpr DQMAxis kRocChannelEdges{CRVDQMRun1::kNChanPerROC, 0., CRVDQMRun1::kNChanPerROC};
  static constexpr DQMAxis kFebChannelEdges{CRVDQMRun1::kNChanPerFEB, 0., CRVDQMRun1::kNChanPerFEB};
  static constexpr DQMAxis kFebPortEdges{CRVDQMRun1::kNFebPorts, 0., CRVDQMRun1::kNFebPorts};
  // DqmCrv / ValCrv* plots. "2" is the full-readout-window view: KPP windows
  // reach ~410 us.
  static constexpr DQMAxis kNPulses{101, -0.5, 100.5};
  static constexpr DQMAxis kNPulses2{200, -0.5, 3999.5};
  static constexpr DQMAxis kBarId{200, -0.5, CRVId::nBars - 0.5};
  static constexpr DQMAxis kSiPM = DQMAxis::Counts(0, 3);
  static constexpr DQMAxis kPulsePE{100, 0., 400.};
  static constexpr DQMAxis kTime{100, 0., 2000.};
  static constexpr DQMAxis kTime2{500, 0., 500.e3};
  static constexpr DQMAxis kChi2{100, 0., 20.};
  static constexpr DQMAxis kLogChi2{100, -3., 5.};
  static constexpr DQMAxis kNClusters{101, -0.5, 100.5};
  static constexpr DQMAxis kClusterPE{250, 0., 5000.};

  // MPV fit: window as a fraction of the spectrum peak, its hard low edge, and
  // the least a spectrum must hold to be fitted.
  static constexpr double kFitRangeStart = 0.7;
  static constexpr double kFitRangeEnd = 2.0;
  static constexpr double kFitMinPE = 15.0;
  static constexpr double kFitMinBinCenter = 10.0;
  static constexpr double kFitMinPeakEntries = 20.0;

  struct LandauGaussResult {
    float mpv{0.f};
    float chi2PerNdf{std::numeric_limits<float>::quiet_NaN()};
    bool fitted{false};  //false when the spectrum was too sparse to fit
  };
  // One channel's spectrum.
  static LandauGaussResult FitSpectrum(TH1* h);
  // Every offline channel of crvPEsVsChannel, merged or not; an empty channel
  // gives an unfitted result.
  static std::vector<LandauGaussResult> FitPEsVsChannel(const TH2& peVsChannel);

  explicit CRVRecoDQM(const DQMHistSet::Config& hists = {});

  // Caller-injected layout, callable any time after Book(). channelToSector:
  // offline channel -> index into the configuration's CRVDQMRun1 sector list,
  // -1 to skip it (how a caller drops notConnected channels).
  void SetConfiguration(int configuration, const std::vector<int>& channelToSector);

  void Fill(const CrvCoincidenceClusterCollection& clusters);
  // The reco-pulse plots need every pulse, not just the ones a coincidence
  // kept, so they are only filled through this overload.
  void Fill(const CrvCoincidenceClusterCollection& clusters,
            const CrvRecoPulseCollection& recoPulses);

  TH1F* nEventsWithClustersHist() const { return h_nEventsWithClusters_; }
  TH1I* coincidenceClusters() const { return h_coincidenceClusters_; }  //by CrvSectorType
  TH2F* PEsVsChannel() const { return h_PEsVsChannel_; }  //the mergeable fit input
  TH2F* PEsMPV() const { return h_PEsMPV_; }  //FEB channel vs FEB port
  TH1F* PEsMPVROC(int roc) const;  //online channel within ROC, roc 1-based
  TH1F* PEsMPVSector(int configuration, int sector) const;

  TH1D* NPulses() const { return h_NPulses_; }
  TH1D* NPulse2() const { return h_NPulse2_; }
  TH1D* BarIdr() const { return h_BarIdr_; }
  TH1D* SiPMr() const { return h_SiPMr_; }
  TH1D* PEr() const { return h_PEr_; }
  TH1D* PEHeight() const { return h_PEHeight_; }
  TH1D* PulseTime() const { return h_PulseTime_; }
  TH1D* PulseTime2() const { return h_PulseTime2_; }
  TH1D* chi2() const { return h_chi2_; }
  TH1D* logchi2() const { return h_logchi2_; }
  TH1D* LeadingTime() const { return h_LeadingTime_; }
  TH1D* LeadingTime2() const { return h_LeadingTime2_; }
  TH1D* NClus() const { return h_NClus_; }
  TH1D* NPc() const { return h_NPc_; }
  TH1D* PEc() const { return h_PEc_; }
  TH1D* tc() const { return h_tc_; }
  TH1D* t2c() const { return h_t2c_; }
  // Cluster position, one set per configuration: X(c), Y(c), Z(c).
  TH1D* X(int configuration) const { return position(configuration, 0); }
  TH1D* Y(int configuration) const { return position(configuration, 1); }
  TH1D* Z(int configuration) const { return position(configuration, 2); }

  std::size_t nEventsWithClusters() const { return nEventsWithClusters_; }
  std::size_t nClusters() const { return nClusters_; }
  std::size_t nRecoPulses() const { return nRecoPulses_; }  //cluster members only
  std::size_t nFits() const { return nFits_; }
  std::size_t nFitsSucceeded() const { return nFitsSucceeded_; }

private:
  void book() override;
  void endJob() override;

  TH1D* position(int configuration, int coordinate) const;
  void fillClusters(const CrvCoincidenceClusterCollection& clusters);
  void fillRecoPulses(const CrvRecoPulseCollection& recoPulses);

  DQMH1<TH1F> h_nEventsWithClusters_;
  DQMH1<TH1I> h_coincidenceClusters_;
  DQMH2<TH2F> h_PEsVsChannel_;
  DQMH2<TH2F> h_PEsMPV_;
  std::vector<DQMH1<TH1F>> h_PEsMPVROC_;
  std::vector<std::vector<DQMH1<TH1F>>> h_PEsMPVSector_;  //[configuration][sector]

  DQMH1<TH1D> h_NPulses_;
  DQMH1<TH1D> h_NPulse2_;
  DQMH1<TH1D> h_BarIdr_;
  DQMH1<TH1D> h_SiPMr_;
  DQMH1<TH1D> h_PEr_;
  DQMH1<TH1D> h_PEHeight_;
  DQMH1<TH1D> h_PulseTime_;
  DQMH1<TH1D> h_PulseTime2_;
  DQMH1<TH1D> h_chi2_;
  DQMH1<TH1D> h_logchi2_;
  DQMH1<TH1D> h_LeadingTime_;
  DQMH1<TH1D> h_LeadingTime2_;
  DQMH1<TH1D> h_NClus_;
  DQMH1<TH1D> h_NPc_;
  DQMH1<TH1D> h_PEc_;
  DQMH1<TH1D> h_tc_;
  DQMH1<TH1D> h_t2c_;
  std::vector<std::vector<DQMH1<TH1D>>> h_position_;  //[configuration][x,y,z]

  int configuration_{-1};
  std::vector<int> channelToSector_;
  // Online channel of each offline channel, learned from the pulses themselves,
  // so the online MPV maps need no channel map. -1 until seen.
  std::vector<int> offlineToOnline_;

  std::size_t nEventsWithClusters_{0};
  std::size_t nClusters_{0};
  std::size_t nRecoPulses_{0};
  std::size_t nFits_{0};
  std::size_t nFitsSucceeded_{0};
};

} // namespace mu2e

#endif /* CRVDQM_inc_CRVRecoDQM_hh */
