#ifndef CRVDQM_inc_CRVDigiDQM_hh
#define CRVDQM_inc_CRVDigiDQM_hh
// CRV digi DQM client: occupancy, waveform and timing histograms, used by the
// otsdaq online monitor and the offline DQM modules alike. No GeometryService:
// the caller injects the sector map and the FEB topology.
//
// Every histogram is booked up front for the full CRVId address space and the
// CRVDQMRun1 configurations, on the fixed axes below, so every job writes the
// same set and files merge across run periods.
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVDQMRun1.hh"
#include "Offline/DQMHelpers/inc/DQMClient.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"

#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"

#include <cstdint>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace mu2e {

class CRVDigiDQM : public DQMClient {
public:
  static constexpr int kBinningVersion = 1;

  // Axes. Each quantity whose range differs by orders of magnitude between
  // on-spill Run 1 and the ~400 us KPP readout window has a short view and a
  // full view ("2").
  static constexpr DQMAxis kDigisPerEvt = DQMAxis::Counts(0, 499);
  static constexpr DQMAxis kDigisPerEvt2{500, -0.5, 4999.5};
  static constexpr DQMAxis kPeakAdc{512, 0., 4096.};
  static constexpr DQMAxis kAdcSample{128, 0., 4096.};
  static constexpr DQMAxis kTdc{200, 0., 200.};
  static constexpr DQMAxis kTdc2{1024, 0., 40960.};
  static constexpr DQMAxis kOnlineChannel = DQMAxis::Counts(0, CRVDQMRun1::kNOnlineChannels - 1);
  static constexpr DQMAxis kFebChannel = DQMAxis::Counts(0, CRVDQMRun1::kNChanPerFEB - 1);
  static constexpr DQMAxis kFebPort = DQMAxis::Counts(0, CRVDQMRun1::kNFebPorts - 1);
  static constexpr DQMAxis kRocChannelEdges{CRVDQMRun1::kNChanPerROC, 0., CRVDQMRun1::kNChanPerROC};
  static constexpr DQMAxis kFebChannelEdges{CRVDQMRun1::kNChanPerFEB, 0., CRVDQMRun1::kNChanPerFEB};
  static constexpr DQMAxis kFebPortEdges{CRVDQMRun1::kNFebPorts, 0., CRVDQMRun1::kNFebPorts};
  static constexpr DQMAxis kOfflineChannel = DQMAxis::Counts(0, CRVDQMRun1::kNOfflineChannels - 1);
  static constexpr DQMAxis kBarId{200, -0.5, CRVId::nBars - 0.5};
  static constexpr DQMAxis kSiPM = DQMAxis::Counts(0, 3);
  static constexpr DQMAxis kFpgaPair =
      DQMAxis::Counts(0, CRVDQMRun1::kNFebPorts * CRVDQMRun1::kNFpgaPairs - 1);
  static constexpr DQMAxis kDtFpga = DQMAxis::Symmetric(50., 0.5);
  static constexpr DQMAxis kDigisPerChannelAndEvent{200, 0., 0.1};
  static constexpr DQMAxis kDigisPerChannelAndEvent2{250, 0., 5.};

  // Partner-FEB timing. One FEB reads two of the four layers on one side of one
  // module, so a muon fires a handful of geometrically related FEBs: these four
  // relations. dt is measured only inside a local coincidence group, never
  // against the detector as a whole, and never across sectors -- a muon going
  // in one side and out the other gives two genuinely separated traversals.
  enum DtClass { kDtSameModuleSameSide = 0, kDtSameModuleOtherSide,
                 kDtAdjModuleSameSide, kDtAdjModuleOtherSide, kNDtClasses };
  static const char* dtClassName(int c);
  // The same-side layer pair is a few-ns core; the cross-side classes carry the
  // light-propagation spread along the bar (~20 ns).
  static constexpr DQMAxis kDtPartner[kNDtClasses] = {
      DQMAxis::Symmetric(50., 0.5), DQMAxis::Symmetric(100., 1.),
      DQMAxis::Symmetric(100., 1.), DQMAxis::Symmetric(100., 1.)};
  static constexpr DQMAxis kLayersPerGroup = DQMAxis::Counts(0, 4);
  // Air showers put many muons through the whole CRV at once: real physics,
  // and what the full view is for.
  static constexpr DQMAxis kGroupsPerEvent = DQMAxis::Counts(0, 49);
  static constexpr DQMAxis kGroupsPerEvent2{200, -0.5, 999.5};
  static constexpr DQMAxis kSectorsPerEvent = DQMAxis::Counts(0, 23);

  // Live series: digis summed over the last kEwtWindow events, and the mean
  // over blocks of kAvgBlockSize events.
  static constexpr std::size_t kEwtWindow = 1000;
  static constexpr std::size_t kAvgBlockSize = 30;
  static constexpr std::size_t kAvgGraphPoints = 1000;

  // Which sector, module and side a FEB reads, indexed by FEB port.
  struct FebTopology {
    int sector{-1};
    int module{-1};
    int side{-1};
    bool valid{false};
  };

  explicit CRVDigiDQM(const DQMHistSet::Config& hists = {});

  // Caller-injected layout; either may be called at any time after Book().
  // channelToSector: offline channel -> index into the configuration's
  // CRVDQMRun1 sector list, -1 to skip the channel.
  void SetConfiguration(int configuration, const std::vector<int>& channelToSector);
  // febTopology indexed by FEB port; channelToLayer by offline channel.
  // Without it the partner-dt histograms stay empty.
  void SetFebTopology(const std::vector<FebTopology>& febTopology,
                      const std::vector<int>& channelToLayer);

  void Fill(const CrvDigiCollection& crvDigis, const CrvStatusCollection& crvStatus);

  TH1F* h1_digisPerEvt() const { return h1_digisPerEvt_; }  //digis per event
  TH1F* h1_digisPerEvt2() const { return h1_digisPerEvt2_; }  //same, full view
  TH1F* h1_peakAdc() const { return h1_peakAdc_; }  //largest ADC sample of a digi
  TH1F* h1_tdc() const { return h1_tdc_; }  //digi start time in 12.5 ns ticks
  TH1F* h1_tdc2() const { return h1_tdc2_; }  //same, full readout window
  TH1F* h1_channels() const { return h1_channels_; }  //occupancy vs online channel
  TH2F* h2_channels() const { return h2_channels_; }  //FEB port vs FEB channel
  TH1D* BarId() const { return hBarId_; }  //ValCrvDigi: scintillator bar index
  TH1D* SiPM() const { return hSiPM_; }  //ValCrvDigi: SiPM number within the bar
  TH1D* ADC() const { return hADC_; }  //ValCrvDigi: every ADC sample

  //CRVId occupancy maps (raw counts; divide by nEvents after hadd)
  TH1F* crvDigisPerChannel() const { return h_crvDigisPerChannel_; }  //vs offline channel
  TH2F* crvDigiRates() const { return h_crvDigiRates_; }  //FEB channel vs FEB port
  const std::vector<DQMH1<TH1F>>& crvDigiRatesROC() const { return h_crvDigiRatesROC_; }

  //intra-FEB timing: x = febPort*kNFpgaPairs + fpgaPairIndex, y = dt
  TH2F* dtFpgaPairs() const { return h2_dtFpgaPairs_; }
  //partner-FEB dt in a local coincidence group: a slipped FEB is a displaced column
  TH2F* dtPartner(int dtClass) const;
  TH1F* layersPerGroup() const { return h_layersPerGroup_; }
  TH1F* groupsPerEvent() const { return h_groupsPerEvent_; }
  TH1F* sectorsPerEvent() const { return h_sectorsPerEvent_; }
  TH1F* febNoGroup() const { return h_febNoGroup_; }  //FEB had hits, no group formed

  std::size_t nDigis() const { return nDigis_; }
  std::size_t nGroups() const { return nGroups_; }
  bool hasEwtWindow() const { return !ewtWindow_.empty(); }
  uint64_t lastEwt() const { return ewtWindow_.empty() ? 0 : ewtWindow_.back().first; }
  const std::set<int>& activeFebPorts() const { return activeFebPorts_; }
  const std::set<uint8_t>& activeROCs() const { return activeROCs_; }
  const std::map<uint8_t, std::set<uint8_t>>& rocFEBMap() const { return rocFEBMap_; }

private:
  struct FpgaHit {
    double time_ns;
    uint8_t channel;
  };
  // One above-threshold hit, with the geography the grouping needs.
  struct PartnerHit {
    double time_ns{0.};
    int febPort{-1};
    int sector{-1};
    int module{-1};
    int side{-1};
    int layer{-1};
  };

  void book() override;
  void endJob() override;
  void resetForNewRun() override;

  void fillFpgaTiming(const std::map<int, std::map<uint8_t, std::vector<FpgaHit>>>& hitTimes);
  void fillPartnerTiming(std::vector<PartnerHit>& hits);
  void fillGroup(const std::vector<const PartnerHit*>& group);
  int dtClassFor(int febA, int febB) const;
  void fillEwtSeries(uint64_t ewt, int nDigis);
  void fillSectorOccupancy();

  DQMH1<TH1F> h1_digisPerEvt_;
  DQMH1<TH1F> h1_digisPerEvt2_;
  DQMH1<TH1F> h1_peakAdc_;
  DQMH1<TH1F> h1_tdc_;
  DQMH1<TH1F> h1_tdc2_;
  DQMH1<TH1F> h1_channels_;
  DQMH2<TH2F> h2_channels_;
  DQMH1<TH1D> hBarId_;
  DQMH1<TH1D> hSiPM_;
  DQMH1<TH1D> hADC_;
  DQMH1<TH1F> h_crvDigisPerChannel_;
  DQMH2<TH2F> h_crvDigiRates_;
  std::vector<DQMH1<TH1F>> h_crvDigiRatesROC_;
  DQMH2<TH2F> h2_dtFpgaPairs_;
  std::vector<DQMH2<TH2F>> h2_dtPartner_;
  DQMH1<TH1F> h_layersPerGroup_;
  DQMH1<TH1F> h_groupsPerEvent_;
  DQMH1<TH1F> h_groupsPerEvent2_;
  DQMH1<TH1F> h_sectorsPerEvent_;
  DQMH1<TH1F> h_febNoGroup_;
  // [configuration][sector], short and full view
  std::vector<std::vector<DQMH1<TH1F>>> h_sectorOccupancy_;
  std::vector<std::vector<DQMH1<TH1F>>> h_sectorOccupancy2_;

  DQMSeries* g_digisVsEwt_{nullptr};
  DQMSeries* g_digisAvgVsEwt_{nullptr};

  int configuration_{-1};
  std::vector<int> channelToSector_;
  std::vector<FebTopology> febTopology_;
  std::vector<int> channelToLayer_;
  std::vector<int> nDigisOffline_;

  std::size_t nDigis_{0};
  std::size_t nGroups_{0};
  std::set<int> activeFebPorts_;
  std::set<uint8_t> activeROCs_;
  std::map<uint8_t, std::set<uint8_t>> rocFEBMap_;

  std::deque<std::pair<uint64_t, int>> ewtWindow_;
  long long ewtWindowSum_{0};
  long long avgBlockSum_{0};
  std::size_t avgBlockCount_{0};
  uint64_t avgBlockFirstEwt_{0};
};

} // namespace mu2e

#endif /* CRVDQM_inc_CRVDigiDQM_hh */
