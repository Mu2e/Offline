// Thin art analyzer that constructs, books, and fills CRVRecoDQM.
//
// Original Author: R. Mina

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/DQMHelpers/inc/CRVRecoDQM.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

class CRVRecoDQMAnalyzer : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<art::InputTag> crvCoincidenceClusterTag{
        Name("crvCoincidenceClusterTag"),
        Comment("CRV coincidence cluster finder"),
        art::InputTag{"CrvCoincidenceClusterFinder"}};
    fhicl::Atom<art::InputTag> crvRecoPulseTag{
        Name("crvRecoPulseTag"),
        Comment("CRV reco pulse producer; empty skips the inclusive pulse plots"),
        art::InputTag{"CrvRecoPulses"}};
    fhicl::Atom<std::string> outputTag{
        Name("outputTag"),
        Comment("TFileService subdirectory; empty books in the module directory"),
        ""};
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic level"), 0};

    fhicl::Atom<int> nBinsPEs{
        Name("nBinsPEs"), Comment("Bins for the PE spectra and the MPV hists"), 75};
    fhicl::Atom<double> minPEs{
        Name("minPEs"), Comment("Low edge for the PE spectra and MPV hists"), 0.0};
    fhicl::Atom<double> maxPEs{
        Name("maxPEs"), Comment("High edge for the PE spectra and MPV hists"), 150.0};
    fhicl::Atom<double> PEfitRangeStart{
        Name("PEfitRangeStart"),
        Comment("Low end of the PE MPV fit range as a fraction of the peak"),
        0.7};
    fhicl::Atom<double> PEfitRangeEnd{
        Name("PEfitRangeEnd"),
        Comment("High end of the PE MPV fit range as a fraction of the peak"),
        2.0};
    fhicl::Atom<double> PEstart{
        Name("PEstart"), Comment("Lowest PE for the fit"), 15.0};
    fhicl::Atom<int> nSectorTypeBins{
        Name("nSectorTypeBins"),
        Comment("Bins for crvCoincidencesClusters (CrvSectorType axis)"),
        10};
    fhicl::Atom<bool> writePerChannelPE{
        Name("writePerChannelPE"),
        Comment("Write the ~50k per-channel PE spectra (expert output)"),
        false};
    fhicl::Atom<bool> fillInclusive{
        Name("fillInclusive"),
        Comment("Also fill the per-event DqmCrv reco-pulse and cluster plots"),
        true};

    // Defaults are the DqmCrv values (full CRV, Mu2e coordinates). The
    // extracted CRV needs y/z moved or those two plots are all overflow.
    fhicl::Atom<int> nBinsTime{
        Name("nBinsTime"), Comment("Bins for PulseTime / LeadingTime / tc"), 100};
    fhicl::Atom<double> minTime{
        Name("minTime"), Comment("Low edge for PulseTime / LeadingTime / tc [ns]"), 0.0};
    fhicl::Atom<double> maxTime{
        Name("maxTime"), Comment("High edge for PulseTime / LeadingTime / tc [ns]"), 2000.0};
    fhicl::Atom<int> nBinsTime2{
        Name("nBinsTime2"),
        Comment("Bins for the full-window PulseTime2 / LeadingTime2 / t2c"), 100};
    fhicl::Atom<double> minTime2{
        Name("minTime2"),
        Comment("Low edge for PulseTime2 / LeadingTime2 / t2c [ns]"), 0.0};
    fhicl::Atom<double> maxTime2{
        Name("maxTime2"),
        Comment("High edge for PulseTime2 / LeadingTime2 / t2c [ns]"), 100000.0};
    fhicl::Atom<int> nBinsPos{
        Name("nBinsPos"), Comment("Bins for the cluster position plots X/Y/Z"), 100};
    fhicl::Atom<double> minX{Name("minX"), Comment("Low edge for X [mm]"), -6904.0};
    fhicl::Atom<double> maxX{Name("maxX"), Comment("High edge for X [mm]"), -904.0};
    fhicl::Atom<double> minY{Name("minY"), Comment("Low edge for Y [mm]"), 0.0};
    fhicl::Atom<double> maxY{Name("maxY"), Comment("High edge for Y [mm]"), 3000.0};
    fhicl::Atom<double> minZ{Name("minZ"), Comment("Low edge for Z [mm]"), -3500.0};
    fhicl::Atom<double> maxZ{Name("maxZ"), Comment("High edge for Z [mm]"), 20000.0};
    fhicl::Atom<bool> fillSectorMPV{
        Name("fillSectorMPV"),
        Comment("Fill per-sector crvPEsMPV_CRVsector* using GeometryService"),
        false};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit CRVRecoDQMAnalyzer(const Parameters& conf);

  void beginJob() override;
  void beginRun(const art::Run& run) override;
  void analyze(const art::Event& event) override;
  void endJob() override;

private:
  static CRVRecoDQM::Config makeHelperConfig(const Config& conf);

  art::InputTag crvCoincidenceClusterTag_;
  art::InputTag crvRecoPulseTag_;
  std::string outputTag_;
  int diagLevel_;
  bool fillSectorMPV_;
  CRVRecoDQM dqm_;
  bool sectorMapSent_{false};
  bool warnedMissingClusters_{false};
  bool warnedMissingPulses_{false};
};

CRVRecoDQM::Config CRVRecoDQMAnalyzer::makeHelperConfig(const Config& conf)
{
  CRVRecoDQM::Config c;
  c.nBinsPEs = conf.nBinsPEs();
  c.minPEs = conf.minPEs();
  c.maxPEs = conf.maxPEs();
  c.PEfitRangeStart = conf.PEfitRangeStart();
  c.PEfitRangeEnd = conf.PEfitRangeEnd();
  c.PEstart = conf.PEstart();
  c.nSectorTypeBins = conf.nSectorTypeBins();
  c.writePerChannelPE = conf.writePerChannelPE();
  c.fillInclusive = conf.fillInclusive();
  c.nBinsTime = conf.nBinsTime();
  c.minTime = conf.minTime();
  c.maxTime = conf.maxTime();
  c.nBinsTime2 = conf.nBinsTime2();
  c.minTime2 = conf.minTime2();
  c.maxTime2 = conf.maxTime2();
  c.nBinsPos = conf.nBinsPos();
  c.minX = conf.minX();
  c.maxX = conf.maxX();
  c.minY = conf.minY();
  c.maxY = conf.maxY();
  c.minZ = conf.minZ();
  c.maxZ = conf.maxZ();
  return c;
}

CRVRecoDQMAnalyzer::CRVRecoDQMAnalyzer(const Parameters& conf) :
    art::EDAnalyzer{conf},
    crvCoincidenceClusterTag_(conf().crvCoincidenceClusterTag()),
    crvRecoPulseTag_(conf().crvRecoPulseTag()),
    outputTag_(conf().outputTag()),
    diagLevel_(conf().diagLevel()),
    fillSectorMPV_(conf().fillSectorMPV()),
    dqm_(makeHelperConfig(conf()))
{}

void CRVRecoDQMAnalyzer::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  if (outputTag_.empty()) {
    // TFileService already gives each module its own directory; book into it.
    art::TFileDirectory dir = *tfs;
    dqm_.Book(dir);
  } else {
    dqm_.Book(tfs->mkdir(outputTag_));
  }
}

void CRVRecoDQMAnalyzer::beginRun(const art::Run&)
{
  if (!fillSectorMPV_ || sectorMapSent_) {
    return;
  }
  sectorMapSent_ = true;

  GeomHandle<CosmicRayShield> CRS;

  //cluster position axes from the CRV envelope, as CrvDQMcollector does. Gated
  //on fillSectorMPV because that is what says this job configured a geometry;
  //without it the helper falls back to the Config ranges on first fill.
  // Built from the counters, not CosmicRayShield::getSectorHalfLengths: that one
  // measures the aluminum sheets, and a countersOnly module (all of the extracted
  // geometry) has none, so it throws. getHalfLengths() is already indexed by world
  // axis -- see CRSScintillatorBarDetail::getHalfThickness(), which reads
  // _halfLengths[_localToWorld[0]] -- so no permutation is needed.
  std::vector<double> crvMin(3, 0.0), crvMax(3, 0.0);
  bool firstBar = true;
  for (const auto& bar : CRS->getAllCRSScintillatorBars()) {
    const std::vector<double>& hl = bar->getHalfLengths();
    if (hl.size() < 3) continue;
    const CLHEP::Hep3Vector& pos = bar->getPosition();
    for (int i = 0; i < 3; ++i) {
      const double lo = pos[i] - hl[i];
      const double hi = pos[i] + hl[i];
      if (firstBar) {
        crvMin[i] = lo;
        crvMax[i] = hi;
      } else {
        crvMin[i] = std::min(crvMin[i], lo);
        crvMax[i] = std::max(crvMax[i], hi);
      }
    }
    firstBar = false;
  }
  dqm_.BookPositionAxes(crvMin, crvMax);  //empty envelope falls back to fcl

  auto const& crvSectors = CRS->getCRSScintillatorShields();
  std::vector<std::string> sectorNames;
  sectorNames.reserve(crvSectors.size());
  for (std::size_t i = 0; i < crvSectors.size(); ++i) {
    sectorNames.emplace_back(crvSectors.at(i).name(""));
  }

  // No Proditions here, so unlike CrvDQMcollector this keeps notConnected
  // channels; they enter the sector MPV hists at zero.
  auto const& crvCounters = CRS->getAllCRSScintillatorBars();
  std::vector<int> channelToSector(crvCounters.size() * CRVId::nChanPerBar, -1);
  for (std::size_t channel = 0; channel < channelToSector.size(); ++channel) {
    CRSScintillatorBarIndex barIndex(
        static_cast<int>(channel / CRVId::nChanPerBar));
    channelToSector[channel] = CRS->getBar(barIndex).id().getShieldNumber();
  }

  dqm_.BookSectorMPV(sectorNames, channelToSector);
}

void CRVRecoDQMAnalyzer::analyze(const art::Event& event)
{
  art::Handle<CrvCoincidenceClusterCollection> clusterHandle;
  event.getByLabel(crvCoincidenceClusterTag_, clusterHandle);
  if (!clusterHandle.isValid() || clusterHandle.product() == nullptr) {
    if (!warnedMissingClusters_) {
      warnedMissingClusters_ = true;
      mf::LogWarning("CRVRecoDQMAnalyzer")
          << "No CrvCoincidenceClusterCollection at "
          << crvCoincidenceClusterTag_ << ". Event skipped. "
          << "Reported once per job.";
    }
    return;
  }

  if (crvRecoPulseTag_.empty()) {
    dqm_.Fill(*clusterHandle);
    return;
  }

  art::Handle<CrvRecoPulseCollection> pulseHandle;
  event.getByLabel(crvRecoPulseTag_, pulseHandle);
  if (!pulseHandle.isValid() || pulseHandle.product() == nullptr) {
    if (!warnedMissingPulses_) {
      warnedMissingPulses_ = true;
      mf::LogWarning("CRVRecoDQMAnalyzer")
          << "No CrvRecoPulseCollection at " << crvRecoPulseTag_
          << ". The inclusive reco-pulse plots stay empty; everything driven "
          << "by the clusters still fills. Reported once per job.";
    }
    dqm_.Fill(*clusterHandle);
    return;
  }

  dqm_.Fill(*clusterHandle, *pulseHandle);
}

void CRVRecoDQMAnalyzer::endJob()
{
  dqm_.WriteGraphs();

  if (diagLevel_ > 0) {
    std::cout << "[CRVRecoDQMAnalyzer] Total events: " << dqm_.nEvents()
              << std::endl;
    std::cout << "[CRVRecoDQMAnalyzer] Events with coincidence clusters: "
              << dqm_.nEventsWithClusters() << std::endl;
    std::cout << "[CRVRecoDQMAnalyzer] Coincidence clusters: "
              << dqm_.nClusters() << std::endl;
    std::cout << "[CRVRecoDQMAnalyzer] Reco pulses in clusters: "
              << dqm_.nRecoPulses() << " of " << dqm_.nInclusiveRecoPulses()
              << " reconstructed" << std::endl;
    std::cout << "[CRVRecoDQMAnalyzer] Channel fits: " << dqm_.nFitsSucceeded()
              << " of " << dqm_.nFits() << " attempted, mean chi2/ndf "
              << dqm_.meanFitChi2PerNdf() << std::endl;
    if (dqm_.nOnlineIdOutOfRange() > 0) {
      std::cout << "[CRVRecoDQMAnalyzer] CRVId-range skips (crvPEsMPV_ROC*): "
                << dqm_.nOnlineIdOutOfRange() << std::endl;
    }
    if (dqm_.nOfflineChannelOutOfRange() > 0) {
      std::cout << "[CRVRecoDQMAnalyzer] Offline-channel-range skips: "
                << dqm_.nOfflineChannelOutOfRange() << std::endl;
    }
    if (dqm_.nNullRecoPulsePtrs() > 0) {
      std::cout << "[CRVRecoDQMAnalyzer] Null CrvRecoPulse Ptrs: "
                << dqm_.nNullRecoPulsePtrs() << std::endl;
    }
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CRVRecoDQMAnalyzer)
