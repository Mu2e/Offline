//
// A module to find clusters of coincidences of CRV pulses
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVReco/inc/CrvHelper.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DQMHelpers/inc/CRVDigiDQM.hh"
#include "Offline/DQMHelpers/inc/CRVRecoDQM.hh"
#include "Offline/DQMHelpers/inc/CRVStatusDQM.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <TH1F.h>
#include <TTree.h>

#include <algorithm>
#include <bitset>
#include <string>
#include <vector>

namespace mu2e
{
  class CrvDQMcollector : public art::EDAnalyzer
  {
    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<bool> useDQMcollector{Name("useDQMcollector"), Comment("fill DQM values, histograms, ...")};
      fhicl::Atom<std::string> crvDigiModuleLabel{Name("crvDigiModuleLabel"), Comment("label of CrvDigi module")};
      fhicl::Atom<std::string> crvStatusModuleLabel{Name("crvStatusModuleLabel"), Comment("label of CrvStatus module; empty uses crvDigiModuleLabel"), ""};
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("label of CrvReco module; empty skips the inclusive reco-pulse histograms"), ""};
      fhicl::Atom<std::string> crvCoincidenceClusterFinderModuleLabel{Name("crvCoincidenceClusterFinderModuleLabel"), Comment("label of CoincidenceClusterFinder module")};
      fhicl::Atom<std::string> crvDaqErrorModuleLabel{Name("crvDaqErrorModuleLabel"), Comment("label of module that found the CRV-DAQ errors")};
      fhicl::Atom<std::string> crvDigiDQMDir{Name("crvDigiDQMDir"), Comment("TFileService subdirectory for CRVDigiDQM histograms"), "CRVDigiDQM"};
      fhicl::Atom<std::string> crvRecoDQMDir{Name("crvRecoDQMDir"), Comment("TFileService subdirectory for CRVRecoDQM histograms; empty books in the module directory"), ""};
      fhicl::Atom<std::string> crvStatusDQMDir{Name("crvStatusDQMDir"), Comment("TFileService subdirectory for CRVStatusDQM histograms"), "CRVStatusDQM"};
      fhicl::Atom<bool> fillInclusiveDigiDQM{Name("fillInclusiveDigiDQM"), Comment("also fill BarId/SiPM/ADC in CRVDigiDQM"), true};
      fhicl::Atom<bool> fillInclusiveRecoDQM{Name("fillInclusiveRecoDQM"), Comment("also fill the per-event DqmCrv reco-pulse and cluster plots in CRVRecoDQM"), true};
      fhicl::Atom<bool> crvDigiDQMkppReadout{Name("crvDigiDQMkppReadout"), Comment("KPP FEB-axis sizing (ROC 1-2); ROC4->ROC2 is the unpacker's job"), true};
      fhicl::Atom<bool> writePerChannelPE{Name("writePerChannelPE"), Comment("write the ~50k per-channel PE spectra (expert output)"), false};

      fhicl::Atom<int>    histPEsBins{Name("histPEsBins"), Comment("number of bins for PE histograms"), 75};
      fhicl::Atom<double> histPEsStart{Name("histPEsStart"), Comment("range start for PE histograms"), 0};
      fhicl::Atom<double> histPEsEnd{Name("histPEsEnd"), Comment("range end for PE histograms"), 150};
      fhicl::Atom<int>    histPedestalsBins{Name("histPedestalsBins"), Comment("number of bins for pedestal histograms"), 200};
      fhicl::Atom<double> histPedestalsStart{Name("histPedestalsStart"), Comment("range start for pedestal histograms"), 1950};
      fhicl::Atom<double> histPedestalsEnd{Name("histPedestalsEnd"), Comment("range end for pedestal histograms"), 2150};
      fhicl::Atom<int>    histCalibConstsBins{Name("histCalibConstsBins"), Comment("number of bins for calib consts histograms"), 100};
      fhicl::Atom<double> histCalibConstsStart{Name("histCalibConstsStart"), Comment("range start for calib consts histograms"), 0};
      fhicl::Atom<double> histCalibConstsEnd{Name("histCalibConstsEnd"), Comment("range end for calib consts histograms"), 2000};
      fhicl::Atom<int>    histDigisBins{Name("histDigisBins"), Comment("number of bins for digis per channel and event histograms"), 200};
      fhicl::Atom<double> histDigisStart{Name("histDigisStart"), Comment("range start for digis per channel and event histograms"), 0};
      fhicl::Atom<double> histDigisEnd{Name("histDigisEnd"), Comment("range end for digis per channel and event histograms"), 0.1};
      fhicl::Atom<double> PEfitRangeStart{Name("PEfitRangeStart"), Comment("low end of the PE MPV fit range as fraction of peak"), 0.7};
      fhicl::Atom<double> PEfitRangeEnd{Name("PEfitRangeEnd"), Comment("high end of the PE MPV fit range as fraction of peak"), 2.0};
      fhicl::Atom<double> PEstart{Name("PEstart"), Comment("lowest PE for fit"), 15.0};

      //inclusive-plot axes that depend on detector position and readout window.
      //defaults are the DqmCrv values (full CRV, Mu2e coordinates); the extracted
      //CRV sits at y~4237-4649, z~21440-23065 and needs minY/maxY, minZ/maxZ moved.
      fhicl::Atom<int>    nBinsTime{Name("nBinsTime"), Comment("bins for PulseTime/LeadingTime/tc"), 100};
      fhicl::Atom<double> minTime{Name("minTime"), Comment("range start for PulseTime/LeadingTime/tc [ns]"), 0.0};
      fhicl::Atom<double> maxTime{Name("maxTime"), Comment("range end for PulseTime/LeadingTime/tc [ns]"), 2000.0};
      fhicl::Atom<int>    nBinsTime2{Name("nBinsTime2"), Comment("bins for PulseTime2/LeadingTime2/t2c"), 100};
      fhicl::Atom<double> minTime2{Name("minTime2"), Comment("range start for PulseTime2/LeadingTime2/t2c [ns]"), 0.0};
      fhicl::Atom<double> maxTime2{Name("maxTime2"), Comment("range end for PulseTime2/LeadingTime2/t2c [ns]"), 100000.0};
      fhicl::Atom<int>    nBinsPos{Name("nBinsPos"), Comment("bins for the cluster position plots X/Y/Z"), 100};
      fhicl::Atom<double> minX{Name("minX"), Comment("range start for X [mm]"), -6904.0};
      fhicl::Atom<double> maxX{Name("maxX"), Comment("range end for X [mm]"), -904.0};
      fhicl::Atom<double> minY{Name("minY"), Comment("range start for Y [mm]"), 0.0};
      fhicl::Atom<double> maxY{Name("maxY"), Comment("range end for Y [mm]"), 3000.0};
      fhicl::Atom<double> minZ{Name("minZ"), Comment("range start for Z [mm]"), -3500.0};
      fhicl::Atom<double> maxZ{Name("maxZ"), Comment("range end for Z [mm]"), 20000.0};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvDQMcollector(const Parameters& config);
    void beginJob() override;
    void analyze(const art::Event& e) override;
    void beginRun(const art::Run &run) override;
    void endSubRun(const art::SubRun &subRun) override;
    void endJob() override;

    private:
    bool        _useDQMcollector;
    std::string _crvDigiModuleLabel;
    std::string _crvStatusModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    std::string _crvCoincidenceClusterFinderModuleLabel;
    std::string _crvDaqErrorModuleLabel;
    std::string _crvDigiDQMDir;
    std::string _crvRecoDQMDir;
    std::string _crvStatusDQMDir;
    bool _warnedMissingDigi{false};
    bool _warnedMissingStatus{false};
    bool _warnedMissingClusters{false};
    bool _warnedMissingPulses{false};
    bool _warnedMissingDaqErrors{false};

    int                _histPedestalsBins;
    double             _histPedestalsStart;
    double             _histPedestalsEnd;
    int                _histCalibConstsBins;
    double             _histCalibConstsStart;
    double             _histCalibConstsEnd;
    int                _histDigisBins;
    double             _histDigisStart;
    double             _histDigisEnd;

    int                _totalEvents;
    int                _totalEventsWithCoincidenceClusters;
    int                _totalEventsWithDAQerrors;
    std::pair<int,int> _firstRunSubrun;
    std::pair<int,int> _lastRunSubrun;

    std::vector<bool>  _notConnected;        //for each channel

    std::vector<TH1F*> _histPedestals;
    std::vector<TH1F*> _histCalibConstants;
    TTree*             _treeMetaData;

    ProditionsHandle<CRVCalib>  _calib;
    ProditionsHandle<CRVStatus> _sipmStatus;

    CRVDigiDQM _digiDQM;
    CRVRecoDQM _recoDQM;
    CRVStatusDQM _statusDQM;
  };

  CrvDQMcollector::CrvDQMcollector(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _useDQMcollector(conf().useDQMcollector()),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _crvStatusModuleLabel(conf().crvStatusModuleLabel().empty() ?
                          conf().crvDigiModuleLabel() : conf().crvStatusModuleLabel()),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _crvCoincidenceClusterFinderModuleLabel(conf().crvCoincidenceClusterFinderModuleLabel()),
    _crvDaqErrorModuleLabel(conf().crvDaqErrorModuleLabel()),
    _crvDigiDQMDir(conf().crvDigiDQMDir()),
    _crvRecoDQMDir(conf().crvRecoDQMDir()),
    _crvStatusDQMDir(conf().crvStatusDQMDir()),
    _histPedestalsBins(conf().histPedestalsBins()),
    _histPedestalsStart(conf().histPedestalsStart()),
    _histPedestalsEnd(conf().histPedestalsEnd()),
    _histCalibConstsBins(conf().histCalibConstsBins()),
    _histCalibConstsStart(conf().histCalibConstsStart()),
    _histCalibConstsEnd(conf().histCalibConstsEnd()),
    _histDigisBins(conf().histDigisBins()),
    _histDigisStart(conf().histDigisStart()),
    _histDigisEnd(conf().histDigisEnd()),
    _totalEvents(0),
    _totalEventsWithCoincidenceClusters(0),
    _totalEventsWithDAQerrors(0),
    _treeMetaData(nullptr),
    _digiDQM([] (bool fillInclusive, bool kppReadout) {
      CRVDigiDQM::Config c;
      c.fillInclusive = fillInclusive;
      c.kppReadout = kppReadout;
      c.fillLivePlots = false;
      return c;
    }(conf().fillInclusiveDigiDQM(), conf().crvDigiDQMkppReadout())),
    _recoDQM([] (const Config &conf) {
      CRVRecoDQM::Config c;
      c.nBinsPEs = conf.histPEsBins();
      c.minPEs = conf.histPEsStart();
      c.maxPEs = conf.histPEsEnd();
      c.PEfitRangeStart = conf.PEfitRangeStart();
      c.PEfitRangeEnd = conf.PEfitRangeEnd();
      c.PEstart = conf.PEstart();
      c.writePerChannelPE = conf.writePerChannelPE();
      c.fillInclusive = conf.fillInclusiveRecoDQM();
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
    }(conf())),
    //binning defaults are the online ones; fillLivePlots stays false so the
    //per-job file hadds
    _statusDQM([] {
      CRVStatusDQM::Config c;
      c.fillLivePlots = false;
      return c;
    }())
  {
  }

  void CrvDQMcollector::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    _digiDQM.Book(tfs->mkdir(_crvDigiDQMDir));
    //the reco histograms keep their historical place in the module directory
    if(_crvRecoDQMDir.empty())
    {
      art::TFileDirectory dir = *tfs;
      _recoDQM.Book(dir);
    }
    else
    {
      _recoDQM.Book(tfs->mkdir(_crvRecoDQMDir));
    }
    _statusDQM.Book(tfs->mkdir(_crvStatusDQMDir));
  }

  void CrvDQMcollector::endSubRun(const art::SubRun &subRun)
  {
    _statusDQM.EndSubRun(subRun.run(), subRun.subRun());
  }

  void CrvDQMcollector::endJob()
  {
    _digiDQM.WriteGraphs();
    _recoDQM.WriteGraphs();  //runs the per-channel PE fits and fills the MPV maps
    _statusDQM.WriteGraphs();

    if(_totalEvents<=0) return;

    _totalEventsWithCoincidenceClusters = static_cast<int>(_recoDQM.nEventsWithClusters());
    _treeMetaData->Fill();
  }

  void CrvDQMcollector::beginRun(const art::Run &run)
  {
    GeomHandle<CosmicRayShield> CRS;

    //The cluster position axes follow the detector, so take them from the CRV
    //envelope rather than a fixed range that is wrong on some geometry.
    //Built from the counters, not CosmicRayShield::getSectorHalfLengths: that
    //one measures the aluminum sheets, and a countersOnly module (all of the
    //extracted geometry) has none, so it throws. Counters are also the right
    //basis here, since a cluster position is derived from them.
    //getHalfLengths() is already indexed by world axis -- see
    //CRSScintillatorBarDetail::getHalfThickness(), which reads
    //_halfLengths[_localToWorld[0]] -- so no permutation is needed. Idempotent.
    std::vector<double> crvMin(3, 0.0), crvMax(3, 0.0);
    bool firstBar = true;
    for(const auto &bar : CRS->getAllCRSScintillatorBars())
    {
      const std::vector<double> &hl = bar->getHalfLengths();
      if(hl.size()<3) continue;
      const CLHEP::Hep3Vector &pos = bar->getPosition();
      for(int i=0; i<3; ++i)
      {
        const double lo = pos[i]-hl[i];
        const double hi = pos[i]+hl[i];
        if(firstBar) { crvMin[i]=lo; crvMax[i]=hi; }
        else         { crvMin[i]=std::min(crvMin[i],lo); crvMax[i]=std::max(crvMax[i],hi); }
      }
      firstBar = false;
    }
    _recoDQM.BookPositionAxes(crvMin, crvMax);  //empty envelope falls back to fcl

    if(_histPedestals.size()>0) return;  //don't initialize again for additional runs

    auto &crvSectors = CRS->getCRSScintillatorShields();
    auto &crvCounters = CRS->getAllCRSScintillatorBars();
    _histPedestals.reserve(crvSectors.size());
    _histCalibConstants.reserve(crvSectors.size());
    _notConnected.resize(crvCounters.size()*CRVId::nChanPerBar);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t i=0; i<crvSectors.size(); ++i)
    {
      _histPedestals.emplace_back(tfs->make<TH1F>(Form("crvPedestals_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvPedestals_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histPedestalsBins,_histPedestalsStart,_histPedestalsEnd));
      _histCalibConstants.emplace_back(tfs->make<TH1F>(Form("crvCalibConstants_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvCalibConstants_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histCalibConstsBins,_histCalibConstsStart,_histCalibConstsEnd));
    }

    _treeMetaData=tfs->make<TTree>("crvMetaData","crvMetaData");
    _treeMetaData->Branch("runNumberStart",&_firstRunSubrun.first);
    _treeMetaData->Branch("runNumberEnd",&_lastRunSubrun.first);
    _treeMetaData->Branch("subrunNumberStart",&_firstRunSubrun.second);
    _treeMetaData->Branch("subrunNumberEnd",&_lastRunSubrun.second);
    _treeMetaData->Branch("nEvents",&_totalEvents);
    _treeMetaData->Branch("nEventsWithCoincidenceClusters",&_totalEventsWithCoincidenceClusters);
    _treeMetaData->Branch("nEventsWithDAQerrors",&_totalEventsWithDAQerrors);
  }

  void CrvDQMcollector::analyze(const art::Event& event)
  {
    ++_totalEvents;

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvDigiCollection> crvDigiCollection;
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;
    art::Handle<CrvDAQerrorCollection> crvDaqErrorCollection;

    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);
    art::Handle<CrvStatusCollection> crvStatusCollection;
    event.getByLabel(_crvStatusModuleLabel,"",crvStatusCollection);
    if(!_crvRecoPulsesModuleLabel.empty()) event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);
    event.getByLabel(_crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusterCollection);
    event.getByLabel(_crvDaqErrorModuleLabel,"",crvDaqErrorCollection);

    auto const& calib = _calib.get(event.id());
    auto const& sipmStatus = _sipmStatus.get(event.id());

    const CrvStatusCollection emptyStatus;
    const CrvDigiCollection emptyDigis;
    const CrvCoincidenceClusterCollection emptyClusters;
    const bool haveDigis =
        crvDigiCollection.isValid() && crvDigiCollection.product()!=nullptr;
    const bool haveStatus =
        crvStatusCollection.isValid() && crvStatusCollection.product()!=nullptr;
    const bool haveClusters =
        crvCoincidenceClusterCollection.isValid() &&
        crvCoincidenceClusterCollection.product()!=nullptr;
    const bool havePulses =
        crvRecoPulseCollection.isValid() &&
        crvRecoPulseCollection.product()!=nullptr;
    const bool haveDaqErrors =
        crvDaqErrorCollection.isValid() &&
        crvDaqErrorCollection.product()!=nullptr;
    if (!haveDigis && !_warnedMissingDigi) {
      _warnedMissingDigi = true;
      mf::LogWarning("CrvDQMcollector")
          << "No CrvDigiCollection at " << _crvDigiModuleLabel
          << " (empty collection used). Reported once per job.";
    }
    if (!haveStatus && !_warnedMissingStatus) {
      _warnedMissingStatus = true;
      mf::LogWarning("CrvDQMcollector")
          << "No CrvStatusCollection at " << _crvStatusModuleLabel
          << " (empty collection used). Reported once per job.";
    }
    if (!haveClusters && !_warnedMissingClusters) {
      _warnedMissingClusters = true;
      mf::LogWarning("CrvDQMcollector")
          << "No CrvCoincidenceClusterCollection at "
          << _crvCoincidenceClusterFinderModuleLabel
          << " (empty collection used). Reported once per job.";
    }
    if (!havePulses && !_crvRecoPulsesModuleLabel.empty() && !_warnedMissingPulses) {
      _warnedMissingPulses = true;
      mf::LogWarning("CrvDQMcollector")
          << "No CrvRecoPulseCollection at " << _crvRecoPulsesModuleLabel
          << ". The inclusive reco-pulse histograms stay empty; everything "
          << "driven by the coincidence clusters still fills. "
          << "Reported once per job.";
    }
    if (!haveDaqErrors && !_warnedMissingDaqErrors) {
      _warnedMissingDaqErrors = true;
      mf::LogWarning("CrvDQMcollector")
          << "No CrvDAQerrorCollection at " << _crvDaqErrorModuleLabel
          << " (empty collection used; nEventsWithDAQerrors stays 0). "
          << "Reported once per job.";
    }
    const CrvDigiCollection& digis = haveDigis ? *crvDigiCollection : emptyDigis;
    const CrvStatusCollection& status = haveStatus ? *crvStatusCollection : emptyStatus;
    const CrvDAQerrorCollection emptyDaqErrors;
    const CrvDAQerrorCollection& daqErrors =
        haveDaqErrors ? *crvDaqErrorCollection : emptyDaqErrors;
    _digiDQM.Fill(digis, status);
    _statusDQM.Fill(status, daqErrors);

    static bool first=true;
    if(first)
    {
      first=false;
      GeomHandle<CosmicRayShield> CRS;
      auto &crvCounters = CRS->getAllCRSScintillatorBars();
      for(size_t channel=0; channel<crvCounters.size()*CRVId::nChanPerBar; ++channel)
      {
        std::bitset<16> status(sipmStatus.status(channel));
        if(status.test(CRVStatus::Flags::notConnected))
        {
          _notConnected.at(channel)=true;
          continue;
        }

        double pedestal = calib.pedestal(channel);
        double calibPulseArea = calib.pulseArea(channel);

        CRSScintillatorBarIndex barIndex(channel/CRVId::nChanPerBar);
        int sectorNumber = CRS->getBar(barIndex).id().getShieldNumber();
        _histPedestals.at(sectorNumber)->Fill(pedestal);
        _histCalibConstants.at(sectorNumber)->Fill(calibPulseArea);
      }
      //helpers own crvDigisPerChannelAndEvent_* and crvPEsMPV_CRVsector*;
      //-1 drops notConnected channels from both
      std::vector<std::string> sectorNames;
      auto &crvSectors = CRS->getCRSScintillatorShields();
      sectorNames.reserve(crvSectors.size());
      for(size_t i=0; i<crvSectors.size(); ++i) sectorNames.emplace_back(crvSectors.at(i).name(""));

      std::vector<int> channelToSector(crvCounters.size()*CRVId::nChanPerBar, -1);
      for(size_t channel=0; channel<channelToSector.size(); ++channel)
      {
        if(_notConnected.at(channel)) continue;
        CRSScintillatorBarIndex barIndex(channel/CRVId::nChanPerBar);
        channelToSector.at(channel)=CRS->getBar(barIndex).id().getShieldNumber();
      }
      _digiDQM.BookSectorOccupancy(sectorNames, channelToSector,
                                   _histDigisBins, _histDigisStart, _histDigisEnd);
      _recoDQM.BookSectorMPV(sectorNames, channelToSector);

      _firstRunSubrun=std::pair<int,int>(event.run(),event.subRun());
    }
    _lastRunSubrun=std::pair<int,int>(event.run(),event.subRun());

    const CrvCoincidenceClusterCollection& clusters =
        haveClusters ? *crvCoincidenceClusterCollection : emptyClusters;
    if(havePulses) _recoDQM.Fill(clusters, *crvRecoPulseCollection);
    else           _recoDQM.Fill(clusters);

    if(haveDaqErrors)
    {
      for(size_t i=0; i<crvDaqErrorCollection->size(); ++i)
      {
        if(crvDaqErrorCollection->at(i).GetErrorCode()!=mu2e::CrvDAQerrorCode::wrongSubsystemID)  //don't count this error
        {
          ++_totalEventsWithDAQerrors;
          break;
        }
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvDQMcollector;
DEFINE_ART_MODULE(CrvDQMcollector)
