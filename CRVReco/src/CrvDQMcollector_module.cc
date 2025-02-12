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
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

#include <TH1F.h>
#include <TH1I.h>
#include <TTree.h>

#include <string>
#include <array>

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
      fhicl::Atom<std::string> crvDigiModuleLabelNZS{Name("crvDigiModuleLabelNZS"), Comment("label of CrvDigi NZS module")};
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("label of CrvReco module")};
      fhicl::Atom<std::string> crvCoincidenceClusterFinderModuleLabel{Name("crvCoincidenceClusterFinderModuleLabel"),
                                                                      Comment("label of CoincidenceClusterFinder module")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvDQMcollector(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run &run);
    void endJob();

    private:
    bool        _useDQMcollector;
    std::string _crvDigiModuleLabel;
    std::string _crvDigiModuleLabelNZS;
    std::string _crvRecoPulsesModuleLabel;
    std::string _crvCoincidenceClusterFinderModuleLabel;

    int                _totalEvents;
    int                _totalEventsWithCoincidenceClusters;
    std::pair<int,int> _firstRunSubrun;
    std::pair<int,int> _lastRunSubrun;

    std::vector<int>   _nCoincidences;       //for each sector
    std::vector<int>   _nDigis, _nDigisNZS;  //for each channel
    std::vector<TH1F*> _histPEs;             //for each channel
    std::vector<bool>  _notConnected;        //for each channel

    std::vector<TH1F*> _histDigisPerChannelAndEvent;
    std::vector<TH1F*> _histDigisPerChannelAndEventNZS;
    std::vector<TH1F*> _histPEsMPV;
    std::vector<TH1F*> _histPedestals;
    std::vector<TH1F*> _histCalibConstants;
    TH1I*              _histCoincidenceClusters;
    TTree*             _treeMetaData;

    ProditionsHandle<CRVCalib>  _calib;
    ProditionsHandle<CRVStatus> _sipmStatus;
  };

  CrvDQMcollector::CrvDQMcollector(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _useDQMcollector(conf().useDQMcollector()),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _crvDigiModuleLabelNZS(conf().crvDigiModuleLabelNZS()),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _crvCoincidenceClusterFinderModuleLabel(conf().crvCoincidenceClusterFinderModuleLabel()),
    _totalEvents(0),
    _totalEventsWithCoincidenceClusters(0)
  {
  }

  void CrvDQMcollector::endJob()
  {
    GeomHandle<CosmicRayShield> CRS;
    auto &crvCounters = CRS->getAllCRSScintillatorBars();
    for(size_t channel=0; channel<crvCounters.size()*CRVId::nChanPerBar; ++channel)
    {
      if(_notConnected.at(channel)) continue;

      CRSScintillatorBarIndex barIndex(channel/CRVId::nChanPerBar);
      int sectorNumber = CRS->getBar(barIndex).id().getShieldNumber();

      _histDigisPerChannelAndEvent.at(sectorNumber)->Fill((float)(_nDigis.at(channel))/_totalEvents);
      _histDigisPerChannelAndEventNZS.at(sectorNumber)->Fill((float)(_nDigisNZS.at(channel))/_totalEvents);

      float PEsMPV = _histPEs.at(channel)->GetMean();  //FIXME: replace by MPV from Gauss-Landau fit
      _histPEsMPV.at(sectorNumber)->Fill(PEsMPV);
    }

    _treeMetaData->Fill();
  }

  void CrvDQMcollector::beginRun(const art::Run &run)
  {
    if(_histPEsMPV.size()>0) return;  //don't initialize again for additional runs

    GeomHandle<CosmicRayShield> CRS;
    auto &crvSectors = CRS->getCRSScintillatorShields();
    auto &crvCounters = CRS->getAllCRSScintillatorBars();
    _histPEsMPV.reserve(crvSectors.size());
    _histPedestals.reserve(crvSectors.size());
    _histCalibConstants.reserve(crvSectors.size());
    _histDigisPerChannelAndEvent.reserve(crvSectors.size());
    _histDigisPerChannelAndEventNZS.reserve(crvSectors.size());
    _nCoincidences.resize(crvSectors.size());
    _nDigis.resize(crvCounters.size()*CRVId::nChanPerBar);
    _nDigisNZS.resize(crvCounters.size()*CRVId::nChanPerBar);
    _histPEs.reserve(crvCounters.size()*CRVId::nChanPerBar);
    _notConnected.resize(crvCounters.size()*CRVId::nChanPerBar);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t i=0; i<crvCounters.size()*CRVId::nChanPerBar; ++i)
    {
      _histPEs.emplace_back(new TH1F(Form("crvPEsMPV_channel%lu",i), Form("crvPEsMPV_channel%lu",i), 150,0,150));
    }
    for(size_t i=0; i<crvSectors.size(); ++i)
    {
      _histPEsMPV.emplace_back(tfs->make<TH1F>(Form("crvPEsMPV_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            Form("crvPEsMPV_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            150,0,150));
      _histPedestals.emplace_back(tfs->make<TH1F>(Form("crvPedestals_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            Form("crvPedestals_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            201,-50.5,150.5));
      _histCalibConstants.emplace_back(tfs->make<TH1F>(Form("crvCalibConstants_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            Form("crvCalibConstants_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            150,0,3000));
      _histDigisPerChannelAndEvent.emplace_back(tfs->make<TH1F>(Form("crvDigisPerChannelAndEvent_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            Form("crvDigisPerChannelAndEvent_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            200,0.0,0.2));
      _histDigisPerChannelAndEventNZS.emplace_back(tfs->make<TH1F>(Form("crvDigisPerChannelAndEventNZS_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            Form("crvDigisPerChannelAndEventNZS_sector%s",crvSectors.at(i).name("CRV_").c_str()),
                                            200,0.0,0.2));
    }
    _histCoincidenceClusters=tfs->make<TH1I>("crvCoincidencesClusters","crvCoincidenceClusters",10,0,10);

    _treeMetaData=tfs->make<TTree>("crvMetaData","crvMetaData");
    _treeMetaData->Branch("runNumberStart",&_firstRunSubrun.first);
    _treeMetaData->Branch("runNumberEnd",&_lastRunSubrun.first);
    _treeMetaData->Branch("subrunNumberStart",&_firstRunSubrun.second);
    _treeMetaData->Branch("subrunNumberEnd",&_lastRunSubrun.second);
    _treeMetaData->Branch("nEvents",&_totalEvents);
    _treeMetaData->Branch("nEventsWithCoincidenceClusters",&_totalEventsWithCoincidenceClusters);
  }

  void CrvDQMcollector::analyze(const art::Event& event)
  {
    ++_totalEvents;

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvDigiCollection> crvDigiCollection;
    art::Handle<CrvDigiCollection> crvDigiCollectionNZS;
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;

    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);
    event.getByLabel(_crvDigiModuleLabelNZS,"NZS",crvDigiCollectionNZS);
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);
    event.getByLabel(_crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusterCollection);

    auto const& calib = _calib.get(event.id());
    auto const& sipmStatus = _sipmStatus.get(event.id());

    for(size_t i=0; i<crvDigiCollection->size(); ++i)
    {
      const CrvDigi &digi = crvDigiCollection->at(i);
      int barIndex = digi.GetScintillatorBarIndex().asUint();
      int SiPM = digi.GetSiPMNumber();
      size_t channel = barIndex*CRVId::nChanPerBar + SiPM;
      ++_nDigis.at(channel);
    }

    for(size_t i=0; i<crvDigiCollectionNZS->size(); ++i)
    {
      const CrvDigi &digi = crvDigiCollectionNZS->at(i);
      int barIndex = digi.GetScintillatorBarIndex().asUint();
      int SiPM = digi.GetSiPMNumber();
      size_t channel = barIndex*CRVId::nChanPerBar + SiPM;
      ++_nDigisNZS.at(channel);
    }

    for(size_t i=0; i<crvRecoPulseCollection->size(); ++i)
    {
      const CrvRecoPulse &recoPulse = crvRecoPulseCollection->at(i);
      int barIndex = recoPulse.GetScintillatorBarIndex().asUint();
      int SiPM = recoPulse.GetSiPMNumber();
      size_t channel = barIndex*CRVId::nChanPerBar + SiPM;
      float PEs =recoPulse.GetPEs();

      if(crvCoincidenceClusterCollection->size()==0) continue;  //TODO: should we do it like this to remove most of the dark counts?

      _histPEs.at(channel)->Fill(PEs);
    }

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
      _firstRunSubrun=std::pair<int,int>(event.run(),event.subRun());
    }
    _lastRunSubrun=std::pair<int,int>(event.run(),event.subRun());

    for(size_t i=0; i<crvCoincidenceClusterCollection->size(); ++i)
    {
      int sectorType = crvCoincidenceClusterCollection->at(i).GetCrvSectorType();
      _histCoincidenceClusters->Fill(sectorType);
    }
    if(crvCoincidenceClusterCollection->size()>0) ++_totalEventsWithCoincidenceClusters;
  }

} // end namespace mu2e

using mu2e::CrvDQMcollector;
DEFINE_ART_MODULE(CrvDQMcollector)
