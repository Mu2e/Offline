//
// A module to find clusters of coincidences of CRV pulses
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
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
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"

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
#include <TGraph.h>

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
      //fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("label of CrvReco module")};
      fhicl::Atom<std::string> crvCoincidenceClusterFinderModuleLabel{Name("crvCoincidenceClusterFinderModuleLabel"), Comment("label of CoincidenceClusterFinder module")};
      fhicl::Atom<std::string> crvDaqErrorModuleLabel{Name("crvDaqErrorModuleLabel"), Comment("label of module that found the CRV-DAQ errors")};

      fhicl::Atom<int>    histPEsBins{Name("histPEsBins"), Comment("number of bins for PE histograms"), 150};
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
      fhicl::Atom<int>    histDigisNZSBins{Name("histDigisNZSBins"), Comment("number of bins for NZS digis per channel and event histograms"), 200};
      fhicl::Atom<double> histDigisNZSStart{Name("histDigisNZSStart"), Comment("range start for NZS digis per channel and event histograms"), 0};
      fhicl::Atom<double> histDigisNZSEnd{Name("histDigisNZSEnd"), Comment("range end for NZS digis per channel and event histograms"), 0.1};
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
    //std::string _crvRecoPulsesModuleLabel;
    std::string _crvCoincidenceClusterFinderModuleLabel;
    std::string _crvDaqErrorModuleLabel;

    int                _histPEsBins;
    double             _histPEsStart;
    double             _histPEsEnd;
    int                _histPedestalsBins;
    double             _histPedestalsStart;
    double             _histPedestalsEnd;
    int                _histCalibConstsBins;
    double             _histCalibConstsStart;
    double             _histCalibConstsEnd;
    int                _histDigisBins;
    double             _histDigisStart;
    double             _histDigisEnd;
    int                _histDigisNZSBins;
    double             _histDigisNZSStart;
    double             _histDigisNZSEnd;

    int                _totalEvents;
    int                _totalEventsWithCoincidenceClusters;
    int                _totalEventsWithDAQerrors;
    std::pair<int,int> _firstRunSubrun;
    std::pair<int,int> _lastRunSubrun;

    std::vector<int>   _nCoincidences;       //for each sector
    std::vector<int>   _nDigis, _nDigisNZS;  //for each channel
    std::vector<int>   _nDigisROC, _nDigisROCNZS;  //for each channel
    std::vector<TH1F*> _histPEs;             //for each channel
    std::vector<TH1F*> _histPEsROC;          //for each channel
    std::vector<bool>  _notConnected;        //for each channel

    std::vector<TH1F*> _histDigisPerChannelAndEvent;
    std::vector<TH1F*> _histDigisPerChannelAndEventNZS;
    std::vector<TH1F*> _histDigiRatesROC;
    std::vector<TH1F*> _histDigiRatesROCNZS;
    std::vector<TH1F*> _histPEsMPV;
    std::vector<TH1F*> _histPEsMPVROC;
    std::vector<TH1F*> _histPedestals;
    std::vector<TH1F*> _histCalibConstants;
    TH1I*              _histCoincidenceClusters;
    TTree*             _treeMetaData;

    ProditionsHandle<CRVCalib>  _calib;
    ProditionsHandle<CRVStatus> _sipmStatus;
    ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;
  };

  CrvDQMcollector::CrvDQMcollector(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _useDQMcollector(conf().useDQMcollector()),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _crvDigiModuleLabelNZS(conf().crvDigiModuleLabelNZS()),
    //_crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _crvCoincidenceClusterFinderModuleLabel(conf().crvCoincidenceClusterFinderModuleLabel()),
    _crvDaqErrorModuleLabel(conf().crvDaqErrorModuleLabel()),
    _histPEsBins(conf().histPEsBins()),
    _histPEsStart(conf().histPEsStart()),
    _histPEsEnd(conf().histPEsEnd()),
    _histPedestalsBins(conf().histPedestalsBins()),
    _histPedestalsStart(conf().histPedestalsStart()),
    _histPedestalsEnd(conf().histPedestalsEnd()),
    _histCalibConstsBins(conf().histCalibConstsBins()),
    _histCalibConstsStart(conf().histCalibConstsStart()),
    _histCalibConstsEnd(conf().histCalibConstsEnd()),
    _histDigisBins(conf().histDigisBins()),
    _histDigisStart(conf().histDigisStart()),
    _histDigisEnd(conf().histDigisEnd()),
    _histDigisNZSBins(conf().histDigisNZSBins()),
    _histDigisNZSStart(conf().histDigisNZSStart()),
    _histDigisNZSEnd(conf().histDigisNZSEnd()),
    _totalEvents(0),
    _totalEventsWithCoincidenceClusters(0),
    _totalEventsWithDAQerrors(0)
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

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t ROC=1; ROC<=CRVId::nROC; ++ROC)
    {
      for(size_t channel=0; channel<CRVId::nFEBPerROC*CRVId::nChanPerFEB; ++channel)
      {
        _histDigiRatesROC.at(ROC-1)->Fill(channel,(float)(_nDigisROC.at((ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB+channel))/_totalEvents);
        _histDigiRatesROCNZS.at(ROC-1)->Fill(channel,(float)(_nDigisROCNZS.at((ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB+channel))/_totalEvents);

        float PEsMPV = _histPEsROC.at((ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB+channel)->GetMean();  //FIXME: replace by MPV from Gauss-Landau fit
        _histPEsMPVROC.at(ROC-1)->Fill(channel,PEsMPV);
      }
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
    _histPEsMPVROC.reserve(CRVId::nROC);
    _histPedestals.reserve(crvSectors.size());
    _histCalibConstants.reserve(crvSectors.size());
    _histDigisPerChannelAndEvent.reserve(crvSectors.size());
    _histDigisPerChannelAndEventNZS.reserve(crvSectors.size());
    _histDigiRatesROC.reserve(CRVId::nROC);
    _histDigiRatesROCNZS.reserve(CRVId::nROC);
    _nCoincidences.resize(crvSectors.size());
    _nDigis.resize(crvCounters.size()*CRVId::nChanPerBar);
    _nDigisNZS.resize(crvCounters.size()*CRVId::nChanPerBar);
    _nDigisROC.resize(CRVId::nROC*CRVId::nFEBPerROC*CRVId::nChanPerFEB);
    _nDigisROCNZS.resize(CRVId::nROC*CRVId::nFEBPerROC*CRVId::nChanPerFEB);
    _histPEs.reserve(crvCounters.size()*CRVId::nChanPerBar);
    _histPEsROC.reserve(CRVId::nROC*CRVId::nFEBPerROC*CRVId::nChanPerFEB);
    _notConnected.resize(crvCounters.size()*CRVId::nChanPerBar);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t i=0; i<crvCounters.size()*CRVId::nChanPerBar; ++i)
    {
      _histPEs.emplace_back(new TH1F(Form("crvPEs_channel%lu",i), Form("crvPEs_channel%lu",i), _histPEsBins,_histPEsStart,_histPEsEnd));
    }
    for(size_t i=0; i<CRVId::nROC*CRVId::nFEBPerROC*CRVId::nChanPerFEB; ++i)
    {
      _histPEsROC.emplace_back(new TH1F(Form("crvPEsROC_channel%lu",i), Form("crvPEsROC_channel%lu",i), _histPEsBins,_histPEsStart,_histPEsEnd));
    }
    for(size_t i=0; i<crvSectors.size(); ++i)
    {
      _histPEsMPV.emplace_back(tfs->make<TH1F>(Form("crvPEsMPV_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvPEsMPV_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histPEsBins,_histPEsStart,_histPEsEnd));
      _histPedestals.emplace_back(tfs->make<TH1F>(Form("crvPedestals_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvPedestals_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histPedestalsBins,_histPedestalsStart,_histPedestalsEnd));
      _histCalibConstants.emplace_back(tfs->make<TH1F>(Form("crvCalibConstants_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvCalibConstants_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histCalibConstsBins,_histCalibConstsStart,_histCalibConstsEnd));
      _histDigisPerChannelAndEvent.emplace_back(tfs->make<TH1F>(Form("crvDigisPerChannelAndEvent_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvDigisPerChannelAndEvent_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histDigisBins,_histDigisStart,_histDigisEnd));
      _histDigisPerChannelAndEventNZS.emplace_back(tfs->make<TH1F>(Form("crvDigisPerChannelAndEventNZS_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            Form("crvDigisPerChannelAndEventNZS_CRVsector%s",crvSectors.at(i).name("").c_str()),
                                            _histDigisNZSBins,_histDigisNZSStart,_histDigisNZSEnd));
    }
    for(size_t ROC=1; ROC<=CRVId::nROC; ++ROC)
    {
      _histPEsMPVROC.emplace_back(tfs->make<TH1F>(Form("crvPEsMPV_ROC%zu",ROC),
                                            Form("crvPEsMPV_ROC%zu",ROC),
                                            CRVId::nFEBPerROC*CRVId::nChanPerFEB,0,CRVId::nFEBPerROC*CRVId::nChanPerFEB-1));
      _histDigiRatesROC.emplace_back(tfs->make<TH1F>(Form("crvDigiRates_ROC%zu",ROC),
                                            Form("crvDigiRates_ROC%zu",ROC),
                                            CRVId::nFEBPerROC*CRVId::nChanPerFEB,0,CRVId::nFEBPerROC*CRVId::nChanPerFEB-1));
      _histDigiRatesROCNZS.emplace_back(tfs->make<TH1F>(Form("crvDigiRatesNZS_ROC%zu",ROC),
                                            Form("crvDigiRatesNZS_ROC%zu",ROC),
                                            CRVId::nFEBPerROC*CRVId::nChanPerFEB,0,CRVId::nFEBPerROC*CRVId::nChanPerFEB-1));
    }
    _histCoincidenceClusters=tfs->make<TH1I>("crvCoincidencesClusters","crvCoincidenceClusters",10,0,10);

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
    art::Handle<CrvDigiCollection> crvDigiCollectionNZS;
    //art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;
    art::Handle<CrvDAQerrorCollection> crvDaqErrorCollection;

    event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection);
    event.getByLabel(_crvDigiModuleLabelNZS,"NZS",crvDigiCollectionNZS);
    //event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);
    event.getByLabel(_crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusterCollection);
    event.getByLabel(_crvDaqErrorModuleLabel,"",crvDaqErrorCollection);

    auto const& calib = _calib.get(event.id());
    auto const& sipmStatus = _sipmStatus.get(event.id());
    auto const& channelMap = _channelMap_h.get(event.id());

    for(size_t i=0; i<crvDigiCollection->size(); ++i)
    {
      const CrvDigi &digi = crvDigiCollection->at(i);
      int barIndex = digi.GetScintillatorBarIndex().asUint();
      int SiPM = digi.GetSiPMNumber();
      size_t channel = barIndex*CRVId::nChanPerBar + SiPM;
      ++_nDigis.at(channel);

      int ROC=digi.GetROC();
      int ROCport=digi.GetFEB();
      int FEBchannel=digi.GetFEBchannel();
      size_t channelOnline = (ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB + (ROCport-1)*CRVId::nChanPerFEB + FEBchannel;
      ++_nDigisROC.at(channelOnline);
    }

    for(size_t i=0; i<crvDigiCollectionNZS->size(); ++i)
    {
      const CrvDigi &digi = crvDigiCollectionNZS->at(i);
      int barIndex = digi.GetScintillatorBarIndex().asUint();
      int SiPM = digi.GetSiPMNumber();
      size_t channel = barIndex*CRVId::nChanPerBar + SiPM;
      ++_nDigisNZS.at(channel);

      int ROC=digi.GetROC();
      int ROCport=digi.GetFEB();
      int FEBchannel=digi.GetFEBchannel();
      size_t channelOnline = (ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB + (ROCport-1)*CRVId::nChanPerFEB + FEBchannel;
      ++_nDigisROCNZS.at(channelOnline);
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

      const std::vector<art::Ptr<CrvRecoPulse> > &recoPulses = crvCoincidenceClusterCollection->at(i).GetCrvRecoPulses();
      for(const auto recoPulse: recoPulses)
      {
        int barIndex = recoPulse->GetScintillatorBarIndex().asUint();
        int SiPM = recoPulse->GetSiPMNumber();
        size_t channel = barIndex*CRVId::nChanPerBar + SiPM;
        float PEs =recoPulse->GetPEs();
        _histPEs.at(channel)->Fill(PEs);

        mu2e::CRVROC onlineChannel = channelMap.online(channel);
        int ROC=onlineChannel.ROC();
        int ROCport=onlineChannel.FEB();
        int FEBchannel=onlineChannel.FEBchannel();
        size_t onlineChannelIndex = (ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB + (ROCport-1)*CRVId::nChanPerFEB + FEBchannel;
        _histPEsROC.at(onlineChannelIndex)->Fill(PEs);
      }
    }
    if(crvCoincidenceClusterCollection->size()>0) ++_totalEventsWithCoincidenceClusters;

    for(size_t i=0; i<crvDaqErrorCollection->size(); ++i)
    {
      if(crvDaqErrorCollection->at(i).GetErrorCode()!=mu2e::CrvDAQerrorCode::wrongSubsystemID)  //don't count this error
      {
        ++_totalEventsWithDAQerrors;
        break;
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvDQMcollector;
DEFINE_ART_MODULE(CrvDQMcollector)
