//
// A module that finds the calibration constants for Crv channels
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include <TMath.h>
#include <TH1F.h>

namespace mu2e
{
  class PrintCrvChannelMap : public art::EDAnalyzer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> offlineChannelMapName{Name("offlineChannelMap"), Comment("map of offline channels"), "offlineChannelMap.txt"};
      fhicl::Atom<std::string> onlineChannelMapName{Name("onlineChannelMap"), Comment("map of online channels"), "onlineChannelMap.txt"};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit PrintCrvChannelMap(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);

    private:
    std::string                              _offlineChannelMapName;
    std::string                              _onlineChannelMapName;
    mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;
  };


  PrintCrvChannelMap::PrintCrvChannelMap(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _offlineChannelMapName(conf().offlineChannelMapName()),
    _onlineChannelMapName(conf().onlineChannelMapName())
  {
  }

  void PrintCrvChannelMap::beginRun(const art::Run&)
  {
  }

  void PrintCrvChannelMap::analyze(const art::Event& event)
  {
    static bool first=true;
    if(first)
    {
      first=false;

      GeomHandle<CosmicRayShield> CRS;
      const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
      auto const& channelMap = _channelMap_h.get(event.id());
      std::set<size_t> onlineChannelSet;

      //offline
      std::ofstream offlineChannelMap;
      offlineChannelMap.open(_offlineChannelMapName);
      offlineChannelMap<<"offlineChannelIndex,barIndex,SiPM,sectorName,moduleNumber,side,ROC,ROCport,FEBchannel"<<std::endl;
      for(size_t barIndex=0; barIndex<counters.size(); ++barIndex)
      {
        const CRSScintillatorBarId &counterId = counters.at(barIndex)->id();
        const std::string &sectorName         = CRS->getCRSScintillatorShield(counterId.getShieldId()).getName();
        const int moduleNumber                = counterId.getModuleNumber();

        for(size_t SiPM=0; SiPM<CRVId::nChanPerBar; ++SiPM)
        {
          int side=SiPM%CRVId::nSidesPerBar;
          if(!counters.at(barIndex)->hasCMB(side)) continue;  //non-existing offline channel

          size_t offlineChannelIndex=barIndex*CRVId::nChanPerBar+SiPM;

          CRVROC onlineChannel = channelMap.online(offlineChannelIndex);
          uint16_t rocID       = onlineChannel.ROC();
          uint16_t rocPort     = onlineChannel.FEB();
          uint16_t febChannel  = onlineChannel.FEBchannel();

          offlineChannelMap<<offlineChannelIndex<<","<<barIndex<<","<<SiPM<<","<<sectorName<<","<<moduleNumber<<","<<(side==0?"-":"+")<<","<<rocID<<","<<rocPort<<","<<febChannel<<std::endl;

          size_t onlineChannelIndex = (rocID-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB + (rocPort-1)*CRVId::nChanPerFEB + febChannel;
          onlineChannelSet.insert(onlineChannelIndex);
        } //SiPMs
      } //scintillator bars
      offlineChannelMap.close();

      //online
      std::ofstream onlineChannelMap;
      onlineChannelMap.open(_onlineChannelMapName);
      onlineChannelMap<<"ROC,ROCport,FEBchannel,offlineChannelIndex,barIndex,SiPM,sectorName,moduleNumber,side"<<std::endl;
      for(auto channelIter=onlineChannelSet.begin(); channelIter!=onlineChannelSet.end(); ++channelIter)
      {
        uint16_t rocID       = (*channelIter)/(CRVId::nFEBPerROC*CRVId::nChanPerFEB)+1;
        uint16_t rocPort     = ((*channelIter)%(CRVId::nFEBPerROC*CRVId::nChanPerFEB))/CRVId::nChanPerFEB+1;
        uint16_t febChannel  = (*channelIter)%(CRVId::nChanPerFEB);
        CRVROC onlineChannel(rocID, rocPort, febChannel);

        size_t offlineChannelIndex = channelMap.offline(onlineChannel);
        size_t barIndex = offlineChannelIndex/CRVId::nChanPerBar;
        int    SiPM     = offlineChannelIndex%CRVId::nChanPerBar;
        int    side     = SiPM%CRVId::nSidesPerBar;

        const CRSScintillatorBarId &counterId = counters.at(barIndex)->id();
        const std::string &sectorName         = CRS->getCRSScintillatorShield(counterId.getShieldId()).getName();
        const int moduleNumber                = counterId.getModuleNumber();

        onlineChannelMap<<rocID<<","<<rocPort<<","<<febChannel<<","<<offlineChannelIndex<<","<<barIndex<<","<<SiPM<<","<<sectorName<<","<<moduleNumber<<","<<(side==0?"-":"+")<<std::endl;
      }
      onlineChannelMap.close();

    } //first event
  } //analyze

} // end namespace mu2e

using mu2e::PrintCrvChannelMap;
DEFINE_ART_MODULE(PrintCrvChannelMap)
