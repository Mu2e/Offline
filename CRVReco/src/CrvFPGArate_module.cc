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
  class CrvFPGArate : public art::EDAnalyzer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvDigiModuleLabel{Name("crvDigiModuleLabel"), Comment("module label for CrvDigis")};
      fhicl::Atom<int> minTDC{Name("minTDC"), Comment("minimum TDC")};
      fhicl::Atom<int> maxTDC{Name("maxTDC"), Comment("maximum TDC")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvFPGArate(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);

    private:
    std::string                              _crvDigiModuleLabel;
    int                                      _minTDC;
    int                                      _maxTDC;
    std::unique_ptr<art::TFileDirectory>     _tfdirFPGAhitRates;
    std::unique_ptr<art::TFileDirectory>     _tfdirFPGAhitMultiplicities;
    std::map<int,TH1F*>                      _FPGArateHists;
    std::map<int,TH1F*>                      _FPGAmultiplicityHists;
    mu2e::ProditionsHandle<mu2e::CRVOrdinal> _channelMap_h;
  };


  CrvFPGArate::CrvFPGArate(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()),
    _minTDC(conf().minTDC()),
    _maxTDC(conf().maxTDC())
  {
    art::ServiceHandle<art::TFileService> tfs;
    _tfdirFPGAhitRates = std::unique_ptr<art::TFileDirectory>(new art::TFileDirectory(tfs->mkdir("FPGAhitRates")));
    _tfdirFPGAhitMultiplicities = std::unique_ptr<art::TFileDirectory>(new art::TFileDirectory(tfs->mkdir("FPGAhitMultiplicities")));
  }

  void CrvFPGArate::beginRun(const art::Run&)
  {
  }

  void CrvFPGArate::analyze(const art::Event& event)
  {
    GeomHandle<CosmicRayShield> CRS;
    art::ServiceHandle<art::TFileService> tfs;

    auto const& channelMap = _channelMap_h.get(event.id());

    art::Handle<CrvDigiCollection> crvDigiCollection;
    if(!event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection)) return;

    //first create an FPGA map with entries for all FPGAs so that the 0-bin can be filled in the histogram for FPGAs without hits
    std::map<int,int> eventFPGAmap;
    std::map<int,std::vector<int>> eventHitMultiplicityMap;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    for(size_t barIndex=0; barIndex<counters.size(); ++barIndex)
    {
      for(size_t SiPM=0; SiPM<4; ++SiPM)
      {
        const auto counter = counters.at(barIndex);
        if(!counter->hasCMB(SiPM%2)) continue;  //non-existing offline channel

        size_t offlineChannel = barIndex*4 + SiPM;

        CRVROC onlineChannel = channelMap.online(offlineChannel);
        uint16_t rocID       = onlineChannel.ROC();
        uint16_t rocPort     = onlineChannel.FEB();
        uint16_t febChannel  = onlineChannel.FEBchannel();
        uint16_t febFPGA     = febChannel/16;

        uint16_t FPGAIndex = 24*4*(rocID-1) + 4*(rocPort-1) + febFPGA;
        eventFPGAmap[FPGAIndex]=0;
        eventHitMultiplicityMap[FPGAIndex].clear();
      }
    }

    //now loop through all digis and update the FPGA map
    //(use this loop to also check for pulses with more than one hits)
    int previousBarIndex = -1;
    int previousSiPM     = -1;
    int previousTDC      = -1;
    for(auto iter=crvDigiCollection->begin(); iter!=crvDigiCollection->end(); ++iter)
    {
      uint16_t TDC  = iter->GetStartTDC();
      if(TDC<_minTDC || TDC>_maxTDC) continue;

      int barIndex            = iter->GetScintillatorBarIndex().asInt();
      int SiPM                = iter->GetSiPMNumber();
      uint16_t offlineChannel = barIndex*4 + SiPM;

      CRVROC onlineChannel = channelMap.online(offlineChannel);
      uint16_t rocID       = onlineChannel.ROC();
      uint16_t rocPort     = onlineChannel.FEB();
      uint16_t febChannel  = onlineChannel.FEBchannel();
      uint16_t febFPGA     = febChannel/16;

      uint16_t FPGAIndex = 24*4*(rocID-1) + 4*(rocPort-1) + febFPGA;
      eventFPGAmap[FPGAIndex]++;

      //check for pulses with multiple hits
      if(previousBarIndex==barIndex && previousSiPM==SiPM && previousTDC+iter->GetADCs().size()==TDC)
        eventHitMultiplicityMap[FPGAIndex].back()++;
      else
        eventHitMultiplicityMap[FPGAIndex].push_back(1);

      previousBarIndex = barIndex;
      previousSiPM     = SiPM;
      previousTDC      = TDC;
    }

    //now collect the hit rates of all FPGAs
    for(auto iter=eventFPGAmap.begin(); iter!=eventFPGAmap.end(); ++iter)
    {
      uint16_t FPGAIndex = iter->first;
      auto FPGAiter = _FPGArateHists.find(FPGAIndex);
      if(FPGAiter ==_FPGArateHists.end())
      {
        uint16_t rocID      = FPGAIndex/(24*4)+1;
        uint16_t rocPort    = (FPGAIndex%(24*4))/4+1;
        uint16_t febFPGA    = FPGAIndex%4;
        uint16_t febChannel = febFPGA*16;  //use the 1st channel of the FPGA
        mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

        uint16_t offlineChannel = channelMap.offline(onlineChannel);
        int barIndex            = offlineChannel/4;
        const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
        const CRSScintillatorBarId &counterId = counters.at(barIndex)->id();
        const std::string &sectorName = CRS->getCRSScintillatorShield(counterId.getShieldId()).getName();
        const int moduleNumber = counterId.getModuleNumber();

        const std::string histName=Form("%i_%i_%i__%s_%i",rocID,rocPort,febFPGA,sectorName.c_str(),moduleNumber);

        FPGAiter = _FPGArateHists.emplace(FPGAIndex,_tfdirFPGAhitRates->make<TH1F>(histName.c_str(),histName.c_str(),100,0,100)).first;  //returns iterator to new element
      }

      FPGAiter->second->Fill(iter->second);
    }

    //now collect the hit mulitplicities of all FPGAs
    for(auto iter=eventHitMultiplicityMap.begin(); iter!=eventHitMultiplicityMap.end(); ++iter)
    {
      uint16_t FPGAIndex = iter->first;
      auto FPGAiter = _FPGAmultiplicityHists.find(FPGAIndex);
      if(FPGAiter ==_FPGAmultiplicityHists.end())
      {
        uint16_t rocID      = FPGAIndex/(24*4)+1;
        uint16_t rocPort    = (FPGAIndex%(24*4))/4+1;
        uint16_t febFPGA    = FPGAIndex%4;
        uint16_t febChannel = febFPGA*16;  //use the 1st channel of the FPGA
        mu2e::CRVROC onlineChannel(rocID, rocPort, febChannel);

        uint16_t offlineChannel = channelMap.offline(onlineChannel);
        int barIndex            = offlineChannel/4;
        const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
        const CRSScintillatorBarId &counterId = counters.at(barIndex)->id();
        const std::string &sectorName = CRS->getCRSScintillatorShield(counterId.getShieldId()).getName();
        const int moduleNumber = counterId.getModuleNumber();

        const std::string histName=Form("%i_%i_%i__%s_%i",rocID,rocPort,febFPGA,sectorName.c_str(),moduleNumber);

        FPGAiter = _FPGAmultiplicityHists.emplace(FPGAIndex,_tfdirFPGAhitMultiplicities->make<TH1F>(histName.c_str(),histName.c_str(),100,0,100)).first;  //returns iterator to new element
      }

      std::vector<int> &multiplicities = iter->second;
      for(size_t i=0; i<multiplicities.size(); ++i) FPGAiter->second->Fill(multiplicities.at(i));
    }
  }

} // end namespace mu2e

using mu2e::CrvFPGArate;
DEFINE_ART_MODULE(CrvFPGArate)
