//
// A module that finds the pedestals for Crv channels
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include <TH1F.h>

namespace mu2e
{
  class CrvPedestalFinder : public art::EDAnalyzer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvDigiModuleLabel{Name("crvDigiModuleLabel"), Comment("module label for CrvDigis")};
      fhicl::Atom<double>      maxADCspread{Name("maxADCspread"), Comment("maximum spread of ADC values within a waveform to be considered for the pedestal")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvPedestalFinder(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);

    private:
    std::string        _crvDigiModuleLabel;
    double             _maxADCspread;
    std::vector<TH1F*> _pedestalHists;
  };


  CrvPedestalFinder::CrvPedestalFinder(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()), _maxADCspread(conf().maxADCspread())
  {
  }

  void CrvPedestalFinder::beginRun(const art::Run&)
  {
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    _pedestalHists.reserve(counters.size()*4);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t barIndex=0; barIndex<counters.size(); ++barIndex)
    {
      for(size_t SiPM=0; SiPM<4; ++SiPM)
      {
        //produce histograms also for non-existing channels to get a continuously running index
        size_t channelIndex=barIndex*4 + SiPM;
        _pedestalHists.emplace_back(tfs->make<TH1F>(Form("crvPedestalHist_%lu",channelIndex),
                                                    Form("crvPedestalHist_%lu",channelIndex),
                                                    200,-50,150));   //TODO: needs to be only between -50 and +50, but Offline currently sets the pedestal at +100
      }
    }
  }

  void CrvPedestalFinder::analyze(const art::Event& event)
  {
    art::Handle<CrvDigiCollection> crvDigiCollection;
    if(!event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection)) return;

    art::ServiceHandle<art::TFileService> tfs;
    for(auto iter=crvDigiCollection->begin(); iter!=crvDigiCollection->end(); ++iter)
    {
      auto minmaxTest = std::minmax_element(iter->GetADCs().begin(),iter->GetADCs().end());
      if(*minmaxTest.second-*minmaxTest.first>_maxADCspread) continue;  //ignore waveforms that fluctuate too much.
                                                                        //they probably contain a signal or noise pulse

      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      int channelIndex=barIndex*4+SiPM;
      auto hist = _pedestalHists.at(channelIndex);

      for(size_t i=0; i<CrvDigi::NSamples; ++i)
      {
        hist->Fill(iter->GetADCs()[i]);
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvPedestalFinder;
DEFINE_ART_MODULE(CrvPedestalFinder)
