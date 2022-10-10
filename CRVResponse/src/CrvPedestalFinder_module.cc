//
// A module that finds the pedestals for Crv channels
//
//
// Original Author: Ralf Ehrlich

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
    void endJob();

    private:
    std::string _crvDigiModuleLabel;
    double      _maxADCspread;

    std::map<std::pair<int,int>,TH1F*> _pedestalHists;
  };


  CrvPedestalFinder::CrvPedestalFinder(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()), _maxADCspread(conf().maxADCspread())
  {
  }

  void CrvPedestalFinder::endJob()
  {
  }

  void CrvPedestalFinder::analyze(const art::Event& event)
  {
    art::Handle<CrvDigiCollection> crvDigiCollection;
    if(!event.getByLabel(_crvDigiModuleLabel,"",crvDigiCollection)) return;

    art::ServiceHandle<art::TFileService> tfs;
    for(auto iter=crvDigiCollection->begin(); iter!=crvDigiCollection->end(); ++iter)
    {
      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      auto histIter = _pedestalHists.find(std::make_pair(barIndex,SiPM));
      if(histIter==_pedestalHists.end())
      {
        histIter=_pedestalHists.emplace(std::make_pair(barIndex,SiPM),tfs->make<TH1F>(Form("crvPedestalHist_%i_%i",barIndex,SiPM),
                                                                                      Form("crvPedestalHist_%i_%i",barIndex,SiPM),
                                                                                      2000,-50,150)).first;
      }
      auto minmaxTest = std::minmax_element(iter->GetADCs().begin(),iter->GetADCs().end());
      if(*minmaxTest.second-*minmaxTest.first>_maxADCspread) continue;
      for(size_t i=0; i<CrvDigi::NSamples; ++i)
      {
        histIter->second->Fill(iter->GetADCs()[i]);
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvPedestalFinder;
DEFINE_ART_MODULE(CrvPedestalFinder)
