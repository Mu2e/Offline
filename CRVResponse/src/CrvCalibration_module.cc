//
// A module that finds the calibration constants for Crv channels
//
//
// Original Author: Ralf Ehrlich

#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"

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
  class CrvCalibration : public art::EDAnalyzer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("module label of the input CrvRecoPulses")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvCalibration(const Parameters& config);
    void analyze(const art::Event& e);
    void endJob();

    private:
    std::string _crvRecoPulsesModuleLabel;

    std::map<std::pair<int,int>,TH1F*> _calibHists;
  };


  CrvCalibration::CrvCalibration(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel())
  {
  }

  void CrvCalibration::endJob()
  {
  }

  void CrvCalibration::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection)) return;

    art::ServiceHandle<art::TFileService> tfs;
    for(auto iter=crvRecoPulseCollection->begin(); iter!=crvRecoPulseCollection->end(); ++iter)
    {
      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      auto histIter = _calibHists.find(std::make_pair(barIndex,SiPM));
      if(histIter==_calibHists.end())
      {
        histIter = _calibHists.emplace(std::make_pair(barIndex,SiPM),tfs->make<TH1F>(Form("crvCalibrationHist_%i_%i",barIndex,SiPM),
                                                                                     Form("crvCalibrationHist_%i_%i",barIndex,SiPM),
                                                                                     150,0,3000)).first;
      }
      if(iter->GetRecoPulseFlags().none()) histIter->second->Fill(iter->GetPulseBeta()*iter->GetPulseHeight()*TMath::E());
    }
  }

} // end namespace mu2e

using mu2e::CrvCalibration;
DEFINE_ART_MODULE(CrvCalibration)
