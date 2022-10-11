//
// A module that finds the calibration constants for Crv channels
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
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
    void beginRun(const art::Run&);

    private:
    std::string        _crvRecoPulsesModuleLabel;
    std::vector<TH1F*> _calibHists;
  };


  CrvCalibration::CrvCalibration(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel())
  {
  }

  void CrvCalibration::beginRun(const art::Run&)
  {
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    _calibHists.reserve(counters.size()*4);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t barIndex=0; barIndex<counters.size(); ++barIndex)
    {
      for(size_t SiPM=0; SiPM<4; ++SiPM)
      {
        //produce histograms also for non-existing channels to get a continuously running index
        size_t channelIndex=barIndex*4 + SiPM;
        _calibHists.emplace_back(tfs->make<TH1F>(Form("crvCalibrationHist_%lu",channelIndex),
                                                 Form("crvCalibrationHist_%lu",channelIndex),
                                                 150,0,3000));
      }
    }
  }

  void CrvCalibration::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection)) return;

    for(auto iter=crvRecoPulseCollection->begin(); iter!=crvRecoPulseCollection->end(); ++iter)
    {
      if(!iter->GetRecoPulseFlags().none()) continue;

      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      int channelIndex=barIndex*4+SiPM;
      _calibHists.at(channelIndex)->Fill(iter->GetPulseBeta()*iter->GetPulseHeight()*TMath::E());
    }
  }

} // end namespace mu2e

using mu2e::CrvCalibration;
DEFINE_ART_MODULE(CrvCalibration)
