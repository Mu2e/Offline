//
// Analyzer module to create a histogram of the STMPHDigi uncalibrated energies
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include <utility>
// root
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TGraph.h"

#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PlotSTMPHSpectrum : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmPHDigisTag{ Name("stmPHDigisTag"), Comment("InputTag for STMPHDigiCollection")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMPHSpectrum(const Parameters& conf);

    private:
    void beginJob() override;
    void analyze(const art::Event& e) override;

    TH2F* _twoDhist; //Histograms of Energy vs binned event
    int eventCount = 0;
    
    TH1D* _phSpectrum;
    art::ProductToken<STMPHDigiCollection> _stmPHDigisToken;
  };

  PlotSTMPHSpectrum::PlotSTMPHSpectrum(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmPHDigisToken(consumes<STMPHDigiCollection>(config().stmPHDigisTag()))
  { }

  void PlotSTMPHSpectrum::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create histograms
    _phSpectrum=tfs->make<TH1D>("phSpectrum", "PH Spectrum", 1000, 0, 1e4);
    _twoDhist=tfs->make<TH2F>("twoDhist","Pulse Height vs events;Event Bins; Pulse Height",
			    1000,0,1000,     // X-axis scale
			    1000,0,1e5);   // Y-axis scale
  }

  void PlotSTMPHSpectrum::analyze(const art::Event& event) {

    auto phDigisHandle = event.getValidHandle(_stmPHDigisToken);
    int binBlock = eventCount/100;
    
    for (const auto& phDigi : *phDigisHandle) {
      auto energy = phDigi.energy();
      _phSpectrum->Fill(energy);
      _twoDhist->Fill(binBlock, energy);
    }
    ++eventCount;
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMPHSpectrum)
