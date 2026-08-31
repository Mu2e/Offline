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
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include <utility>
#include <map>
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
        fhicl::Atom<art::InputTag> stmPHDigisTag{ Name("stmPHDigisTag"), Comment("InputTag for STMPHDigiCollectionMap")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMPHSpectrum(const Parameters& conf);

    private:
    void beginJob() override;
    void analyze(const art::Event& e) override;

    TH2F* _twoDhist; //Histograms of Energy vs binned event
    int artEventCount = 0;

    TH1D* _phSpectrum;
    art::ProductToken<STMPHDigiCollectionMap> _stmPHDigisMapToken;
    STMChannel _channel;
  };

  PlotSTMPHSpectrum::PlotSTMPHSpectrum(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmPHDigisMapToken(consumes<STMPHDigiCollectionMap>(config().stmPHDigisTag())),
    _channel(STMUtils::getChannel(config().stmPHDigisTag()))
  { }

  void PlotSTMPHSpectrum::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create histograms
    std::string phSpectrumTitle = "PH Spectrum (" +_channel.name() + ")"; // Builds title for PH Spectrum
    _phSpectrum=tfs->make<TH1D>("phSpectrum", (phSpectrumTitle + ";Pulse Height;Count").c_str(), 10000, -10000, 0); //bins,min,max

    _twoDhist=tfs->make<TH2F>("phEvent",("Pulse Height vs Art Events (" + _channel.name() + ");Event Bins; Pulse Height").c_str(), // (name, title;xtitle;ytitle, nbinsX, xlow, xup, nbinsY, ylow, yup)
                              10,0,10,     // X-axis scale
                              10000,-10000,0);   // Y-axis scale
  }

  void PlotSTMPHSpectrum::analyze(const art::Event& event) {
    // get map handle
    auto phDigisMapHandle = event.getValidHandle(_stmPHDigisMapToken);
    const auto& phDigiMap = *phDigisMapHandle;
    int binBlock = artEventCount/100;
    for (const auto& i_phDigiMap : phDigiMap){
        // get PH Digi Collection
        const auto& phDigis = i_phDigiMap.second;
        for (const auto& phDigi : phDigis) {
            // get uncalibrated energy (adc)
            auto energy = phDigi.energy();
            _phSpectrum->Fill(energy);
            _twoDhist->Fill(binBlock, energy);
            }
        }
    ++artEventCount;
    }
}

DEFINE_ART_MODULE(mu2e::PlotSTMPHSpectrum)
