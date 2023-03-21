//
// Analyzer module to create a histogram of the STMMWDDigi uncalibrated energies
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
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TGraph.h"

#include "Offline/RecoDataProducts/inc/STMHit.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PlotSTMEnergySpectrum : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmHitsTag{ Name("stmHitsTag"), Comment("InputTag for STMHitCollection")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMEnergySpectrum(const Parameters& conf);

    private:
    void beginJob() override;
    void analyze(const art::Event& e) override;

    TH1D* _energySpectrum;
    art::ProductToken<STMHitCollection> _stmHitsToken;
  };

  PlotSTMEnergySpectrum::PlotSTMEnergySpectrum(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmHitsToken(consumes<STMHitCollection>(config().stmHitsTag()))
  { }

  void PlotSTMEnergySpectrum::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create histograms
    double min_energy = 0;
    double max_energy = 10;
    double energy_bin_width = 0.001;
    int n_bins = (max_energy - min_energy) / energy_bin_width;
    _energySpectrum=tfs->make<TH1D>("energySpectrum", "Energy Spectrum", n_bins, min_energy, max_energy);
  }

  void PlotSTMEnergySpectrum::analyze(const art::Event& event) {

    auto stmHitsHandle = event.getValidHandle(_stmHitsToken);

    for (const auto& stmHit : *stmHitsHandle) {
      auto energy = stmHit.energy();
      _energySpectrum->Fill(energy);
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMEnergySpectrum)
