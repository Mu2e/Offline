//
// Analyzer module to create a histogram of calibrated STMHit energies
// Stand alone module for quick energy-spectrum checks
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
#include <map>
// root
#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "TGraph.h"

#include "Offline/RecoDataProducts/inc/STMHit.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/DataProducts/inc/STMChannel.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PlotSTMEnergySpectrum : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmHitsMapTag{ Name("stmHitsMapTag"), Comment("InputTag for STMHitCollectionMap")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMEnergySpectrum(const Parameters& conf);

    private:
    void beginJob() override;
    void analyze(const art::Event& e) override;

    TH1D* _energySpectrum;
    art::ProductToken<STMHitCollectionMap> _stmHitsCollectionMapToken;
    STMChannel _channel;
  };

  PlotSTMEnergySpectrum::PlotSTMEnergySpectrum(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmHitsCollectionMapToken(consumes<STMHitCollectionMap>(config().stmHitsMapTag())),
    _channel(STMUtils::getChannel(config().stmHitsMapTag()))
  { }

  void PlotSTMEnergySpectrum::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create histograms
    double min_energy = 0;
    double max_energy = 10;
    double energy_bin_width = 0.001;
    int n_bins = (max_energy - min_energy) / energy_bin_width;

    std::string energySpectrumTitle = "Energy Spectrum (" + _channel.name() + ")" ;
    _energySpectrum=tfs->make<TH1D>("energySpectrum",
                                    (energySpectrumTitle +";Energy;Count").c_str(),
                                    n_bins, min_energy, max_energy);
  }

  void PlotSTMEnergySpectrum::analyze(const art::Event& event) {

    auto stmHitsCollectionMapHandle = event.getValidHandle(_stmHitsCollectionMapToken);
    for (const auto& mu2e_evt: *stmHitsCollectionMapHandle) {
      // Can get eventHeader if you need it here

      // Get hit collection
      const auto& stmHits = mu2e_evt.second;
      for (const auto& stmHit : stmHits) {
        auto energy = stmHit.energy();
        _energySpectrum->Fill(energy);
      }
    }
  } // analyze
} // namespace

DEFINE_ART_MODULE(mu2e::PlotSTMEnergySpectrum)
