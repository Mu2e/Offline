//
//  This module creates the StrawGasStep objects used in downstream digi simulation, using the
//  G4 StepPointMCs
//
//  Original author: David Brown (LBNL), Krzysztof Genser 19 Aug. 2019
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
#include "TTree.h"

#include "Offline/RecoDataProducts/inc/STMHit.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PlotSTMSpectrum : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<art::InputTag> stmHitsTag{ Name("stmHitsTag"), Comment("InputTag for STMHitCollection")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMSpectrum(const Parameters& conf);

    private:
    void beginJob() override; 
      void analyze(const art::Event& e) override;

    art::InputTag _stmHitsTag;
    TH1D* _energySpectrum;
  };

  PlotSTMSpectrum::PlotSTMSpectrum(const Parameters& config )  : 
    art::EDAnalyzer{config},
    _stmHitsTag(config().stmHitsTag())
  {
    consumes<STMHitCollection>(_stmHitsTag);
  }

  void PlotSTMSpectrum::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create TTree
    _energySpectrum=tfs->make<TH1D>("energySpectrum", "Energy Spectrum", 2000,0,2.0);
  }

  void PlotSTMSpectrum::analyze(const art::Event& event) {

    auto hitsHandle = event.getValidHandle<STMHitCollection>(_stmHitsTag);

    for (const auto& hit : *hitsHandle) {
      //      float time = hit.time();
      float energy = hit.energy();
      _energySpectrum->Fill(energy);
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMSpectrum)
