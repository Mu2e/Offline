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

#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PlotSTMMWDSpectrum : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmMWDDigisTag{ Name("stmMWDDigisTag"), Comment("InputTag for STMMWDDigiCollection")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMMWDSpectrum(const Parameters& conf);

    private:
    void beginJob() override;
    void analyze(const art::Event& e) override;

    art::InputTag _stmMWDDigisTag;
    TH1D* _mwdSpectrum;
  };

  PlotSTMMWDSpectrum::PlotSTMMWDSpectrum(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmMWDDigisTag(config().stmMWDDigisTag())
  {
    consumes<STMMWDDigiCollection>(_stmMWDDigisTag);
  }

  void PlotSTMMWDSpectrum::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create histograms
    _mwdSpectrum=tfs->make<TH1D>("mwdSpectrum", "MWD Spectrum", 1000, 0, 1e4);
  }

  void PlotSTMMWDSpectrum::analyze(const art::Event& event) {

    auto mwdDigisHandle = event.getValidHandle<STMMWDDigiCollection>(_stmMWDDigisTag);

    for (const auto& mwdDigi : *mwdDigisHandle) {
      auto energy = mwdDigi.energy();
      _mwdSpectrum->Fill(energy);
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMMWDSpectrum)
