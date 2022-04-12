//
// Analyzer module to create a histogram of the STMDigi energies
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

#include "Offline/RecoDataProducts/inc/STMDigi.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PlotSTMDigis : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<art::InputTag> stmDigisTag{ Name("stmDigisTag"), Comment("InputTag for STMDigiCollection")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMDigis(const Parameters& conf);

    private:
    void beginJob() override; 
      void analyze(const art::Event& e) override;

    void endJob() override;
 
    art::InputTag _stmDigisTag;
    TH1D* _energySpectrum;
  };

  PlotSTMDigis::PlotSTMDigis(const Parameters& config )  : 
    art::EDAnalyzer{config},
    _stmDigisTag(config().stmDigisTag())
  {
    consumes<STMDigiCollection>(_stmDigisTag);
  }

  void PlotSTMDigis::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    // create TTree
    _energySpectrum=tfs->make<TH1D>("energySpectrum", "Energy Spectrum", 5000,0,5000.0);
  }

  void PlotSTMDigis::analyze(const art::Event& event) {

    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);

    for (const auto& digi : *digisHandle) {
      //      float time = hit.time();
      float baselineMean = digi.baselineMean();
      _energySpectrum->Fill(baselineMean);
    }
  }
  void PlotSTMDigis::endJob() {
    // Insert pain (fits) here
    // float peaks[] = [0.123,0.162,0.245,0.344]; 
    // double peak_energy = 0.344;
    // int bin_number = _energySpectrum->GetXaxis()->FindBin(peak_energy);
    // float bin_content = _energyContent->GetBinContent(bin_number);
    //TF1* fitGaus = new TF1("fit_Gaus", "[0]*TMath::Gaus(x,[1],[2])",0.340,0.350);
    //fitGaus->SetParameters(1700,0.344,0.001);
    //_energySpectrum->Fit(fitGaus,"R");
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMDigis)
