//
// Analyzer module to create a histogram of the STMHit energies
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
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "Offline/RecoDataProducts/inc/STMHit.hh"

// c++ headers
#include <fstream>

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

    void endJob() override;

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
  void PlotSTMSpectrum::endJob() {
    // Insert pain (fits) here
    double peaks[4] = {0.123,0.162,0.245,0.344};
    // double peak_energy = 0.344;
    //int bin_number = _energySpectrum->GetXaxis()->FindBin(peak_energy);
    //float bin_content = _energyContent->GetBinContent(bin_number);
    ofstream out;
    out.open("out.log", ios::out | ios::trunc);
    out << "par_name, par_value, par_error" << std::endl;
    for (int j = 0; j < 4; j++)
      {
       TString fname(Form("fgaus_%d", j));
       TF1* fitGaus = new TF1(fname, "[0]*TMath::Gaus(x,[1],[2])",peaks[j]-0.005,peaks[j]+0.005);
       fitGaus->SetParameters(1700,peaks[j],0.001);
       _energySpectrum->Fit(fitGaus,"R+");
       // Export fit parameters to a file
       // Write out fit parameters to screen separate from the MSG
       TFitResultPtr fitresult = _energySpectrum->Fit(fitGaus, "QRS");
       int n_par = fitresult->NPar();

       for (int i_par = 0; i_par < n_par; ++i_par)
         {
           out   << j << ", "
                 << fitresult->GetParameterName(i_par)
                 << ", " << fitresult->Parameter(i_par)
                 << ", " << fitresult->ParError(i_par)
                 << std::endl;
         }
      }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMSpectrum)
