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
    TH1D* _baselineMean;
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
    _baselineMean=tfs->make<TH1D>("baselineMean", "Energy Spectrum", 5000,0,5000.0);
  }

  void PlotSTMDigis::analyze(const art::Event& event) {
    art::ServiceHandle<art::TFileService> tfs;
    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);
    int j = 0;
    for (const auto& digi : *digisHandle) {
      float baselineMean = digi.baselineMean();
      _baselineMean->Fill(baselineMean);
      TString fname(Form("digiSpectrum_%d_%d", event.event(), j));
      // Create a histogram
      TH1D* hWaveform = tfs->make<TH1D>(fname, "Digi Spectrum", digi.adcs().size(), 0, digi.adcs().size());
      // Loop through the adcs
      int i_bin = 1;
      for (const auto& sample : digi.adcs())
         {
           hWaveform->SetBinContent(i_bin, sample);
           i_bin++;
         }
      j++;
    }
  }
  void PlotSTMDigis::endJob() {
    // Insert pain (fits) here
    // Goal is to automatically find peaks, these peaks are from the Eu152 source, so...
    /*
      double peaks[4] = [0.123,0.162,0.245,0.344];
    // double peak_energy = 0.344;
      for (int j = 0; j < 4; j++)
      {
       TString fname(Form("fgaus_%d", j));
       TF1* fitGaus = new TF1(fname, "[0]*TMath::Gaus(x,[1],[2])",peaks[j]-0.005,peaks[j]+0.005);
       fitGaus->SetParameters(1700,peaks[j],0.001);
       _baselineMean->Fit(fitGaus,"R+");
      }
    */
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMDigis)
