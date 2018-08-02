// Histogram a proton bunch intensity distribution.
//
// Andrei Gaponenko, 2018

#include "fhiclcpp/types/Atom.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "MCDataProducts/inc/ProtonBunchIntensity.hh"

#include "TH1.h"

namespace mu2e {

  class ProtonBunchIntensityAnalyzer : public art::EDAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> input{ Name("input"), Comment("Tag of the ProtonBunchIntensity object to analyze.")};

      fhicl::Atom<unsigned> nbins{Name("nbins"), Comment("Number of bins in the histogram"), 250u};
      fhicl::Atom<double> hmin{Name("hmin"), Comment("Histogram min"), 0.};
      fhicl::Atom<double> hmax{Name("hmax"), Comment("Histogram max"), 1.e8};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit ProtonBunchIntensityAnalyzer(const Parameters& conf);
    void analyze(const art::Event& evt) override;
  private:
    Config conf_;
    TH1 *hh_;
  };

  //================================================================
  ProtonBunchIntensityAnalyzer::ProtonBunchIntensityAnalyzer(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , conf_(conf())
    , hh_(art::ServiceHandle<art::TFileService>()->
          make<TH1D>("pbi", "Proton bunch intensity",
                     conf().nbins(), conf().hmin(), conf().hmax()))
  {
    hh_->StatOverflows();
  }

  //================================================================
  void ProtonBunchIntensityAnalyzer::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<ProtonBunchIntensity>(conf_.input());
    hh_->Fill(ih->intensity());
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensityAnalyzer);
