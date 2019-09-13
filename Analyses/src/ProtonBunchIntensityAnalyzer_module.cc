// Histogram a proton bunch intensity distribution.
//
// Andrei Gaponenko, 2018

#include "fhiclcpp/types/Atom.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
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
      fhicl::Atom<art::InputTag> meanPBItag{Name("MeanBeamIntensity"), Comment("Tag for MeanBeamIntensity"), art::InputTag()};

      fhicl::Atom<unsigned> nbins{Name("nbins"), Comment("Number of bins in the histogram"), 250u};
      fhicl::Atom<double> hmin{Name("hmin"), Comment("Histogram min"), 0.};
      fhicl::Atom<double> hmax{Name("hmax"), Comment("Absolute PBI histogram max"), 1.e8};
      fhicl::Atom<double> rmax{Name("rmax"), Comment("Relative PBI histogram max"), 3.0};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit ProtonBunchIntensityAnalyzer(const Parameters& conf);
    void analyze(const art::Event& evt) override;
    void beginSubRun(const art::SubRun& subrun) override;
  private:
    Config conf_;
    TH1 *hh_, *hhr_;
    double meanPBI_;
  };

  //================================================================
  ProtonBunchIntensityAnalyzer::ProtonBunchIntensityAnalyzer(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , conf_(conf())
    , hh_(art::ServiceHandle<art::TFileService>()->
          make<TH1D>("pbi", "Absolute proton bunch intensity",
                     conf().nbins(), conf().hmin(), conf().hmax()))
    , hhr_(art::ServiceHandle<art::TFileService>()->
          make<TH1D>("rpbi", "Relative proton bunch intensity",
                     conf().nbins(), 0.0, conf().rmax()))
    , meanPBI_(-1.0)
  {
    hh_->StatOverflows();
    hhr_->StatOverflows();
  }

  void ProtonBunchIntensityAnalyzer::beginSubRun(const art::SubRun & subrun ) {
    // mean number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    subrun.getByLabel(conf_.meanPBItag(), PBIHandle);
    if(PBIHandle.isValid())
      meanPBI_ = PBIHandle->intensity();
  }

  //================================================================
  void ProtonBunchIntensityAnalyzer::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<ProtonBunchIntensity>(conf_.input());
    hh_->Fill(ih->intensity());
    if(meanPBI_>0.0)hhr_->Fill(ih->intensity()/meanPBI_);
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensityAnalyzer);
