// Histogram time offsets in a time map object.
//
// Andrei Gaponenko, 2016

#include "fhiclcpp/types/Atom.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "MCDataProducts/inc/SimParticleTimeMap.hh"

#include "TH1.h"

namespace mu2e {

  class SimParticleTimeMapAnalyzer : public art::EDAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> input{ Name("input"), Comment("Tag of the SimParticleTimeMap to analyze.")};

      fhicl::Atom<unsigned> nbins{Name("nbins"), Comment("Number of bins in the histogram"), 250u};
      fhicl::Atom<double> tmin{Name("tmin"), Comment("Histogram tmin"), -500.};
      fhicl::Atom<double> tmax{Name("tmax"), Comment("Histogram tmax"), +2000.};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit SimParticleTimeMapAnalyzer(const Parameters& conf);
    void analyze(const art::Event& evt) override;
  private:
    Config conf_;
    TH1 *histTime_;
  };

  //================================================================
  SimParticleTimeMapAnalyzer::SimParticleTimeMapAnalyzer(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , conf_(conf())
    , histTime_(art::ServiceHandle<art::TFileService>()->
                make<TH1D>("timeOffset", "Time offset",
                           conf().nbins(), conf().tmin(), conf().tmax()))
  {}

  //================================================================
  void SimParticleTimeMapAnalyzer::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<SimParticleTimeMap>(conf_.input());
    for(const auto& entry: *ih) {
      histTime_->Fill(entry.second);
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleTimeMapAnalyzer);
