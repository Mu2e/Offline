// TrackSummaryDataAnalyzer: no MC dependencies.
//
// Andrei Gaponenko, 2014

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <limits>

#include "cetlib_except/exception.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "RecoDataProducts/inc/TrackSummary.hh"
#include "Mu2eUtilities/inc/HistTrackSum.hh"
#include "Mu2eUtilities/inc/TrackCuts.hh"

#include "TH1.h"

namespace mu2e {

  class TrackSummaryDataAnalyzer : public art::EDAnalyzer {
  public:
    explicit TrackSummaryDataAnalyzer(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    TrackCuts an_;
    art::InputTag trackInput_;
  };

  //================================================================
  TrackSummaryDataAnalyzer::TrackSummaryDataAnalyzer(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , an_(pset.get<fhicl::ParameterSet>("cuts"), *art::ServiceHandle<art::TFileService>(), "")
    , trackInput_(pset.get<std::string>("trackInput"))
  {}

  //================================================================
  void TrackSummaryDataAnalyzer::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<TrackSummaryCollection>(trackInput_);
    for(unsigned i=0; i<ih->size(); ++i) {
      const auto& track = ih->at(i);
      const double weight = 1.; // no weights for real data
      if(an_.accepted(track, weight)) {
        // Fill more histograms for selected tracks here...
      }
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackSummaryDataAnalyzer);
