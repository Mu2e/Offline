// TrackSummaryMCAnalyzer: start with tracks.
// Use MC weights and MC process info to evaluate
// backgrounds without double counting.
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
#include "canvas/Persistency/Common/FindOneP.h"

#include "RecoDataProducts/inc/TrackSummary.hh"
#include "MCDataProducts/inc/TrackSummaryTruthAssns.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"

#include "Mu2eUtilities/inc/HistTrackSum.hh"
#include "Mu2eUtilities/inc/TrackCuts.hh"

#include "TH1.h"

namespace mu2e {

  class TrackSummaryMCAnalyzer : public art::EDAnalyzer {
  public:
    explicit TrackSummaryMCAnalyzer(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    TrackCuts signalAnalysis_;
    TrackCuts otherAnalysis_;
    TrackCuts fakeAnalysis_;

    art::InputTag trackInput_;
    art::InputTag truthMapInput_;
    art::InputTag signalSimParticleCollection_;

    std::vector<art::InputTag> eventWeights_;
  };

  //================================================================
  TrackSummaryMCAnalyzer::TrackSummaryMCAnalyzer(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , signalAnalysis_(pset.get<fhicl::ParameterSet>("cuts"), *art::ServiceHandle<art::TFileService>(), "signal")
    , otherAnalysis_(pset.get<fhicl::ParameterSet>("cuts"), *art::ServiceHandle<art::TFileService>(), "other")
    , fakeAnalysis_(pset.get<fhicl::ParameterSet>("cuts"), *art::ServiceHandle<art::TFileService>(), "fake")
    , trackInput_(pset.get<std::string>("trackInput"))
    , truthMapInput_(pset.get<std::string>("truthMapInput"))
    , signalSimParticleCollection_(pset.get<std::string>("signalSimParticleCollection"))
    , eventWeights_(pset.get<std::vector<art::InputTag> >("eventWeights"))
  {}

  //================================================================
  void TrackSummaryMCAnalyzer::analyze(const art::Event& event) {
    // Can't use ValidHandle: redmine #6422
    // auto ih = event.getValidHandle<TrackSummaryCollection>(trackInput_);
    art::Handle<TrackSummaryCollection> ih;
    event.getByLabel(trackInput_, ih);

    art::Handle<SimParticleCollection> signalHandle;
    event.getByLabel(signalSimParticleCollection_, signalHandle);

    for(unsigned i=0; i<ih->size(); ++i) {
      const auto& track = ih->at(i);

      // To evaluate backgrounds without double counting we should
      // accept only tracks corresponding to the simulated process of
      // interest.

      art::FindOneP<SimParticle, mu2e::TrackSummaryMatchInfo> f1(ih, event, truthMapInput_);
      if(f1.size() == 1) {
        const art::Ptr<SimParticle>& sim = f1.at(i);

        if(sim.id() == signalHandle.id()) {

          // Use the specified weights.
          double weight = 1.;
          for (const auto& iwt : eventWeights_ ) {
            weight *= event.getValidHandle<EventWeight>(iwt)->weight();
          }

          if(signalAnalysis_.accepted(track, weight)) {
            // Fill more histograms for selected tracks...
          }

        }
        else {
          // Known MC source, different than the current process of interest
          otherAnalysis_.accepted(track, 1.);
        }
      }
      else {
        // Fake tracks. (Or no MC info.)
        fakeAnalysis_.accepted(track, 1.);
      }

    } // for(tracks)
  } // analyze()

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::TrackSummaryMCAnalyzer);
