// Ntuple dumper for StepPointMCs.
//
// Andrei Gaponenko, 2013

#include <string>
#include <vector>
#include <sstream>

#include "cetlib_except/exception.h"

#include "TH1.h"
#include "TH2.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {

  //================================================================
  class G4ReactionAnalyzer : public art::EDAnalyzer {
    art::InputTag inputs_;
    int pdgId_;

    TH1* hStoppingCodes_;
    TH1* hDaugherCreationCodes_;
    TH2* hDaugherMultiplicity_;

    art::ServiceHandle<art::TFileService> tfs() { return art::ServiceHandle<art::TFileService>(); }

  public:
    explicit G4ReactionAnalyzer(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  G4ReactionAnalyzer::G4ReactionAnalyzer(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , inputs_(pset.get<std::string>("inputs"))
    , pdgId_(pset.get<int>("pdgId"))
    , hStoppingCodes_(tfs()->make<TH1D>("stoppingCodes", "stopping codes", 1, 0., 1.))
    , hDaugherCreationCodes_(tfs()->make<TH1D>("creationCodes", "daughter creation codes", 1, 0., 1.))
    , hDaugherMultiplicity_(tfs()->make<TH2D>("multiplicity", "daughter multiplicity", 1, 0., 1., 25, 0.5, 25.5))
  {
    hDaugherMultiplicity_->SetOption("colz");
  }

  //================================================================
  void G4ReactionAnalyzer::analyze(const art::Event& event) {
    typedef std::map<std::string, int> StringStats;

    const auto& coll = event.getValidHandle<SimParticleCollection>(inputs_);
    for(const auto& i : *coll) {
      const auto& particle = i.second;
      if(particle.pdgId() == pdgId_) {
        hStoppingCodes_->Fill(particle.stoppingCode().name().c_str(), 1.);
        StringStats ccount;
        for(const auto& daughter: particle.daughters()) {
          std::ostringstream pp;
          if(static_cast<int>(daughter->pdgId()) < 1000000000) {
            pp<<daughter->pdgId()<<" "<<daughter->creationCode();
          }
          else {
            pp<<"nucleus "<<daughter->creationCode();
          }

          hDaugherCreationCodes_->Fill(pp.str().c_str(), 1.);
          ++ccount[pp.str()];
        }
        for(const auto& i : ccount) {
          hDaugherMultiplicity_->Fill(i.first.c_str(), i.second, 1.);
        }
      }
    }
  } // analyze(event)

    //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::G4ReactionAnalyzer);
