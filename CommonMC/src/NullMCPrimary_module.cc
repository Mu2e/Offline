//
//  Trivial module to insert an empty PrimaryParticle object and related empty DetStep collections in the event, used for NoPrimary production
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/types/Sequence.h"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include <vector>
#include <string>
namespace mu2e {
  class NullMCPrimary : public art::EDProducer {
    public:
      struct Config {
        fhicl::Sequence<std::string> extraSteps { fhicl::Name("ExtraSteps"), fhicl::Comment("2ndary keys of extra StepPointMC collections") };
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit NullMCPrimary(const Parameters& conf);
      void produce(art::Event& evt) override;
    private:
      std::vector<std::string> extrasteps_;
  };

  NullMCPrimary::NullMCPrimary(const Parameters& config )  :
    art::EDProducer{config},
    extrasteps_(config().extraSteps())
    {
      produces <PrimaryParticle>();
      produces <StatusG4>();
      produces <GenParticleCollection>();
      produces <SimParticleCollection>();
      produces <MCTrajectoryCollection>();
      for(auto extrastep : extrasteps_)
        produces <StepPointMCCollection>(extrastep);
      produces <StrawGasStepCollection>();
      produces <CrvStepCollection>();
      produces <CaloShowerStepCollection>();
      produces <SurfaceStepCollection>();
    }

  void NullMCPrimary::produce(art::Event& event) {
    // create empty output objects
    PrimaryParticle pp;
    MCTrajectoryCollection mctc;
    event.put(std::make_unique<PrimaryParticle>(pp));
    event.put(std::make_unique<StatusG4>());
    event.put(std::make_unique<GenParticleCollection>());
    event.put(std::make_unique<SimParticleCollection>());
    event.put(std::make_unique<MCTrajectoryCollection>(mctc));
    StepPointMCCollection empty;
    for(auto extrastep : extrasteps_)
      event.put(std::make_unique<StepPointMCCollection>(empty),extrastep);
    event.put(std::make_unique<StrawGasStepCollection>());
    event.put(std::make_unique<CrvStepCollection>());
    event.put(std::make_unique<CaloShowerStepCollection>());
    event.put(std::make_unique<SurfaceStepCollection>());
  }
}
DEFINE_ART_MODULE(mu2e::NullMCPrimary)
