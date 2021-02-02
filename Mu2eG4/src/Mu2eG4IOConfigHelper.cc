// Andrei Gaponenko, 2021

#include "Mu2eG4/inc/Mu2eG4IOConfigHelper.hh"

#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"

#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Core/ConsumesCollector.h"

namespace mu2e {

  Mu2eG4IOConfigHelper::Mu2eG4IOConfigHelper(const Mu2eG4Config::Top& conf,
                                             art::ProducesCollector& pc,
                                             art::ConsumesCollector& cc)
    : generatorModuleLabel_{conf.generatorModuleLabel()}
    , multiStagePars_{conf}
    , trajectoryControl_(conf.TrajectoryControl())
    , timeVD_enabled_(conf.SDConfig().TimeVD().enabled())
    , extMonPixelsEnabled_{false}
    , mu2elimits_{conf.ResourceLimits()}
    , stackingCutsConf_{conf.Mu2eG4StackingOnlyCut.get<fhicl::ParameterSet>()}
    , steppingCutsConf_{conf.Mu2eG4SteppingOnlyCut.get<fhicl::ParameterSet>()}
    , commonCutsConf_{conf.Mu2eG4CommonCut.get<fhicl::ParameterSet>()}
  {
    // Use temporary local objects to parse config and declare i/o below
    SensitiveDetectorHelper sd(conf.SDConfig());
    extMonPixelsEnabled_ = sd.extMonPixelsEnabled();

    const art::InputTag invalid_tag;

    if((generatorModuleLabel_ == invalid_tag) && !multiStagePars_.multiStage()) {
      throw cet::exception("CONFIG")
        << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
    }

    //----------------------------------------------------------------
    // Declare which products we will read.
    auto const& inputPhysVolTag = multiStagePars_.inputPhysVolumeMultiInfo();
    if (inputPhysVolTag != invalid_tag) {
      cc.consumes<PhysicalVolumeInfoMultiCollection, art::InSubRun>(inputPhysVolTag);
    }
    auto const& inputSimParticlesTag = multiStagePars_.inputSimParticles();
    if (inputSimParticlesTag != invalid_tag) {
      cc.consumes<SimParticleCollection>(inputSimParticlesTag);
    }
    auto const& inputMCTrajectoryTag = multiStagePars_.inputMCTrajectories();
    if (inputMCTrajectoryTag != invalid_tag) {
      cc.consumes<MCTrajectoryCollection>(inputMCTrajectoryTag);
    }
    if (generatorModuleLabel_ != invalid_tag) {
      cc.consumes<GenParticleCollection>(generatorModuleLabel_);
    }
    for (auto const& tag : multiStagePars_.genInputHits()) {
      cc.consumes<StepPointMCCollection>(tag);
    }

    //----------------------------------------------------------------
    pc.produces<StatusG4>();
    pc.produces<SimParticleCollection>();
    pc.produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    if(multiStagePars_.multiStage()) {
      pc.produces<SimParticleRemapping>();
    }

    if(timeVD_enabled_) {
      static const StepInstanceName timeVD(StepInstanceName::timeVD);
      pc.produces<StepPointMCCollection>(timeVD.name());
    }

    if(trajectoryControl_.produce()) {
      pc.produces<MCTrajectoryCollection>();
    }

    // Use temporary local objects to call the declare methods
    sd.declareProducts(pc); // FIXME: also handle SD consumes

    auto stackingCuts{createMu2eG4Cuts(stackingCutsConf_, mu2elimits_)};
    stackingCuts->declareProducts(pc, cc);

    auto steppingCuts{createMu2eG4Cuts(steppingCutsConf_, mu2elimits_)};
    steppingCuts->declareProducts(pc, cc);

    auto commonCuts{createMu2eG4Cuts(commonCutsConf_, mu2elimits_)};
    commonCuts->declareProducts(pc, cc);
  }

} // end namespace mu2e
