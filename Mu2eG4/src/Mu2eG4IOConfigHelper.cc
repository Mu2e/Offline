// Andrei Gaponenko, 2021

#include "Offline/Mu2eG4/inc/Mu2eG4IOConfigHelper.hh"

#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Offline/Mu2eG4/inc/IMu2eG4Cut.hh"

#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/ScorerSummary.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticleRemapping.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

#include "art/Framework/Core/ProducesCollector.h"
#include "art/Framework/Core/ConsumesCollector.h"

#include "cetlib_except/exception.h"

namespace mu2e {

  Mu2eG4IOConfigHelper::Mu2eG4IOConfigHelper(const Mu2eG4Config::Top& conf,
                                             art::ProducesCollector& pc,
                                             art::ConsumesCollector& cc)
    : inputs_{conf.inputs()}
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

    switch(inputs_.primaryType().id()) {
    default: throw cet::exception("CONFIG")
        << "Error: unknown Mu2eG4 primaryType id = "<<inputs_.primaryType().id()<<std::endl;

    case Mu2eG4PrimaryType::GenParticles:
      cc.consumes<GenParticleCollection>(inputs_.primaryTag());
      break;

    case Mu2eG4PrimaryType::StepPoints:
      cc.consumes<StepPointMCCollection>(inputs_.primaryTag());
      break;

    case Mu2eG4PrimaryType::SimParticleLeaves:
      cc.consumes<SimParticleCollection>(inputs_.primaryTag());
      break;

    case Mu2eG4PrimaryType::StageParticles:
      cc.consumes<StageParticleCollection>(inputs_.primaryTag());
      break;
    }

    //----------------------------------------------------------------
    // Declare which products we will read.
    const art::InputTag invalid_tag;
    auto const& inputPhysVolTag = inputs_.inputPhysVolumeMultiInfo();
    if (inputPhysVolTag != invalid_tag) {
      cc.consumes<PhysicalVolumeInfoMultiCollection, art::InSubRun>(inputPhysVolTag);
    }

    if(inputs_.updateEventLevelVolumeInfos()) {
      cc.consumes<PhysicalVolumeInfoMultiCollection>(inputs_.updateEventLevelVolumeInfos()->input);
      pc.produces<PhysicalVolumeInfoMultiCollection>(inputs_.updateEventLevelVolumeInfos()->outInstance);
    }

    auto const& inputMCTrajectoryTag = inputs_.inputMCTrajectories();
    if (inputMCTrajectoryTag != invalid_tag) {
      cc.consumes<MCTrajectoryCollection>(inputMCTrajectoryTag);
    }

    //----------------------------------------------------------------
    pc.produces<StatusG4>();
    pc.produces<SimParticleCollection>();
    pc.produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    if(inputs_.multiStage()) {
      pc.produces<SimParticleRemapping>();
    }

    if(timeVD_enabled_) {
      static const StepInstanceName timeVD(StepInstanceName::timeVD);
      pc.produces<StepPointMCCollection>(timeVD.name());
    }

    if(trajectoryControl_.produce()) {
      pc.produces<MCTrajectoryCollection>();
    }

    //----------------------------------------------------------------
    auto const& scoringTable = conf.scoring();
    if (scoringTable.enabled()) {
      for (const auto& meshName : scoringTable.meshNames()){
        for (const auto& scorerName: scoringTable.scorerNames()){
          std::string instanceName = meshName+scorerName;
          pc.produces<ScorerSummaryCollection,art::InSubRun>(instanceName);
        }
      }
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
