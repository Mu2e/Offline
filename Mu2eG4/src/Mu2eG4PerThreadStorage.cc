// Mu2eG4PerThreadStorage holds the pointers to the UserActions and the various
// Manager classes, and other thread local objects for running G4 in MT mode.
//
// Author: Lisa Goodeough
// Date: 2019/07/09
//

#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"

namespace mu2e {

  Mu2eG4PerThreadStorage::Mu2eG4PerThreadStorage(const art::InputTag &genLabel)
    : generatorModuleLabel(genLabel)
    , tvd{"", nullptr}
 {}

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  initializeEventInfo(art::Event* evt,
                      SimParticleHelper* sim_part_helper,
                      SimParticlePrimaryHelper* sim_part_primary_helper,
                      HitHandles* gen_input_hits) {
    artEvent = evt;
    simParticleHelper = sim_part_helper;
    simParticlePrimaryHelper = sim_part_primary_helper;
    genInputHits = gen_input_hits;

    if(!(generatorModuleLabel == art::InputTag())) {
      artEvent->getByLabel(generatorModuleLabel, gensHandle);
    }

    if ( !gensHandle.isValid() && genInputHits == nullptr ) {
      throw cet::exception("CONFIG")
        << "Error in PerThreadStorage::initializeEventInfo.  You are trying to run a G4job without an input for G4.\n";
    }
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  putSensitiveDetectorData(art::EDProductGetter const* sim_product_getter) {

    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > steps_map =
      std::move(sensitiveDetectorSteps);

    for (std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = steps_map.begin();
         i != steps_map.end(); ++i) {

      for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){

        StepPointMC& step = *j;

        if ( step.simParticle().isNonnull() ){
          step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                     step.simParticle().key(),
                                                     sim_product_getter );
        }//if
      }//for StepPointMCCollection::iterator

      artEvent->put(std::move(i->second), i->first);
    }//for (std::unordered_map...
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  putCutsData(art::EDProductGetter const* sim_product_getter) {

    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > cuts_map = std::move(cutsSteps);

    for (std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = cuts_map.begin();
         i != cuts_map.end(); ++i) {

      for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){
        StepPointMC& step = *j;

        if ( step.simParticle().isNonnull() ){
          step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                     step.simParticle().key(),
                                                     sim_product_getter );
        }//if
      }//for StepPointMCCollection::iterator

      artEvent->put(std::move(i->second), i->first);
    }
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::clearData() {

    artEvent = nullptr;
    simParticleHelper = nullptr;
    simParticlePrimaryHelper = nullptr;
    genInputHits = nullptr;
    gensHandle.clear();

    statG4 = nullptr;
    simPartCollection = nullptr;
    tvd.first = "";
    tvd.second = nullptr;
    mcTrajectories = nullptr;
    simRemapping = nullptr;
    extMonFNALHits = nullptr;
    sensitiveDetectorSteps.clear();
    cutsSteps.clear();
  }

  //----------------------------------------------------------------

} // end namespace mu2e
