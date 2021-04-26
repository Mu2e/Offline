#include "Mu2eG4/inc/Mu2eG4Inputs.hh"

#include <optional>
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {
  Mu2eG4Inputs::Mu2eG4Inputs(const Mu2eG4Config::Inputs_& conf)
    : primaryType_{conf.primaryType()}
    , primaryTag_{conf.primaryTag()}
    , inputMCTrajectories_{conf.inputMCTrajectories()}
    , inputPhysVolumeMultiInfo_{conf.inputPhysVolumeMultiInfo()}
    , multiStage_{ primaryType_.id() != Mu2eG4PrimaryType::GenParticles }
  {
    unsigned simStage{0};
    if(conf.simStageOverride(simStage)) {
      simStageOverride_.emplace(simStage);
    }

    Mu2eG4Config::EventLevelVolInfos evconf;
    if(conf.updateEventLevelVolumeInfos(evconf)) {
      elvi_.emplace(EventLevelVolInfos{evconf.input(), evconf.outInstance()});
    }
  }

  //================================================================
  art::Handle<SimParticleCollection>
  Mu2eG4Inputs::inputSimParticles(const art::Event& evt) const {
    art::Handle<SimParticleCollection> res;

    switch(primaryType_.id()) {
    default: throw cet::exception("CONFIG")
        << "Error: Mu2eG4Inputs::inputSimParticles(): unknown Mu2eG4 primaryType id = "
        <<primaryType_.id()
        <<std::endl;

    case Mu2eG4PrimaryType::GenParticles:
      break; // No input sims for GenParticle jobs

    case Mu2eG4PrimaryType::StepPoints: {
      std::optional<art::ProductID> simid;

      auto const h = evt.getValidHandle<StepPointMCCollection>(primaryTag_);
      for(const auto& hit : *h) {
        if(!simid) {
          simid = hit.simParticle().id();

          // We have to read the pointed-to product into memory by hand
          // for the get(ProductID, &handle) call below to work.
          // See Kyle's e-mail on the art-users list on 2021-03-10.
          *hit.simParticle();
        }
        else {
          if(*simid != hit.simParticle().id()) {
            throw cet::exception("BADINPUT")
              <<"Mu2eG4Inputs::inputSimParticles(): inconsistent SimParticleCollection product ID "
              <<"in input StepPointMCCollection "
              <<primaryTag_
              <<std::endl;
          }
        }
      }

      if(simid) { // nothing to retrieve if primary inputs are emtpy
        evt.get(*simid, res);
        if(!res.isValid()) {
          throw cet::exception("BADINPUT")
            <<"Mu2eG4Inputs::inputSimParticles(): could not get SimParticleCollection "
            <<"pointed to by step points in "
            <<primaryTag_
            <<std::endl;
        }
      }
    }
      break; // StepPoints

    case Mu2eG4PrimaryType::SimParticleLeaves: {
      evt.getByLabel(primaryTag_, res);
      if(!res.isValid()) {
        throw cet::exception("BADINPUT")
          <<"Mu2eG4Inputs::inputSimParticles(): could not get SimParticleCollection "
          <<primaryTag_
          <<std::endl;
      }

    }
      break; // SimParticles
    }

    return res;
  }

  //================================================================
}
