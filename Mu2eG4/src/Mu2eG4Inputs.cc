#include "Offline/Mu2eG4/inc/Mu2eG4Inputs.hh"

#include <optional>
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"

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
  Mu2eG4Inputs::InputSimsInfo
  Mu2eG4Inputs::inputSimParticles(const art::Event& evt) const {
    InputSimsInfo res;

    switch(primaryType_.id()) {
    default: throw cet::exception("CONFIG")
        << "Error: Mu2eG4Inputs::inputSimParticles(): unknown Mu2eG4 primaryType id = "
        <<primaryType_.id()
        <<std::endl;

    case Mu2eG4PrimaryType::GenParticles:
      break; // No input sims for GenParticle jobs

    case Mu2eG4PrimaryType::StepPoints: {
      art::Ptr<SimParticle> simPtr;
      std::optional<art::ProductID> simid;

      auto const h = evt.getValidHandle<StepPointMCCollection>(primaryTag_);
      for(const auto& hit : *h) {
        if(!simid) {
          simid = hit.simParticle().id();
          simPtr = hit.simParticle();
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
        res.reseat(simPtr.parentAs<cet::map_vector<SimParticle>>(), *simid);
       }
    }
      break; // StepPoints

    case Mu2eG4PrimaryType::SimParticleLeaves: {
      art::Handle<SimParticleCollection> sph;
      evt.getByLabel(primaryTag_, sph);
      if(!sph.isValid()) {
        throw cet::exception("BADINPUT")
          <<"Mu2eG4Inputs::inputSimParticles(): could not get SimParticleCollection "
          <<primaryTag_
          <<std::endl;
      }
      res.reseat(*sph, sph.id() );

    }
      break; // SimParticles

    case Mu2eG4PrimaryType::StageParticles: {
      art::Ptr<SimParticle> simPtr;
      std::optional<art::ProductID> simid;

      auto const h = evt.getValidHandle<StageParticleCollection>(primaryTag_);
      for(const auto& s : *h) {
        if(!simid) {
          simid = s.parent().id();
          simPtr = s.parent();
        }
        else {
          if(*simid != s.parent().id()) {
            throw cet::exception("BADINPUT")
              <<"Mu2eG4Inputs::inputSimParticles(): inconsistent SimParticleCollection product ID "
              <<"in input StageParticleCollection "
              <<primaryTag_
              <<std::endl;
          }
        }
      }

      if(simid) { // nothing to retrieve if primary inputs are emtpy
        res.reseat(simPtr.parentAs<cet::map_vector<SimParticle>>(), *simid);
      }
    }
      break; // StageParticles

    }

    return res;
  }

  //================================================================
}
