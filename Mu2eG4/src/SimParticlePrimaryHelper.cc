// Andrei Gaponenko, 2013

#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StageParticle.hh"

#include "cetlib_except/exception.h"

namespace mu2e {

  //================================================================
  SimParticlePrimaryHelper::SimParticlePrimaryHelper(const art::ProductID& simProdID,
                                                     const art::EDProductGetter* sim_prod_getter):
    simProdID_(simProdID),
    simProductGetter_(sim_prod_getter)
  {}

  //================================================================
  art::Ptr<GenParticle> SimParticlePrimaryHelper::genParticlePtr(int g4TrkID) const {
    art::Ptr<GenParticle> res;
    auto orig = getEntry(g4TrkID);
    auto pgen = std::get_if<art::Ptr<GenParticle> >(&orig);
    if(pgen) {
      res = *pgen;
    }
    return res;
  }

  //================================================================
  art::Ptr<SimParticle>  SimParticlePrimaryHelper::simParticlePrimaryPtr(int g4TrkID) const {
    auto orig = getEntry(g4TrkID);

    if(std::holds_alternative<art::Ptr<GenParticle>>(orig)) {
      return art::Ptr<SimParticle>();
    }
    else {

      SimParticle::key_type id;
      if(std::holds_alternative<const SimParticle*>(orig)) {
        std::cout<<"simParticlePrimaryPtr from SimParticle*"<<std::endl;
        id = std::get<const SimParticle*>(orig)->id();
      }
      else if(std::holds_alternative<const  StepPointMC*>(orig)) {
        std::cout<<"simParticlePrimaryPtr from StepPointMC*"<<std::endl;
        id = std::get<const StepPointMC*>(orig)->simParticle()->id();
      }
      else if(std::holds_alternative<const StageParticle*>(orig)) {
        std::cout<<"simParticlePrimaryPtr from StageParticle*"<<std::endl;
        id = std::get<const StageParticle*>(orig)->parent()->id();
      }

      try { id.ensure_valid(); }
      catch(cet::exception& e) {
        throw cet::exception("ImplementationError", "SimParticlePrimaryHelper::simParticlePrimaryPtr(int) missed a variant alternative", e);
      }

      return art::Ptr<SimParticle>(simProdID_, id.asUint(), simProductGetter_);
    }
  }

  //================================================================
  SimParticlePrimaryHelper::InputParticle SimParticlePrimaryHelper::getEntry(int g4TrkID) const {
    return entries_.at(g4TrkID-1);
  }

  //================================================================

}
