#include <memory>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include "MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "MCDataProducts/inc/GenSimParticleLink.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib_except/exception.h"

#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

namespace mu2e { class GenParticle; }

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  art::Ptr<SimParticle> SimParticleParentGetter::nullPtr_;

  //================================================================
  SimParticleParentGetter::SimParticleParentGetter(const art::Event& evt)
    : evt_(&evt)
  {}

  //================================================================
  const art::Ptr<SimParticle>& SimParticleParentGetter::parent(const art::Ptr<SimParticle>& particle) const {
    if(particle->hasParent()) {
      return particle->parent();
    }
    else {
      const art::Ptr<GenParticle>& gp = particle->genParticle();
      if(gp.isNull()) {
        throw cet::exception("BADINPUTS")<<"SimParticleParentGetter::parent(): missing info: particle does not have parent and not from generator";
      }

      // try to find associated hits
      if(stepPointMap_.empty()) { // need to load the associations

        typedef std::vector<art::Handle<GenParticleSPMHistory> > StepPointHandles;
        StepPointHandles stepPointResults = evt_->getMany<GenParticleSPMHistory>();

	if (stepPointResults.empty() ) { //no associations with StepPoints. Trying with simparticle end points

	  if (simParticleMap_.empty()) { //need to load the associations

	    typedef std::vector<art::Handle<GenSimParticleLink> > SimParticleHandles;
	    SimParticleHandles simParticleResults = evt_->getMany<GenSimParticleLink>();

	    for(SimParticleHandles::const_iterator h = simParticleResults.begin(); h != simParticleResults.end(); ++h) {
	      AGDEBUG("In loop over simParticleHandles");
	      // could filter out entries here
	      for(GenSimParticleLink::const_iterator i = (*h)->begin(); i != (*h)->end(); ++i) {
		AGDEBUG("In loop over sim particle assns");
		if(!simParticleMap_.insert(std::make_pair(i->first, i->second)).second) {
		  throw cet::exception("BADINPUTS")<<"Non-unique GenParticle ptr in GenSimParticleLink";
		}
	      }
	    }
	  }

	  SimParticleMapType::const_iterator gensimpair = simParticleMap_.find(gp);

	  return (gensimpair == simParticleMap_.end()) ? nullPtr_ : gensimpair->second;

	}

        for(StepPointHandles::const_iterator h = stepPointResults.begin(); h != stepPointResults.end(); ++h) {
          AGDEBUG("In loop over stepPointHandles");
          // could filter out entries here
          for(GenParticleSPMHistory::const_iterator i = (*h)->begin(); i != (*h)->end(); ++i) {
            AGDEBUG("In loop over step point assns");
            if(!stepPointMap_.insert(std::make_pair(i->first, i->second)).second) {
              throw cet::exception("BADINPUTS")<<"Non-unique GenParticle ptr in GenParticleSPMHistory";
            }
          }
        }
      }

      StepPointMapType::const_iterator hit = stepPointMap_.find(gp);

      return (hit == stepPointMap_.end()) ? nullPtr_ : hit->second->simParticle();

    } // else  (!hasParent)

  } // parent()

} // namespace mu2e
