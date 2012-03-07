#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/FindOne.h"
#include "art/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/GenParticleSPMHistory.hh"

#include "cetlib/exception.h"

#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
//#define AGDEBUG(stuff)

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
      if(map_.empty()) { // need to load the associations

        typedef std::vector<art::Handle<GenParticleSPMHistory> > Handles;
        Handles results;
        evt_->getManyByType(results);

        for(Handles::const_iterator h = results.begin(); h != results.end(); ++h) {
          AGDEBUG("In loop over handles");
          // could filter out entries here
          for(GenParticleSPMHistory::assn_iterator i = (*h)->begin(); i != (*h)->end(); ++i) {
            AGDEBUG("In loop over assns");
            if(!map_.insert(std::make_pair(i->first, i->second)).second) {
              throw cet::exception("BADINPUTS")<<"Non-unique GenParticle ptr in GenParticleSPMHistory";
            }
          }
        }
      }

      MapType::const_iterator hit = map_.find(gp);

      return (hit == map_.end()) ? nullPtr_ : hit->second->simParticle();

    } // else  (!hasParent)

  } // parent()

} // namespace mu2e
