#include <algorithm>

#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
// Framework includes
#include "canvas/Persistency/Common/Ptr.h"
// cetlib includes
#include "cetlib_except/exception.h"

#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"

namespace mu2e {

  //==========================================================================
  double SimParticleGetTau::calculate( const art::Ptr<SimParticle>& p,
                                       const VspMC& hitColls,
                                       const std::vector<int>& decayOffCodes,
                                       const PhysicsParams& gc ){

    double tau = p->endProperTime() / gc.getParticleLifetime(p->pdgId());
    tau += getMultiStageTau( p, hitColls, decayOffCodes, gc );
    return tau;

  }

  //==========================================================================
  double SimParticleGetTau::calculate( const StepPointMC& sp,
                                       const VspMC& hitColls,
                                       const std::vector<int>& decayOffCodes,
                                       const PhysicsParams& gc ){

    double tau = sp.properTime() / gc.getParticleLifetime(sp.simParticle()->pdgId());
    tau += getMultiStageTau( sp.simParticle(), hitColls, decayOffCodes, gc );
    return tau;

  }

  //==========================================================================
  double SimParticleGetTau::calculate( const art::Ptr<SimParticle>& p,
                                       const std::vector<int>& decayOffCodes,
                                       const PhysicsParams& gc ){
    double tau = 0.;
    const auto pdgId = p->pdgId();
    // check if this code was turned off
    if (std::find(decayOffCodes.begin(), decayOffCodes.end(), pdgId) != decayOffCodes.end())
      tau = p->endProperTime() / gc.getParticleLifetime(pdgId);
     // continue through the history
    if(p->parent().isNonnull()) tau += calculate(p->parent(), decayOffCodes, gc);
    return tau;
  }

  //==========================================================================
  double SimParticleGetTau::getMultiStageTau( const art::Ptr<SimParticle>& p,
                                              const VspMC& hitColls,
                                              const std::vector<int>& decayOffCodes,
                                              const PhysicsParams& gc ) {

    double tau(0.);

    // The mu2ePrimary code means that G4 track was created by our
    // PrimaryGeneratorAction.  If the particle is "mu2ePrimary" but
    // still has a parent, it is a continuation of a particle from the
    // previous simulation stage, and their proper times should be
    // combined.

    art::Ptr<SimParticle> part (p);
    while(part->parent().isNonnull()) {

      if((part->creationCode() == ProcessCode::mu2ePrimary)) {

        // The current particle is a continuation from the previous
        // stage, not a physically different particle.  We need to
        // find its record in the StepPointMC collections to get the
        // correct proper time.
        part = part->parent();

        // Find matches in hit collections
        unsigned counter (0);
        const StepPointMC* spMC(nullptr);
        for ( const auto& hitColl : hitColls ) {
          std::for_each( hitColl.begin(), hitColl.end(),
                         [&](const StepPointMC& sp){
                           if ( sp.simParticle().key() == part.key() ) {
                             spMC = &sp;
                             counter++;
                           }
                         } );

        }

        if      ( counter == 0 ) throw cet::exception("StepPointMC") << " Non-existent StepPointMC-SimParticle assignment! " ;
        else if ( counter  > 1 ) throw cet::exception("StepPointMC") << " Ambiguous StepPointMC-SimParticle assignment! " ;
        else  tau += spMC->properTime() / gc.getParticleLifetime(part->pdgId());

      }
      else {
        // The current particle was produced by a G4 physics process.
        // See if proper time of its ancestor should be included.
        part = part->parent();
        if ( std::binary_search( decayOffCodes.begin(), decayOffCodes.end(), int(part->pdgId()) ) ) {
          tau += part->endProperTime() / gc.getParticleLifetime(part->pdgId());
        }
      }
    } // loop up to the primary

    return tau;
  }

} // namespace mu2e
