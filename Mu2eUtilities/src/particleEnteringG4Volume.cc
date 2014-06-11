#include "Mu2eUtilities/inc/particleEnteringG4Volume.hh"

art::Ptr<mu2e::SimParticle> mu2e::particleEnteringG4Volume(const StepPointMC& step) {
    art::Ptr<SimParticle> p = step.simParticle();
    // In Mu2e we save all ancestors of a particle, so navigating up is OK.
    while(p.isNonnull() && (step.volumeId() == p->startVolumeIndex())) {
      p = p->parent();
    }
    return p;
}
