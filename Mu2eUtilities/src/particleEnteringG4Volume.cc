#include "Mu2eUtilities/inc/particleEnteringG4Volume.hh"

art::Ptr<mu2e::SimParticle> mu2e::particleEnteringG4Volume(const StepPointMC& step) {
    art::Ptr<SimParticle> p = step.simParticle();

    // No checks on p here.  StepPoint partricles must be resolvable.
    // If p is null or not available a bug need be fixed upstream.
    // Similarly, if a parent isNonnull but not resolvable that's
    // a bug elsewhere.

    while((step.volumeId() == p->startVolumeIndex())&&(p->parent().isNonnull())) {
      p = p->parent();
    }

    return p;
}
