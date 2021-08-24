#ifndef Mu2eUtilities_SortedStepPoints_hh
#define Mu2eUtilities_SortedStepPoints_hh
//
// Original author Rob Kutschke
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  class StepPointMC;

  class SortedStepPoints{

  public:
    SortedStepPoints ( SimParticleCollection::key_type trackId,
                       StepPointMCCollection const&    allSteps );

    // Accept compiler supplied d'tor, copy c'tor and assignment operator.

    // Accessors

    std::vector<StepPointMC const*> const& steps() const { return steps_;}

    StepPointMC const& middleByZ() const { return *midByZ_; }

  private:

    SimParticleCollection::key_type trackId_;

    std::vector<StepPointMC const*> steps_;

    StepPointMC const* midByZ_;

  };

} // namespace mu2e

#endif /* Mu2eUtilities_SortedStepPoints_hh */
