#ifndef Mu2eUtilities_SortedStepPoints_hh
#define Mu2eUtilities_SortedStepPoints_hh
//
//
// $Id: SortedStepPoints.hh,v 1.5 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/ThreeVector.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include <vector>

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
