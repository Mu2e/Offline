#ifndef Mu2eUtilities_SortedStepPoints_hh
#define Mu2eUtilities_SortedStepPoints_hh
//
//
// $Id: SortedStepPoints.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/StepPointMCCollection.hh"

#include "ToyDP/inc/SimParticleCollection.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

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

    //MapVectorKey trackId_;
    SimParticleCollection::key_type trackId_;

    std::vector<StepPointMC const*> steps_;

    StepPointMC const* midByZ_;


  };

} // namespace mu2e

#endif /* Mu2eUtilities_SortedStepPoints_hh */
