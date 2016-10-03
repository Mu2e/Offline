#ifndef MCDataProducts_PtrStepPointMCVector_hh
#define MCDataProducts_PtrStepPointMCVector_hh
//
// Define a type for a collection of art::Ptr<StepPointMC> objects.
//
// Original author Ivan Logashenko
//

#include "MCDataProducts/inc/StepPointMC.hh"
#include "canvas/Persistency/Common/Ptr.h"

#include <vector>

namespace mu2e {
  typedef std::vector<art::Ptr<StepPointMC> > PtrStepPointMCVector;
}

#endif /* MCDataProducts_PtrStepPointMCVector_hh */
