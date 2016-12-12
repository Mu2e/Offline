#ifndef MCDataProducts_PtrStepPointMCVectorCollection_hh
#define MCDataProducts_PtrStepPointMCVectorCollection_hh
//
// Define a type for a collection of art::Ptr<StepPointMC> objects.
//
// Original author Ivan Logashenko
//

#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "canvas/Persistency/Common/Ptr.h"

#include <vector>

namespace mu2e {
  typedef std::vector<PtrStepPointMCVector> PtrStepPointMCVectorCollection;
}

#endif /* MCDataProducts_PtrStepPointMCVectorCollection_hh */
