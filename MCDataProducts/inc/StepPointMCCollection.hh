#ifndef MCDataProducts_StepPointMCCollection_hh
#define MCDataProducts_StepPointMCCollection_hh

//
// Define a type for a collection of StepPointMC objects.
//
// $Id: StepPointMCCollection.hh,v 1.1 2011/05/24 17:16:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:44 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::StepPointMC> StepPointMCCollection;
}

#endif /* MCDataProducts_StepPointMCCollection_hh */
