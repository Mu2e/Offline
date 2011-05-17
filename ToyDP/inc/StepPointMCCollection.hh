#ifndef ToyDP_StepPointMCCollection_hh
#define ToyDP_StepPointMCCollection_hh

//
// Define a type for a collection of StepPointMC objects.
//
// $Id: StepPointMCCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/StepPointMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::StepPointMC> StepPointMCCollection;
}

#endif /* ToyDP_StepPointMCCollection_hh */
