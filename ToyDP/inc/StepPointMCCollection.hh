#ifndef ToyDP_StepPointMCCollection_hh
#define ToyDP_StepPointMCCollection_hh

//
// Define a type for a collection of StepPointMC objects.
//
// $Id: StepPointMCCollection.hh,v 1.1 2009/10/06 23:19:59 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/06 23:19:59 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "ToyDP/inc/StepPointMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::StepPointMC> StepPointMCCollection;
}

#endif
