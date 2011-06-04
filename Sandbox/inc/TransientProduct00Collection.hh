#ifndef Sandbox_TransientProduct00Collection_hh
#define Sandbox_TransientProduct00Collection_hh

//
// Define a type for a collection of TransientProduct objects.
// Used for tests of making transient-only data products.
//
// $Id: TransientProduct00Collection.hh,v 1.1 2011/06/04 18:00:36 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/04 18:00:36 $
//
// Original author Rob Kutschke
//

#include <vector>

#include "Sandbox/inc/TransientProduct00.hh"

namespace mu2e {
   typedef std::vector<mu2e::TransientProduct00> TransientProduct00Collection;
}

#endif /* Sandbox_TransientProduct00Collection_hh */
