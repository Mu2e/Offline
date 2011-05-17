#ifndef ToyDP_StrawClusterCollection_hh
#define ToyDP_StrawClusterCollection_hh

//
// Define a type for a collection of StrawCluster objects.
//
// $Id: StrawClusterCollection.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Hans Wenzel
//

#include <vector>

#include "ToyDP/inc/StrawCluster.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawCluster> StrawClusterCollection;
}

#endif /* ToyDP_StrawClusterCollection_hh */
