#ifndef ToyDP_StrawClusterCollection_hh
#define ToyDP_StrawClusterCollection_hh

//
// Define a type for a collection of StrawCluster objects.
//
// $Id: StrawClusterCollection.hh,v 1.1 2011/01/11 20:42:28 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/11 20:42:28 $
//
// Original author Hans Wenzel
//

#include <vector>

#include "ToyDP/inc/StrawCluster.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawCluster> StrawClusterCollection;
}

#endif
