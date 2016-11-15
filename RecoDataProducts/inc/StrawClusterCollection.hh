#ifndef RecoDataProducts_StrawClusterCollection_hh
#define RecoDataProducts_StrawClusterCollection_hh

//
// Define a type for a collection of StrawCluster objects.
//
// $Id: StrawClusterCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Hans Wenzel
//

#include <vector>

#include "RecoDataProducts/inc/StrawCluster.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawCluster> StrawClusterCollection;
}

#endif /* RecoDataProducts_StrawClusterCollection_hh */
