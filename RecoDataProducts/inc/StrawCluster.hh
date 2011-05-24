#ifndef RecoDataProducts_StrawCluster_hh
#define RecoDataProducts_StrawCluster_hh
//
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Hans Wenzel
//

// C++ includes
#include <vector>

// Mu2e includes
#include "TrackerGeom/inc/StrawId.hh"

namespace mu2e {

      typedef std::vector<mu2e::StrawId> StrawCluster;
} // namespace mu2e

#endif /* RecoDataProducts_StrawCluster_hh */
