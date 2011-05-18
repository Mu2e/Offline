#ifndef ToyDP_StrawCluster_hh
#define ToyDP_StrawCluster_hh
//
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.3 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
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

#endif /* ToyDP_StrawCluster_hh */
