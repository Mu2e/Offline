#ifndef ToyDP_StrawCluster_hh
#define ToyDP_StrawCluster_hh
// 
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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
