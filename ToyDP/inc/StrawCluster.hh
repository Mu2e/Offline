#ifndef StrawCluster_H
#define StrawCluster_H
// 
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.1 2011/01/11 20:42:08 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/01/11 20:42:08 $
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

#endif
