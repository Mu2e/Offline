#ifndef growCluster_hh
#define growCluster_hh
//
// Free function to grow a cluster.
//
// $Id: growCluster.hh,v 1.1 2009/10/28 14:14:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/28 14:14:13 $
//
//
// Consider all straws in the cluster in the range,
// int i=startIndex, i<cluster.size()
//
// Loop over the unused hits.  If there are any unused hits that
// are nearest neighbours of a straw in the cluster, add that 
// hit to the cluster.
//
// Arguments:
// 1 - reference to cluster to be grown.
// 2 - index within the cluster of the first hit to be considered.
//     ( usually the first hit added in the previous call ).
// 3 - index into the hit list of the first hit to be considered.
//     All hits before this are known to be used.  Hits
//     after may be used or unused.
// 4 - the container of this
// 5 - A used hit map: 0= unused, 1=used.
//     This must be maintained by this code.
// 6 - The LTracker geometry.
//

// C++ includes
#include <vector>

// Framework includes
#include "DataFormats/Common/interface/Handle.h"

// Mu2e includes
#include "ToyDP/inc/CrudeStrawHitPData.hh"

namespace mu2e {

  class ProtoStrawCluster;
  class LTracker;
  
  int growCluster ( ProtoStrawCluster&              cluster,
		    int                             startCluster,
		    int                             startHit,
		    edm::Handle<CrudeStrawHitPData> pdataHandle,
		    std::vector<int>&               used,
		    LTracker const&                 ltracker
		    );
}

#endif
