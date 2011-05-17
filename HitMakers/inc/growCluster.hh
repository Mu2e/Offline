#ifndef HitMakers_growCluster_hh
#define HitMakers_growCluster_hh
//
// Free function to grow a cluster.
//
// $Id: growCluster.hh,v 1.4 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
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
#include "art/Persistency/Common/Handle.h"

// Mu2e includes
#include "ToyDP/inc/CrudeStrawHitPData.hh"

namespace mu2e {

  class ProtoStrawCluster;
  class LTracker;
  
  int growCluster ( ProtoStrawCluster&              cluster,
                    int                             startCluster,
                    int                             startHit,
                    art::Handle<CrudeStrawHitPData> pdataHandle,
                    std::vector<int>&               used,
                    LTracker const&                 ltracker
                    );
}

#endif /* HitMakers_growCluster_hh */
