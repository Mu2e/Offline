#ifndef RecoDataProducts_StrawCluster_hh
#define RecoDataProducts_StrawCluster_hh
//
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.5 2011/06/06 22:51:25 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/06 22:51:25 $
//
// Original author Hans Wenzel
//

// C++ includes
#include <vector>
// Mu2e includes:
#include "DataProducts/inc/DPIndex.hh"
namespace art {
  class ProductID;
}

namespace mu2e {
  class StrawCluster{

  public:

    StrawCluster() {}

    StrawCluster(std::vector<DPIndex> & hitIndices)  {   // Geometry info for the TTracker.
    // Get a reference to one of the T trackers.
    // Throw exception if not successful.
    _StrawHitIndices=hitIndices;
  };

    // Accessors
    std::vector<DPIndex> const & StrawHitIndices() const { return _StrawHitIndices; }
    //    StrawCluster& add(art::ProductID const & CollId, CaloHit const & hit);
  private:
    //   const Tracker& tracker;
    std::vector<DPIndex> _StrawHitIndices;
  };
} // namespace mu2e

#endif /* RecoDataProducts_StrawCluster_hh */
