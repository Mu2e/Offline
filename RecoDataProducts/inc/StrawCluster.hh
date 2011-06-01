#ifndef RecoDataProducts_StrawCluster_hh
#define RecoDataProducts_StrawCluster_hh
//
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.2 2011/06/01 21:12:23 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/01 21:12:23 $
//
// Original author Hans Wenzel
//

// C++ includes
#include <vector>

// Mu2e includes
#include "DataProducts/inc/DPIndex.hh"

namespace art {
  class ProductID;
}

namespace mu2e {
  class StrawCluster{

  public:

    StrawCluster() {}

    StrawCluster(std::vector<DPIndex> & hitIndices);

    // Accessors
    std::vector<DPIndex> const & StrawHitIndices() const { return _StrawHitIndices; }

    //    StrawCluster& add(art::ProductID const & CollId, CaloHit const & hit);
  private:

    std::vector<DPIndex> _StrawHitIndices;
  };
} // namespace mu2e

#endif /* RecoDataProducts_StrawCluster_hh */
