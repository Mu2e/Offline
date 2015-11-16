#ifndef HitMakers_StrawClusterUtilities_hh
#define HitMakers_StrawClusterUtilities_hh
//
// First version of a Cluster.
//
// Original author Hans Wenzel
//

// C++ includes
#include <vector>
#include <map>
// Framework includes:
#include "art/Framework/Principal/Event.h"
// Mu2e includes:
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "GeneralUtilities/inc/LineSegmentPCA.hh"
#include "DataProducts/inc/DeviceId.hh"
#include "DataProducts/inc/SectorId.hh"

// CLHEP includes:
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
namespace art {
  class ProductID;
}

namespace mu2e {
  class StrawClusterUtilities{

  public:

    StrawClusterUtilities() {}

    // Accessors
    CLHEP::Hep3Vector midX(StrawCluster const& cluster,art::Event const& event) const;
    CLHEP::Hep3Vector dirX(StrawCluster const& cluster,art::Event const& event) const;
    CLHEP::Hep3Vector dTX(StrawCluster const& cluster,art::Event const& event) const;
    double Energy(StrawCluster const& cluster,art::Event const& event) const;
    double Halflength(StrawCluster const& cluster,art::Event const& event) const;
    double averageT(StrawCluster const& cluster,art::Event const& event) const;
    double averagedT(StrawCluster const& cluster,art::Event const & event) const;
    DeviceId did(StrawCluster const& cluster,art::Event const & event) const;
    SectorId secid(StrawCluster const& cluster,art::Event const & event) const;
    int Station(StrawCluster const& cluster,art::Event const & event) const;

    LineSegmentPCA linesegment(StrawCluster const& cluster,art::Event const& event) const;
    std::multimap<int,StrawCluster> clusterbydid(StrawClusterCollection const& clusters,art::Event const& event) const;
    std::multimap<int,StrawCluster> clusterbystation(StrawClusterCollection const& clusters,art::Event const& event) const;
  };
} // namespace mu2e

#endif /* HitMakers_StrawClusterUtilities_hh */
