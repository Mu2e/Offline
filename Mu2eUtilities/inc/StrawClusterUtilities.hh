#ifndef HitMakers_StrawClusterUtilities_hh
#define HitMakers_StrawClusterUtilities_hh
//
// First version of a Cluster.
//
// $Id: StrawClusterUtilities.hh,v 1.2 2011/06/24 11:44:59 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/24 11:44:59 $
//
// Original author Hans Wenzel
//

// C++ includes
#include <vector>
#include <map>
// Framework includes:
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
// Mu2e includes:
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "Mu2eUtilities/inc/LineSegmentPCA.hh"
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
