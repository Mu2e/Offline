#ifndef RecoDataProducts_StrawCluster_hh
#define RecoDataProducts_StrawCluster_hh
//
// First version of a Cluster.
//
// $Id: StrawCluster.hh,v 1.3 2011/06/02 22:50:54 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/02 22:50:54 $
//
// Original author Hans Wenzel
//

// C++ includes
#include <vector>



#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveTransients.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "Mu2eUtilities/inc/LineSegmentPCA.hh"
// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
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
    CLHEP::Hep3Vector X(art::Event const& event) const;
    double Energy(art::Event const& event) const;
    double Halflength(art::Event const& event) const;
    double averageT(art::Event const& event) const;
    double averagedT(art::Event const & event) const;
    DeviceId did(art::Event const & event) const;
    SectorId secid(art::Event const & event) const;
    CLHEP::Hep3Vector dirX(art::Event const& event) const; 
    LineSegmentPCA linesegment(art::Event const& event) const;
    //    StrawCluster& add(art::ProductID const & CollId, CaloHit const & hit);
  private:
    //   const Tracker& tracker;
    std::vector<DPIndex> _StrawHitIndices;
  };
} // namespace mu2e

#endif /* RecoDataProducts_StrawCluster_hh */
