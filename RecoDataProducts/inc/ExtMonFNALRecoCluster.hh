// Reconstructed pixel clusters are calibrated objects with alignment
// corrections applied.  The coordinates of a cluster are given in the
// respective SensorStack system.
//
//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALRecoCluster_hh
#define RecoDataProducts_ExtMonFNALRecoCluster_hh

#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "canvas/Persistency/Common/PtrVector.h"

#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"

namespace mu2e {

  class ExtMonFNALRecoCluster {
  public:

    // Int the spectrometer numbering (global)
    unsigned int plane() const { return plane_; }

    // in SensorStack coordinates of the cluster plane
    const CLHEP::Hep3Vector& position() const { return position_; }

    int xWidth() const { return xWidth_; }
    int yWidth() const { return yWidth_; }

    int clock() const { return clock_; }

    // The corresponding raw cluster
    const art::Ptr<ExtMonFNALRawCluster>& raw() const { return raw_; }

    explicit ExtMonFNALRecoCluster(const art::Ptr<ExtMonFNALRawCluster>& raw,
                                   unsigned int plane,
                                   const CLHEP::Hep3Vector& position,
                                   int xWidth,
                                   int yWidth,
                                   int clock);

    // Default constructor for ROOT persistency
    ExtMonFNALRecoCluster() : raw_(), position_(), plane_(), xWidth_(), yWidth_(), clock_() {}

  private:
    art::Ptr<ExtMonFNALRawCluster> raw_;
    CLHEP::Hep3Vector position_;
    unsigned int plane_;
    int xWidth_;
    int yWidth_;
    int clock_;
  };

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALRecoCluster& c);

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALRecoCluster_hh */
