// "Raw" pixel clusters are groups of adjacent hits.
//
//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALRawCluster_hh
#define RecoDataProducts_ExtMonFNALRawCluster_hh

#include <ostream>

#include "canvas/Persistency/Common/PtrVector.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"

namespace mu2e {

  class ExtMonFNALRawCluster {
  public:

    typedef art::PtrVector<ExtMonFNALRawHit> Hits;

    const Hits&  hits() const { return hits_; }

    explicit ExtMonFNALRawCluster(const  Hits&  hits) :  hits_(hits) {}

    // Default constructor for ROOT persistency
    ExtMonFNALRawCluster() : hits_() {}

  private:
    Hits hits_;
  };

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALRawCluster& c);

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALRawCluster_hh */
