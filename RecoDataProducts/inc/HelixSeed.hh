#ifndef RecoDataProducts_HelixSeed_hh
#define RecoDataProducts_HelixSeed_hh
//
// Seed for a track using a robust helix representation
// Original author Dave Brown (LBNL) 19 Aug 2016
//

// Mu2e includes
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"

namespace mu2e {

  struct HelixSeed {
    TimeCluster        _timeCluster; // timing and hits associated with this helix
    RobustHelix        _helix;	     // robust helix created from these hits
    TrkFitFlag	       _status;      // status of processes used to create this seed
  };

} // namespace mu2e

#endif /* RecoDataProducts_HelixSeed_hh */
