#ifndef RecoDataProducts_HelixSeed_hh
#define RecoDataProducts_HelixSeed_hh
//
// Seed for a track using a robust helix representation
// Original author Dave Brown (LBNL) 19 Aug 2016
//

// Mu2e includes
#include "BTrk/TrkBase/TrkT0.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  class CaloCluster;
  class TimeCluster;

  struct HelixSeed {

    TrkT0 const& t0() const { return _t0; }
    ComboHitCollection const& hits() const { return _hhits; }
    RobustHelix const& helix() const { return _helix; }
    TrkFitFlag const& status() const { return _status; }
    art::Ptr<CaloCluster> const& caloCluster() const { return _timeCluster->caloCluster(); }
    art::Ptr<TimeCluster> const& timeCluster() const { return _timeCluster; }

    TrkT0	             _t0;	      // t0 for this helix
    ComboHitCollection       _hhits;	      // hits potentially used for this helix
    RobustHelix              _helix;	     // robust helix created from these hits
    TrkFitFlag	             _status;      // status of processes used to create this seed
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
  };
   typedef std::vector<mu2e::HelixSeed> HelixSeedCollection;
} // namespace mu2e

#endif /* RecoDataProducts_HelixSeed_hh */
