#ifndef RecoDataProducts_HelixSeed_hh
#define RecoDataProducts_HelixSeed_hh
//
// Seed for a track using a robust helix representation
// Original author Dave Brown (LBNL) 19 Aug 2016
//

// Mu2e includes
#include "BTrk/TrkBase/TrkT0.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "art/Persistency/Common/Ptr.h"

namespace mu2e {
  class CaloCluster;

  struct HelixSeed {

    TrkT0 const& t0() const { return _t0; }
    HelixHitCollection const& hits() const { return _hhits; }
    RobustHelix const& helix() const { return _helix; }
    TrkFitFlag const& status() const { return _status; }
    art::Ptr<CaloCluster> const& caloCluster() const { return _caloCluster; }

    TrkT0	       _t0;	      // t0 for this helix
    HelixHitCollection _hhits;	      // hits potentially used for this helix
    RobustHelix        _helix;	     // robust helix created from these hits
    TrkFitFlag	       _status;      // status of processes used to create this seed
    art::Ptr<CaloCluster>    _caloCluster; // associated calorimeter cluster: can be null
  };

} // namespace mu2e

#endif /* RecoDataProducts_HelixSeed_hh */
