#ifndef RecoDataProducts_HelixSeed_hh
#define RecoDataProducts_HelixSeed_hh
//
// Seed for a track using a robust helix representation
// Original author Dave Brown (LBNL) 19 Aug 2016
//

// Mu2e includes
#include "BTrk/TrkBase/TrkT0.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/HelixRecoDir.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  class CaloCluster;
  class TimeCluster;
  using FitDir = TrkFitDirection::FitDirection;
  struct HelixSeed {

    TrkT0                 const& t0()          const { return _t0; }
    ComboHitCollection    const& hits()        const { return _hhits; }
    RobustHelix           const& helix()       const { return _helix; }
    TrkFitFlag            const& status()      const { return _status; }
    HelixRecoDir          const& recoDir()     const { return _recoDir; }
    FitDir                const& propDir()     const { return _propDir; }
    float                 const& eDepAvg()     const { return _eDepAvg; }
    art::Ptr<CaloCluster> const& caloCluster() const { return _timeCluster->caloCluster(); }
    art::Ptr<TimeCluster> const& timeCluster() const { return _timeCluster; }

    TrkT0                    _t0;          // t0 for this helix
    ComboHitCollection       _hhits;       // hits potentially used for this helix
    RobustHelix              _helix;       // robust helix created from these hits
    TrkFitFlag               _status;      // status of processes used to create this seed
    HelixRecoDir             _recoDir;     // sign of the longitudinal velocity (z-axis) derived from a T vs Z linear fit
    FitDir                    _propDir = FitDir::unknown; // direction of propagation (upstream, downstream, or ambiguous)
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
    float                    _eDepAvg =0;  // average energy deposition from helix hits
  };
   typedef std::vector<mu2e::HelixSeed> HelixSeedCollection;
   typedef art::Ptr<mu2e::HelixSeed> HelixSeedPtr;
   typedef std::vector<mu2e::HelixSeedPtr> HelixSeedPtrCollection;
} // namespace mu2e

#endif /* RecoDataProducts_HelixSeed_hh */
