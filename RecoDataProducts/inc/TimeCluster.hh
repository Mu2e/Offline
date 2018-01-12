#ifndef RecoDataProducts_TimeCluster_hh
#define RecoDataProducts_TimeCluster_hh
//
// Defin a 
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>
// art includes
#include "canvas/Persistency/Common/Ptr.h"
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "CLHEP/Vector/ThreeVector.h"
// BTrk includes
#include "BTrk/TrkBase/TrkT0.hh"
// c++
#include <vector>
namespace mu2e {
  class CaloCluster;
 
  struct TimeCluster{
    int                               nhits      () const { return _strawHitIdxs.size(); }
    const std::vector<StrawHitIndex>& hits       () const { return _strawHitIdxs; }
    const TrkT0&                      t0         () const { return _t0; }
    const CLHEP::Hep3Vector&          position   () const { return _pos; }
    const art::Ptr<CaloCluster>&      caloCluster() const { return _caloCluster; }
    int                               caloFastIdx() const { return _caloFastIdx; }

    std::vector<StrawHitIndex> _strawHitIdxs; // associated straw hits: can be empty
    TrkT0		       _t0;           // t0 and associated error
    CLHEP::Hep3Vector          _pos;          // position of the time cluster   
    art::Ptr<CaloCluster>      _caloCluster;  // associated calorimeter cluster: can be null
    int                        _caloFastIdx;
  };
   typedef std::vector<mu2e::TimeCluster> TimeClusterCollection;

} // namespace mu2e

#endif /* RecoDataProducts_TimeCluster_hh */
