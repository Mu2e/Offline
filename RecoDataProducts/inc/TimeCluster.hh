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
    std::vector<StrawHitIndex> const& hits() const { return _strawHitIdxs;}
    TrkT0 const& t0() const { return _t0; }
    CLHEP::Hep3Vector const& position() const { return _pos; }
    art::Ptr<CaloCluster> const& caloCluster() const { return _caloCluster; }

    std::vector<StrawHitIndex>    _strawHitIdxs; // associated straw hits: can be empty
    TrkT0		     _t0; // t0 and associated error
    CLHEP::Hep3Vector        _pos; // position of the time cluster   
    art::Ptr<CaloCluster>    _caloCluster; // associated calorimeter cluster: can be null
  };
   typedef std::vector<mu2e::TimeCluster> TimeClusterCollection;

} // namespace mu2e

#endif /* RecoDataProducts_TimeCluster_hh */
