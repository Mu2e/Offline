#ifndef RecoDataProducts_TimeCluster_hh
#define RecoDataProducts_TimeCluster_hh
//
// Defin a 
//
//
// Original author G. Tassielli
//

// C++ includes
#include <vector>
// art includes
#include "canvas/Persistency/Common/Ptr.h"
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "DataProducts/inc/XYZVec.hh"
// BTrk includes
#include "BTrk/TrkBase/TrkT0.hh"
// c++
#include <vector>
namespace mu2e {
  class CaloCluster;
 
  struct TimeCluster{
    TimeCluster() : _nsh(0)  {}
    size_t                               nhits      () const { return _strawHitIdxs.size(); }
    uint16_t	  nStrawHits() const { return _nsh; }
    const std::vector<StrawHitIndex>& hits       () const { return _strawHitIdxs; }
    const TrkT0&                      t0         () const { return _t0; }
    const XYZVec&          position   () const { return _pos; }
    const art::Ptr<CaloCluster>&      caloCluster() const { return _caloCluster; }
    bool hasCaloCluster() const { return _caloCluster.isNonnull(); }
    std::vector<StrawHitIndex> _strawHitIdxs; // associated straw hits: can be empty
    TrkT0		       _t0;           // t0 and associated error
    XYZVec          _pos;          // position of the time cluster   
    uint16_t	    _nsh;	  // number of straw hits
    art::Ptr<CaloCluster>      _caloCluster;  // associated calorimeter cluster: can be null
  };
   typedef std::vector<mu2e::TimeCluster> TimeClusterCollection;

} // namespace mu2e

#endif /* RecoDataProducts_TimeCluster_hh */
