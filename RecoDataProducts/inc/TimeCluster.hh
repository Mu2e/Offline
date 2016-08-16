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
#include "art/Persistency/Common/Ptr.h"
// Mu2e includes
//#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CLHEP/Vector/ThreeVector.h"
// BTrk includes
#include "BTrk/TrkBase/TrkT0.hh"

namespace mu2e {
  class CaloCluster;
  struct hitIndex {
    size_t _index; // index into the straw hit container
    int _ambig; // hit ambiguity.  0 means no ambiguity
    hitIndex() : _index(0),_ambig(0) {}
    hitIndex(size_t index,int ambig=0) : _index(index),_ambig(ambig) {}
    hitIndex& operator = (size_t index) { _index = index; return *this; }
    bool operator == (hitIndex const& other) const { return _index == other._index; }
  };
 
  struct TimeCluster{

    std::vector<hitIndex>    _strawHitIdxs; // associated straw hits: can be empty
    TrkT0		     _t0; // t0 and associated error
    CLHEP::Hep3Vector        _pos; // position of the time cluster   
    art::Ptr<CaloCluster>    _caloCluster; // associated calorimeter cluster: can be null

  };

} // namespace mu2e

#endif /* RecoDataProducts_TimeCluster_hh */
