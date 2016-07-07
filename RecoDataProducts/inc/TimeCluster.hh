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
#include <ostream>

// Mu2e includes
#include "TrkReco/inc/TrkDef.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  struct hitIndex {
    size_t _index; // index into the straw hit container
    int _ambig; // hit ambiguity.  0 means no ambiguity
    hitIndex() : _index(0),_ambig(0) {}
    hitIndex(size_t index,int ambig=0) : _index(index),_ambig(ambig) {}
    hitIndex& operator = (size_t index) { _index = index; return *this; }
    bool operator == (hitIndex const& other) const { return _index == other._index; }
  };
 
  struct TimeCluster{

    std::vector<hitIndex>    _strawHitIdxs;
    double                   _t0;   //time at z=z0
    double                   _errt0; //error asssociated to t0 (in case of CalPatRec is by default set to 1 ns)
    CLHEP::Hep3Vector        _pos; // position of the time cluster   
    art::Ptr<CaloCluster>    _caloCluster;

  public:

    TimeCluster():
     _t0   (0),
     _errt0(0)
     { }

  };

} // namespace mu2e

#endif /* RecoDataProducts_TimeCluster_hh */
