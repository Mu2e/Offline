#ifndef RecoDataProducts_BkgCluster_hh
#define RecoDataProducts_BkgCluster_hh
//
// Class to describe a cluster of straw hits in (transverse) space and time
// caused by low-energy electron (Compton, Delta-ray) background
// 
// Original author David Brown (LBNL), with contributions from David Ding 
//
// Mu2e includes
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/BkgClusterFlag.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include <vector>

namespace mu2e 
{
  // hit reference for us in the cluster.  The principle payload is the 'distance', which measures
  // how close this hit comes to the cluster in multiple dimensions and whose specific 
  // interpretation is algorithm dependent.
  
  struct BkgClusterHit 
  {
      BkgClusterHit() : _dist(0.0), _index(0) {}
      BkgClusterHit(uint16_t index, StrawHitFlag const& flag) : _dist(10000.0), _index(index), _flag(flag) {
        _flag.merge(StrawHitFlag::active); } // initial hits are active
      BkgClusterHit(float dist, uint16_t index, StrawHitFlag const& flag) : _dist(dist), _index(index), _flag(flag) {
        _flag.merge(StrawHitFlag::active); }

      float distance() const           { return _dist; }
      void  distance(float val)        {_dist=val; }
      uint16_t const& index() const    { return _index; }
      StrawHitFlag const& flag() const { return _flag; }

      float        _dist; // generalized distance from the cluster center
      uint16_t     _index; // reference to the straw hit
      StrawHitFlag _flag; // flag of this specific hit: only some bits are relevant
  };


  struct BkgCluster 
  {
      //Default hit count chosen for compuational efficiency
      BkgCluster() : _time(0.0)  {_hits.reserve(16);}   
      BkgCluster(XYZVec const& pos, float time) : _pos(pos), _time(time) {_hits.reserve(16);}

      XYZVec const&                     pos()  const { return _pos; }
      float                             time() const { return _time; }
      std::vector<BkgClusterHit> const& hits() const { return _hits; }
      BkgClusterFlag const&             flag() const { return _flag; }

      XYZVec                     _pos;  // ideally should be a 2d vec - FIXME
      float                      _time; 
      std::vector<BkgClusterHit> _hits; 
      BkgClusterFlag             _flag; 
  };

  typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;
} 
#endif

