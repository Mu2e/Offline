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
// c++ includes
#include <vector>
namespace mu2e 
{
// hit reference for us in the cluster.  The principle payload is the 'distance', which measures
// how close this hit comes to the cluster in multiple dimensions and whose specific 
// interpretation is algorithm dependent.
  struct BkgClusterHit {
    BkgClusterHit() : _dist(0.0), _index(0) {}
    BkgClusterHit(float dist, uint16_t index, StrawHitFlag const& flag) : _dist(dist), _index(index), _flag(flag) {
      _flag.merge(StrawHitFlag::active); }// initial hits are active
    float distance() const { return _dist; }
    uint16_t const& index() const { return _index; }
    StrawHitFlag const& flag() const { return _flag; }
    float _dist; // generalized distance from the cluster center
    uint16_t _index; // reference to the straw hit
    StrawHitFlag _flag; // flag of this specific hit: only some bits are relevant
  };
// the cluster itself
  struct BkgCluster {
  // default constructor.  Default hit count chosen for compuational efficiency
    BkgCluster() : _time(0.0)  { _hits.reserve(16); }
// construct from a position and time
    BkgCluster(XYZVec const& pos, float time) : _pos(pos), _time(time) { _hits.reserve(16); }
// accessors
    XYZVec const& pos() const { return _pos; }
    float time() const { return _time; }
    std::vector<BkgClusterHit> const& hits() const { return _hits; }
    BkgClusterFlag const& flag() const { return _flag; }
// data members
    XYZVec _pos; // nominal position.  Only the xy values have meaning, ideally this would be a dedicated 2-d vector object FIXME!
    float _time;  //nominal time
    std::vector<BkgClusterHit> _hits; // associated hits
    BkgClusterFlag _flag; // characterize this cluster
  };

  typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;
} 
#endif

