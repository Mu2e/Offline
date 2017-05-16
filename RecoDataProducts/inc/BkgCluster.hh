#ifndef RecoDataProducts_BkgCluster_hh
#define RecoDataProducts_BkgCluster_hh
//
// Class to describe a cluster of straw hits in (transverse) space and time
// caused by low-energy electron (Compton, Delta-ray) background
// 
// Original author David Brown (LBNL), with contributions from David Ding 
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitIndex.hh"
namespace mu2e 
{
// hit reference for us in the cluster.  The only added payload is the 'distance', which measures
// how close this hit comes to the cluster in multiple dimensions and whose specific 
// interpretation is algorithm dependent.
  struct BkgClusterHit {
    float distance() const { return _dist; }
    StrawHitIndex const& index() const { return _index; }
    float _dist; // generalized distance from the cluster center
    StrawHitIndex _index; // reference to the straw hit, 
  };
// the cluster itself
  struct BkgCluster {
// construct from a position and time
    BkgCluster(CLHEP::Hep3Vector const& pos, double time);
// clear the hits and position information
    void clearHits() { _hits.clear(); _pos = CLHEP::Hep3Vector();}
// accessors
    CLHEP::Hep3Vector const& pos() const { return _pos; }
    float time() const { return _time; }
    std::vector<BkgClusterHit> const& hits() const { return _hits; }
// data members
    CLHEP::Hep3Vector _pos; // nominal position.  Only the xy values have meaning, ideally this would be a dedicated 2-d vector object FIXME!
    float _time;  //nominal time
    std::vector<BkgClusterHit> _hits; // vector of associated hits
  };
  // define collection type
  typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;
} 
#endif

