//
// Base class for algorithm Objects used to cluster straw hits for background removal
//
//  David Brown (LBNL) 2013
//
#ifndef BkgClusterer_HH
#define BkgClusterer_HH
// data
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
//
namespace mu2e 
{
  class BkgClusterer
  {
  public:
    virtual ~BkgClusterer() {}
// initialize
    virtual void init() = 0;
// main function: given the straw hits and associated data, find the clusters.
    virtual void findClusters(BkgClusterCollection& clusters,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol) const = 0;
  };
}
#endif
