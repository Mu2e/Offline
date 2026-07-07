#ifndef BkgClusterer_HH
#define BkgClusterer_HH

#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e
{
  class BkgClusterer
  {
    public:
      virtual ~BkgClusterer() {};
      virtual void  init        () = 0;
      virtual void  findClusters(BkgClusterCollection& clusters, const ComboHitCollection& shcol) = 0;
      virtual void  classifyCluster(BkgCluster& clusters, const ComboHitCollection& shcol) = 0;

      virtual float distance    (const BkgCluster& cluster, const ComboHit& hit) const = 0;


  };
}

#endif
