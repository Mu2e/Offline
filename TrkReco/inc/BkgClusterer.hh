//
// Base class for algorithm Objects used to cluster straw hits for background removal
//
//  David Brown (LBNL) 2013
//
#ifndef BkgClusterer_HH
#define BkgClusterer_HH

#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

namespace mu2e 
{
    class BkgClusterer
    {
        public:
            virtual ~BkgClusterer() {}
            virtual void  init() = 0;
            virtual void  findClusters(BkgClusterCollection& preFilterClusters, BkgClusterCollection& postFilterClusters, 
                                       const ComboHitCollection& shcol, float mbtime, int iev) = 0;
            virtual float distance(const BkgCluster& cluster, const ComboHit& hit) const = 0; 
    };
}

#endif
