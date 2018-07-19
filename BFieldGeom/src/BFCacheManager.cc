// Andrei Gaponenko, 2012

#include "BFieldGeom/inc/BFCacheManager.hh"

namespace mu2e {

    BFCacheManager::BFCacheManager():
    innerForLastInner(0),
    innerForLastOuter(0),
    counter(0)
    {
        outerCache.insert( std::make_pair<std::shared_ptr<const BFMap>>(0, CacheElement(0, MapList())) );
        CacheType::const_iterator p = outerCache.find(0);
        assert(p != outerCache.end());
        innerForLastOuter = &p->second;
    }

    void BFCacheManager::setMaps(const MapContainerType& innerMaps,
                                 const MapContainerType& outerMaps) {
        typedef MapContainerType::const_iterator Iter;

        // All inner maps in the input order
        MapList defaultInnerList;
        for (Iter i = innerMaps.begin(); i != innerMaps.end(); ++i) {
            defaultInnerList.push_back(*i);
        }

        // The fixed-order outer map list
        for (Iter i = outerMaps.begin(); i != outerMaps.end(); ++i) {
            this->outer.push_back(*i);
        }

        // Now populate the cache lookup structures

        for (Iter i = innerMaps.begin(); i != innerMaps.end(); ++i) {
            // or can assign a dedicated innerList for this map, e.g. using hints from FHICL
            innerCache.insert(std::make_pair(*i, CacheElement(*i, defaultInnerList)));
        }

        for (Iter i = outerMaps.begin(); i != outerMaps.end(); ++i) {
            // or can assign a dedicated innerList for this map, e.g. using hints from FHICL
            outerCache.insert(std::make_pair(*i, CacheElement(*i, defaultInnerList)));
        }

        CacheType::iterator p = outerCache.find(0);
        assert(p != outerCache.end());
        p->second = CacheElement(0, defaultInnerList);
        innerForLastOuter = &p->second;
    }
}  // namespace mu2e
