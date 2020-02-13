//
// Class to describe the content of a bkg cluster hit - B. Echenard
//

#ifndef RecoDataProducts_BkgClusterHit_hh
#define RecoDataProducts_BkgClusterHit_hh

#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include <vector>

namespace mu2e 
{
   struct BkgClusterHit 
   {
       BkgClusterHit() : _dist(0.0) {}

       BkgClusterHit(float dist, const StrawHitFlag& flag) : _dist(dist), _flag(flag) 
       {_flag.merge(StrawHitFlag::active);}

       float distance() const           { return _dist; }
       const StrawHitFlag& flag() const { return _flag; }

       float        _dist;  
       StrawHitFlag _flag; 
   };

   typedef std::vector<mu2e::BkgClusterHit> BkgClusterHitCollection;

} 
#endif

