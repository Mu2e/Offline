//
// Class to describe a cluster of straw hits in (transverse) space and time
// caused by low-energy electron (Compton, Delta-ray) background - D. Brown (LBL)
//
#ifndef RecoDataProducts_BkgCluster_hh
#define RecoDataProducts_BkgCluster_hh

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterFlag.hh"
#include <vector>

namespace mu2e
{
   struct BkgCluster
   {
       //Default hit count chosen for compuational efficiency
       BkgCluster() {_hits.reserve(16);}
       BkgCluster(XYZVectorF const& pos, float time) : _pos(pos), _time(time) {_hits.reserve(16);}


       float                        getKerasQ() const {return _kerasQ; }
       auto const&        flag() const {return _flag; }
       auto const&                pos()  const {return _pos;  }
       auto const& time() const {return _time; }
       auto const& hits() const {return _hits; }
       auto &       hits()       {return _hits; }

       void pos(XYZVectorF const& pos)               {_pos = pos;}
       void time(float time)                     {_time = time;}
       void addHit(unsigned val)                 {_hits.emplace_back(val);}
       void clearHits()                          {_hits.clear();}
       void setKerasQ(float kerasQ)                  {_kerasQ = kerasQ;}

       XYZVectorF            _pos;  // ideally should be a 2d vec - FIXME
       float                 _time = 0.0;
       std::vector<unsigned> _hits;
       BkgClusterFlag        _flag = BkgClusterFlag(BkgClusterFlag::update);
       float                 _kerasQ = -0.5; //result of keras result for the cluster
   };

   typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;

}
#endif

