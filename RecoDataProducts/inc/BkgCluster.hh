//
// Class to describe a cluster of straw hits in (transverse) space and time
// caused by low-energy electron (Compton, Delta-ray) background - D. Brown (LBL)
// 
#ifndef RecoDataProducts_BkgCluster_hh
#define RecoDataProducts_BkgCluster_hh

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/BkgClusterFlag.hh"
#include <vector>

namespace mu2e 
{
   struct BkgCluster 
   {
       //Default hit count chosen for compuational efficiency
       BkgCluster()                              : _pos(),    _time(0.0),  _hits(), _flag(BkgClusterFlag::update)  {_hits.reserve(16);}   
       BkgCluster(XYZVec const& pos, float time) : _pos(pos), _time(time), _hits(), _flag(BkgClusterFlag::update)  {_hits.reserve(16);}

       BkgClusterFlag const&        flag() const {return _flag; }
       XYZVec const&                pos()  const {return _pos;  }
       float                        time() const {return _time; }
       std::vector<unsigned> const& hits() const {return _hits; }
       std::vector<unsigned>&       hits()       {return _hits; }
       
       void pos(XYZVec const& pos)               {_pos = pos;}
       void time(float time)                     {_time = time;}
       void addHit(unsigned val)                 {_hits.emplace_back(val);}
       void clearHits()                          {_hits.clear();}

       XYZVec                _pos;  // ideally should be a 2d vec - FIXME
       float                 _time; 
       std::vector<unsigned> _hits; 
       BkgClusterFlag        _flag; 
   };

   typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;

} 
#endif

