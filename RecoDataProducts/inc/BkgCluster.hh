//
// Class to describe a cluster of straw hits in (transverse) space and time
// caused by low-energy electron (Compton, Delta-ray) background - D. Brown (LBL)
//
#ifndef RecoDataProducts_BkgCluster_hh
#define RecoDataProducts_BkgCluster_hh

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterFlag.hh"
#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include <vector>

namespace mu2e
{
   struct BkgCluster
   {
       //Default hit count chosen for compuational efficiency
       BkgCluster() {_hits.reserve(16);}
       BkgCluster(TwoDPoint point, float time) : _point(point), _time(time) {_hits.reserve(16);_cpoints.addPoint(_point,0);}

       float                        getKerasQ() const {return _kerasQ; }
       auto const&                  flag() const {return _flag; }
       auto const                   pos()  const {return _cpoints.point().pos3();  }
       auto &                       points()  {return _cpoints;  }
       auto const&                  points() const  {return _cpoints;  }
       auto const&                  time() const {return _time; }
       auto const&                  hits() const {return _hits; }
       auto &                       hits()       {return _hits; }

       void pos(XYZVectorF const& pos)               {_pos = pos;}
       void time(float time)                         {_time = time;}
       void addHit(unsigned val)                     {_hits.emplace_back(val);}
       void clearHits()                              {_hits.clear();}
       void setKerasQ(float kerasQ)                  {_kerasQ = kerasQ;}

       XYZVectorF               _pos;// ideally should be a 2d vec - FIXME
       TwoDPoint                _point;//initial point
       CombineTwoDPoints        _cpoints;//combined points
       float                    _time = 0.0;
       std::vector<unsigned>    _hits;
       BkgClusterFlag           _flag = BkgClusterFlag(BkgClusterFlag::update);
       float                    _kerasQ = -0.5; //result of keras result for the cluster
   };

   typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;

}
#endif

