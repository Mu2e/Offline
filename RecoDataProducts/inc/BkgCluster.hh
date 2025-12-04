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
       enum  distMethod {spatial,chi2};

       //Default hit count chosen for compuational efficiency
       BkgCluster() {_hits.reserve(16);}
       BkgCluster(TwoDPoint point, distMethod method) : _point(point), _distMethod(method){_hits.reserve(16);_cpoints.addPoint(_point,0);_pos = point.pos3();}
       BkgCluster(XYZVectorF const& pos, float time, distMethod method) : _pos(pos), _time(time), _distMethod(method) {_hits.reserve(16);}

       float                        getKerasQ() const {return _kerasQ; }
       auto const&                  getDistMethod() const {return _distMethod; }
       auto const&                  flag() const {return _flag; }
       auto const                   pos()  const {return _pos;  }
       auto &                       points()  {return _cpoints;  }
       auto const&                  points() const  {return _cpoints;  }
       auto const&                  time() const {return _time; }
       auto const&                  hits() const {return _hits; }
       auto const&                  hitposition() const {return _hitpositions; }
       auto &                       hitposition() {return _hitpositions; }
       auto &                       hits()       {return _hits; }

       void pos(XYZVectorF const& pos)                                      {_pos = pos;}
       void time(float time)                                                {_time = time;}
       void clearHits()                                                     {_hits.clear();}
       void setKerasQ(float kerasQ)                                         {_kerasQ = kerasQ;}
       void addHit(unsigned val)                                            {_hits.emplace_back(val);}
       void addHit(unsigned val, TwoDPoint point) {
        _hits.push_back(val);
        _cpoints.addPoint(point,_cpoints.nPoints());
        _pos = _cpoints.point().pos3();
       }
       void addHitPosition(XYZVectorF const& pos) {
        _hitpositions.push_back(pos);
       }

       XYZVectorF               _pos;// ideally should be a 2d vec - FIXME
       TwoDPoint                _point;//initial point
       CombineTwoDPoints        _cpoints;//combined points
       float                    _time = 0.0;//cluster time
       std::vector<unsigned>    _hits;
       std::vector<XYZVectorF>    _hitpositions;
       BkgClusterFlag           _flag = BkgClusterFlag(BkgClusterFlag::update);
       float                    _kerasQ = -0.5; //result of keras result for the cluster
       distMethod               _distMethod;//which distMethod used to create the cluster
   };

   typedef std::vector<mu2e::BkgCluster> BkgClusterCollection;

}
#endif
