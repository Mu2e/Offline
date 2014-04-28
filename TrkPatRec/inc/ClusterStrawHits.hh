//
// Object to cluster straw hits, used in background removal and track fitting
//
// $Id: ClusterStrawHits.hh,v 1.6 2014/04/28 13:51:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/04/28 13:51:26 $
//
#ifndef ClusterStrawHits_HH
#define ClusterStrawHits_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/mean.hpp>

//root
class TH1F;
// C+

namespace mu2e 
{
// struct to define a hit used in a cluster
  struct ClusterHit {
    ClusterHit(CLHEP::Hep3Vector const& pos,double time,StrawHitFlag const& flag,size_t index) : _pos(pos),_time(time),
      _flag(flag),_index(index){}
    ClusterHit(StrawHitPosition const& shpos, StrawHit const& sh,StrawHitFlag const& flag,size_t index) : _pos(shpos.pos()),
    _time(sh.time()),_flag(flag),_index(index) {
      _flag.merge(shpos.flag());
    }
    // allow sorting
    bool operator < (ClusterHit const& other) const { return _index < other._index; }
    bool operator == (ClusterHit const& other) const { return _index == other._index; }
    CLHEP::Hep3Vector _pos;
    double _time;
    double _dist;
    StrawHitFlag _flag;
    size_t _index;
  };
// struct to definea cluster
  class StrawHitCluster {
    public:
// construct from a position
    StrawHitCluster(CLHEP::Hep3Vector const& pos, double time);
// construct from an initial seed hit
    StrawHitCluster(ClusterHit const& hit);
// construct from the collection and a list of hits
    StrawHitCluster(std::vector<ClusterHit> const& hits);
// copy, etc
    StrawHitCluster(StrawHitCluster const& other);
    StrawHitCluster& operator =(StrawHitCluster const& other);
    virtual ~StrawHitCluster();
// identity functions
    bool operator ==(StrawHitCluster const& other) const { return _id == other._id; }
    bool operator <(StrawHitCluster const& other) const { return _id < other._id; }
// append a hit, optionally updating the cache
    void addHit(ClusterHit const&, double dist=0.0,bool update=false);
// clear the hits; this leaves the position information
    void clearHits();
  // merge with another cluster
    void merge(StrawHitCluster const& other);
// accessors
    CLHEP::Hep3Vector const& pos() const { return _pos; }
    double time() const { return _time; }
    unsigned id() const { return _id; }
    std::vector<ClusterHit> const& hits() const { return _hits; }
    std::vector<ClusterHit>&  hits() { return _hits; }
// evaluate the accumulators to find the position and time.
    void updateCache();
// reset static for id assignment
    static void resetId() { _currentid=0; }
    private:
// cache values
    CLHEP::Hep3Vector _pos;
    double _time;
    unsigned _id;
    static unsigned _currentid;
// use a vector for storage.  Removal is not allowed
    std::vector<ClusterHit> _hits; 
    void accumulate(ClusterHit const& hit);
// accumulators
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::median(boost::accumulators::with_p_square_quantile) > > _xacc, _yacc, _tacc;
//    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::mean > > _xacc, _yacc, _tacc;
  };

  struct StrawHitClusterList {
    unsigned _nhits; // number of hits used as input (after filter)
    unsigned _niter; // number of clustering iterations
    unsigned _nchanged; // number of changed hits at last iteration
    std::list<StrawHitCluster> _clist; // list of clusters
    std::vector<int> _cids; // map of which hits go to which clusters (parallel with StrawHitCollection)
  };
 
  class ClusterStrawHits
  {
  public:
// define the mode: distance based on cluster-hit separation, or hit-hit separation
    enum cmode {hitcluster=0,hithit=1};
// parameter set should be passed in on construction
    explicit ClusterStrawHits(fhicl::ParameterSet const&);
    virtual ~ClusterStrawHits();
// main function: given the straw hits and associated data, find the clusters.
     void findClusters(StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol,
      StrawHitClusterList& clusters) const;
// functions to define the distance between a hit and cluster and 2 clusters; this is the variable against which
// the clustering thresholds are applied
    double distance(StrawHitCluster const&, ClusterHit const&) const;
    double distance(StrawHitCluster const&, StrawHitCluster const&) const;
  protected:
    unsigned mergeClusters(StrawHitClusterList& clusters) const;
    unsigned formClusters(std::vector<ClusterHit> const& chits,StrawHitClusterList& clusters) const;
// utlity functions
  private:
// configuration parameters
    int _diag,_debug;
    StrawHitFlag _bkgmask; // mask for background hits
    StrawHitFlag _sigmask; // mask for selecting signals
    double _dseed; // Minimum separation to seed a new cluster
    double _dhit; // Maximum separation to include a hit in a cluster
    double _dmerge; // distance to merge 2 clusters
    double _dlarge; // largest distance
    double _dd; // cluster diameter
    double _dt; // natural time spread
    double _maxdt; // maximum time difference
    double _trms; // time RMS
    double _srms,_nsrms; // Spatial RMS for stereo, non-stereo hits
    double _trms2, _srms2, _nsrms2; // squares of rms
    unsigned _maxniter; // maximum number of iterations
    unsigned _maxnchanged; // maximum # of changed hits to consider 'converged'
    cmode _mode; // clustering mode
    bool _updatefirst; // update positions continuously on the 1st iteration
    static std::vector<std::string> _nullsvec;
  };
}
#endif
