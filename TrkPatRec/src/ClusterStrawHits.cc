//
// Object to perform helix fit to straw hits
//
#include "TrkPatRec/inc/ClusterStrawHits.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost

#include <vector>
#include <algorithm>
using namespace boost::accumulators;
using namespace std;
namespace mu2e
{
  unsigned StrawHitCluster::_currentid(0);

  struct ncomp : public binary_function<StrawHitCluster, StrawHitCluster, bool> {
    bool operator()(StrawHitCluster const& c1, StrawHitCluster const& c2) { return c2.hits().size() < c1.hits().size(); }
  };

  StrawHitCluster::StrawHitCluster(ClusterHit const& hit) : _pos(hit._pos),_time(hit._time), _id(++_currentid) {
    addHit(hit);
  }
// construct from the collection and a list of hits
  StrawHitCluster::StrawHitCluster(vector<ClusterHit> const& hits) : _id(++_currentid),_hits(hits) {
    for(vector<ClusterHit>::iterator ih=_hits.begin();ih!= _hits.end();++ih){
      accumulate(*ih);
    }
    updateCache();
  }

  // copy, etc
  StrawHitCluster::StrawHitCluster(StrawHitCluster const& other) : _pos(other._pos),_time(other._time),_id(other._id),
  _hits(other._hits){
    for(vector<ClusterHit>::iterator ih=_hits.begin();ih!= _hits.end();++ih){
      accumulate(*ih);
    }
  }


  StrawHitCluster& StrawHitCluster::operator =(StrawHitCluster const& other){
    if(&other != this){
      _pos = other._pos;
      _time = other._time;
      _id = other._id;
      _hits.clear();
      _hits = other._hits;
      for(vector<ClusterHit>::iterator ih=_hits.begin();ih!= _hits.end();++ih){
        accumulate(*ih);
      }
    }
    return *this;
  }
  StrawHitCluster::~StrawHitCluster() {}

  void StrawHitCluster::addHit(ClusterHit const& hit,double dist,bool update) {
    _hits.push_back(hit);
    _hits.back()._dist = dist;
    accumulate(hit);
    if(update)updateCache();
  }

  void StrawHitCluster::accumulate(ClusterHit const& hit) {
    _xacc(hit._pos.x());
    _yacc(hit._pos.y());
    _tacc(hit._time);
  }

  void StrawHitCluster::merge(StrawHitCluster const& other){
    vector<ClusterHit> const& ohits = other.hits();
    for(size_t oh=0;oh<ohits.size();++oh){
      addHit(ohits[oh]);
    }
    updateCache();
  }

  void StrawHitCluster::updateCache() {
//    _pos.x() = mean(_xacc);
 //   _pos.y() = mean(_yacc);
//    _time  = mean(_tacc);
// there is a bug in boost, the median is not calculated correctly if there are only a few entries
    if( boost::accumulators::extract::count(_xacc) > 2){
      double x = median(_xacc);
      double y = median(_yacc);
      _pos = CLHEP::Hep3Vector(x,y,0.0);
      _time  = median(_tacc);
    } else if(_hits.size()>0){
      _pos = _hits[0]._pos;
      _time = _hits[0]._time;
    }
  }

  void StrawHitCluster::clearHits() {
    _hits.clear();
// reset the accumulators by brute force
    _xacc = accumulator_set<double, stats<tag::median(with_p_square_quantile) > >();
    _yacc = accumulator_set<double, stats<tag::median(with_p_square_quantile) > >();
    _tacc = accumulator_set<double, stats<tag::median(with_p_square_quantile) > >();
//    _xacc = accumulator_set<double, stats<tag::mean > >();
//    _yacc = accumulator_set<double, stats<tag::mean > >();
//    _tacc = accumulator_set<double, stats<tag::mean > >();
  }
  ClusterStrawHits::ClusterStrawHits(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _bkgmask(pset.get<vector<string> >("BackgroundMask",vector<string>())),
    _sigmask(pset.get<vector<string> >("SignalMask",vector<string>())),
    _dseed(pset.get<double>("SeedDistance",10.0)), // # of sigma to define a new cluster
    _dhit(pset.get<unsigned>("HitDistance",5.0)), // # of sigma to add hits
    _dmerge(pset.get<double>("MergeDistance",3.0)), // # of sigma to merge clusters
    _dlarge(pset.get<double>("LargestDistance",10000.0)),
    _dd(pset.get<double>("ClusterDiameter",10.0)), // mm: the natural cluster size
    _dt(pset.get<double>("TimeDifference",10.0)), // nsec: the natural time spread
    _maxdt(pset.get<double>("MaxTimeDifference",30.0)), // Maximum time difference
    _trms(pset.get<double>("TimeRMS",10.0)), //nsec: individual hit time resolution
    _srms(pset.get<double>("StereoRMS",6.0)), //mm: individual stereo hit resolution
    _nsrms(pset.get<double>("NonStereoRMS",20.)), // mmm: individual non-stereo hit resolution
    _maxniter(pset.get<unsigned>("MaxNIterations",50)),
    _maxnchanged(pset.get<unsigned>("MaxNChanged",10)),
    _mode(static_cast<cmode>(pset.get<int>("ClusterMode",hitcluster))),
    _updatefirst(pset.get<bool>("UpdateFirstIteration",false))
  {
  // compute the squares
    _srms2 = _srms*_srms;
    _nsrms2 = _nsrms*_nsrms;
    _trms2 = _trms*_trms;
  }

  ClusterStrawHits::~ClusterStrawHits()
  {}

  double ClusterStrawHits::distance(StrawHitCluster const& cluster, ClusterHit const& hit) const {
    static StrawHitFlag stflag(StrawHitFlag::stereo);
    double retval = _dlarge;
    double dt = fabs(hit._time-cluster.time());
    if( dt < _maxdt){
      bool stereo = hit._flag.hasAllProperties(stflag);
      double psig2 = stereo ? _srms2 : _nsrms2;
      if(_mode == hitcluster){
        double dperp = (cluster.pos()-hit._pos).perp();
        retval = sqrt(pow(fmax(0.0,dperp-_dd),2)/psig2 + pow(fmax(0.0,dt-_dt),2)/_trms2);
      } else  {
// loop over all hits in the cluster, and determine the distance to the closests one
        double mindp2(FLT_MAX);
        for(size_t ih=0;ih<cluster.hits().size();++ih){
          double dp2 = (cluster.hits()[ih]._pos-hit._pos).perp();
          if(dp2 < mindp2)mindp2 = dp2;
        }
        retval = sqrt(pow(fmax(0.0,sqrt(mindp2)-_dd),2)/psig2 + pow(fmax(0.0,dt-_dt),2)/_trms2);
      }
    }
    return retval;
  }

  double ClusterStrawHits::distance(StrawHitCluster const& c1, StrawHitCluster const& c2) const {
    double retval = _dlarge;
    double dt = fabs(c1.time()-c2.time());
    if( dt < _maxdt){
      if(_mode == hitcluster) {
        double dperp = (c1.pos()-c2.pos()).perp();
        retval = sqrt(pow(fmax(0.0,dperp-_dd),2)/_srms2 + pow(fmax(0.0,dt-_dt),2)/_trms2);
      } else {
        double mindp2(FLT_MAX);
        for(size_t ih=0;ih<c1.hits().size();++ih){
          for(size_t jh=0;jh<c2.hits().size();++jh){
            double dp2 = (c1.hits()[ih]._pos-c2.hits()[jh]._pos).perp();
            if(dp2 < mindp2)mindp2 = dp2;
          }
        }
        retval = sqrt(pow(fmax(0.0,sqrt(mindp2)-_dd),2)/_srms2 + pow(fmax(0.0,dt-_dt),2)/_trms2);
      }
    }
    return retval;
  }

  void ClusterStrawHits::findClusters(StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol,
      StrawHitClusterList& clusters) const {
// require consistency
    if(shcol.size() != shpcol.size() || shcol.size() != shfcol.size()){
      ostringstream os;
      os <<  " ClusterStrawHits: inconsistent collection lengths ";
      throw out_of_range( os.str() );
    }
// reset
    clusters._clist.clear();
    StrawHitCluster::resetId();
    clusters._cids = vector<int>(shcol.size(),-2);
    // loop over the straw hits and create ClusterHits
    vector<ClusterHit> chits;
    chits.reserve(shcol.size());
    for(size_t ish=0;ish<shcol.size();++ish){
      if(shfcol[ish].hasAllProperties(_sigmask) && !shfcol[ish].hasAnyProperty(_bkgmask)){
        chits.push_back(ClusterHit(shpcol[ish],shcol[ish],shfcol[ish],ish));
      }
    }
// The 1st hit becomes the first cluster
    clusters._niter = 0;
    clusters._nchanged = shcol.size();
    clusters._nhits = chits.size();
    if(chits.size() > 0){
      clusters._clist.push_back(StrawHitCluster(chits[0]));
      clusters._cids[chits[0]._index] = clusters._clist.begin()->id();
      // now iterate
      while (clusters._nchanged > _maxnchanged && clusters._niter < _maxniter){
        // update the cluster positions, and record any changes in hit assignments to clusters
        clusters._nchanged = formClusters(chits,clusters);
        ++clusters._niter;
      }
    }
  }

  unsigned ClusterStrawHits::formClusters(vector<ClusterHit> const& chits, StrawHitClusterList& clusters) const {
    unsigned  nchanged(0);
    bool update=_updatefirst && clusters._niter==0;
// clear out the hits
    for(list<StrawHitCluster>::iterator ic=clusters._clist.begin();ic!=clusters._clist.end();++ic){
      ic->clearHits();
    }
 //loop over clusterhits, seeding clusters or appending hits to them, as appropriate
    for(size_t ich=0;ich<chits.size();++ich){
      ClusterHit const& chit = chits[ich];
      double mindist(FLT_MAX);
      list<StrawHitCluster>::iterator minc = clusters._clist.end();
      for(list<StrawHitCluster>::iterator ic=clusters._clist.begin();ic!=clusters._clist.end();++ic){
        double dist = distance(*ic,chit);
        if(dist < mindist){
          mindist = dist;
          minc = ic;
        }
      }
      if(mindist < _dhit){
// append this hit to the cluster, and update the cache
        minc->addHit(chit,mindist,update);
      } else if(mindist > _dseed) {
// seed a new cluster
        clusters._clist.push_back(StrawHitCluster(chit));
        minc = --clusters._clist.end();
      }
// see if the cluster assignment of this hit has changed
      if(minc != clusters._clist.end()){
        if(clusters._cids[chit._index] != (int)minc->id())++nchanged;
        clusters._cids[chit._index] = (int)minc->id();
      } else {
// unassigned hit
        if(clusters._cids[chit._index] >=0)++nchanged;
        clusters._cids[chit._index] = -1;
      }
    }
// final cluster update
    for(list<StrawHitCluster>::iterator ic=clusters._clist.begin();ic!=clusters._clist.end();++ic){
      ic->updateCache();
    }
// sort the list
//    sort(clusters.begin(),clusters.end(),ncomp());
// merge the clusters
//    unsigned nbefore = clusters.size();
    nchanged += mergeClusters(clusters);
//    unsigned nafter = clusters.size();
    return nchanged;
  }

  unsigned ClusterStrawHits::mergeClusters(StrawHitClusterList& clusters) const {
    unsigned nmerged(0);
    for(list<StrawHitCluster>::iterator ic=clusters._clist.begin();ic!=clusters._clist.end();++ic){
      list<StrawHitCluster>::iterator jc=ic;jc++;
      double mindist(FLT_MAX);
      list<StrawHitCluster>::iterator imin=clusters._clist.end();
      for(;jc!=clusters._clist.end();++jc){
        double dist = distance(*ic,*jc);
        if(dist < mindist){
          mindist = dist;
          imin = jc;
        }
      }
      if(mindist < _dmerge){
      // merge the smaller cluster into the bigger
        ic->merge(*imin);
        // reassign the index
        for(size_t ih=0;ih<imin->hits().size();++ih){
          clusters._cids[imin->hits()[ih]._index] = ic->id();
        }
        // erase the 2nd cluster
        nmerged += imin->hits().size();
        clusters._clist.erase(imin);
      }
    }
    return nmerged;
  }

}

