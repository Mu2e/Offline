//
// Object to perform helix fit to straw hits
//
// $Id: ClusterStrawHits.cc,v 1.1 2013/03/08 04:33:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2013/03/08 04:33:26 $
//
//
#include "TrkPatRec/inc/ClusterStrawHits.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "art/Framework/Services/Optional/TFileService.h"
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
// Root
#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <algorithm>
using namespace boost::accumulators;
namespace mu2e 
{
  unsigned StrawHitCluster::_currentid(0);

  struct ncomp : public std::binary_function<StrawHitCluster, StrawHitCluster, bool> {
    bool operator()(StrawHitCluster const& c1, StrawHitCluster const& c2) { return c2.hits().size() < c1.hits().size(); }
  };

  StrawHitCluster::StrawHitCluster(ClusterHit const& hit) : _pos(hit._pos),_time(hit._time), _id(++_currentid) {
    addHit(hit);
  }
// construct from the collection and a list of hits
  StrawHitCluster::StrawHitCluster(std::vector<ClusterHit> const& hits) : _id(++_currentid),_hits(hits) {
    for(std::vector<ClusterHit>::iterator ih=_hits.begin();ih!= _hits.end();++ih){
      accumulate(*ih);
    }
    updateCache();
  }

  // copy, etc
  StrawHitCluster::StrawHitCluster(StrawHitCluster const& other) : _pos(other._pos),_time(other._time),_id(other._id),
  _hits(other._hits){
    for(std::vector<ClusterHit>::iterator ih=_hits.begin();ih!= _hits.end();++ih){
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
      for(std::vector<ClusterHit>::iterator ih=_hits.begin();ih!= _hits.end();++ih){
	accumulate(*ih);
      }
    }
    return *this;
  }
  StrawHitCluster::~StrawHitCluster() {}
  
  void StrawHitCluster::addHit(ClusterHit const& hit,bool update) {
    _hits.push_back(hit);
    accumulate(hit);
    if(update)updateCache();
  }

  void StrawHitCluster::accumulate(ClusterHit const& hit) {
    _xacc(hit._pos.x());
    _yacc(hit._pos.y());
    _tacc(hit._time);
  }

  void StrawHitCluster::merge(StrawHitCluster const& other){
    std::vector<ClusterHit> const& ohits = other.hits();
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
    if(count(_xacc) > 2){
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
  std::vector<std::string> ClusterStrawHits::_nullsvec;

  ClusterStrawHits::ClusterStrawHits(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _bkgmask(pset.get<std::vector<std::string> >("BackgroundMask",_nullsvec)),
    _sigmask(pset.get<std::vector<std::string> >("SignalMask",_nullsvec)),
    _dseed(pset.get<double>("SeedDistance",10.0)), // # of sigma to define a new cluster
    _dhit(pset.get<unsigned>("HitDistance",6.0)), // # of sigma to add hits
    _dmerge(pset.get<double>("MergeDistance",6.0)), // # of sigma to merge clusters
    _dlarge(pset.get<double>("LargestDistance",10000.0)),
    _dd(pset.get<double>("ClusterDiameter",10.0)), // mm: the natural cluster size
    _dt(pset.get<double>("TimeDiference",60.0)), // nsec: the natural time spread
    _trms(pset.get<double>("TimeRMS",5.0)), //nsec: individual hit time resolution
    _srms(pset.get<double>("StereoRMS",5.0)), //mm: individual hit resolution
    _nsrms(pset.get<double>("NonStereoRMS",25.)), // mmm: twice the individual hit resolution
    _maxniter(pset.get<unsigned>("MaxNIterations",10)),
    _mode(static_cast<cmode>(pset.get<int>("ClusterMode",hitcluster)))
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
    if( dt < 2*_dt){
      bool stereo = hit._flag.hasAllProperties(stflag);
      double psig2 = stereo ? _srms2 : _nsrms2;
      double dperp = (cluster.pos()-hit._pos).perp();
      retval = sqrt(pow(std::max(0.0,dperp-_dd),2)/psig2 + pow(std::max(0.0,dt-_dt),2)/_trms2);
    }
    return retval;
  } 

  double ClusterStrawHits::distance(StrawHitCluster const& c1, StrawHitCluster const& c2) const {
    double retval = _dlarge;
    double dt = fabs(c1.time()-c2.time());
    if( dt < 2*_dt){
      double dperp = (c1.pos()-c2.pos()).perp();
      retval = sqrt(pow(std::max(0.0,dperp-_dd),2)/_srms2 + pow(std::max(0.0,dt-_dt),2)/_trms2);
    }
    return retval;
  } 

  void ClusterStrawHits::findClusters(StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol,
      std::list<StrawHitCluster>& clusters,
      std::vector<int>& clusterids) const {
// require consistency
    if(shcol.size() != shpcol.size() || shcol.size() != shfcol.size()){
      std::ostringstream os;
      os <<  " ClusterStrawHits: inconsistent collection lengths ";
      throw std::out_of_range( os.str() );
    }
// reset
    clusters.clear();
    StrawHitCluster::resetId();
    clusterids = std::vector<int>(shcol.size(),-2);
    // loop over the straw hits and create ClusterHits
    std::vector<ClusterHit> chits;
    chits.reserve(shcol.size());
    for(size_t ish=0;ish<shcol.size();++ish){
      if(shfcol[ish].hasAllProperties(_sigmask) && !shfcol[ish].hasAnyProperty(_bkgmask)){
	chits.push_back(ClusterHit(shpcol[ish],shcol[ish],shfcol[ish],ish));
      }
    }
// The 1st hit becomes the first cluster
    clusters.push_back(StrawHitCluster(chits[0]));
// now iterate
    bool changed(true);
    unsigned niter(0);
    while (changed && niter < _maxniter){
// update the positions after each added hit in the first iteration (only)
//      bool update = niter==0;
      bool update(false);
// update the cluster positions, and record any changes in hit assignments to clusters
      changed = formClusters(chits,clusters,clusterids,update);
      ++niter;
    }
  }

  bool ClusterStrawHits::formClusters(std::vector<ClusterHit> const& chits, std::list<StrawHitCluster>& clusters,
      std::vector<int>& clusterids, bool update) const {
    bool changed(false);
// clear out the hits
    for(std::list<StrawHitCluster>::iterator ic=clusters.begin();ic!=clusters.end();++ic){
      ic->clearHits();
    }
 //loop over clusterhits, seeding clusters or appending hits to them, as appropriate
    for(size_t ich=0;ich<chits.size();++ich){
      ClusterHit const& chit = chits[ich];
      double mindist(FLT_MAX);
      std::list<StrawHitCluster>::iterator minc = clusters.end();
      for(std::list<StrawHitCluster>::iterator ic=clusters.begin();ic!=clusters.end();++ic){
	double dist = distance(*ic,chit);
	if(dist < mindist){
	  mindist = dist;
	  minc = ic;
	}
      }
      if(mindist < _dhit){
// append this hit to the cluster, and update the cache
	minc->addHit(chit,update);
      } else if(mindist > _dseed) {
// seed a new cluster
	clusters.push_back(StrawHitCluster(chit));
	minc = --clusters.end();
      }
// see if the cluster assignment of this hit has changed
      if(minc != clusters.end()){
	changed |= clusterids[chit._index] != minc->id();
	clusterids[chit._index] = minc->id();
      } else {
	changed |= clusterids[chit._index] >=0;
	clusterids[chit._index] = -1;
      }
    }
// final cluster update
    for(std::list<StrawHitCluster>::iterator ic=clusters.begin();ic!=clusters.end();++ic){
      ic->updateCache();
    }
// sort the list
//    std::sort(clusters.begin(),clusters.end(),ncomp());
// merge the clusters
//    unsigned nbefore = clusters.size();
    changed |= mergeClusters(clusters);
//    unsigned nafter = clusters.size();
    return changed;
  }

  bool ClusterStrawHits::mergeClusters(std::list<StrawHitCluster>& clusters) const {
    bool merged(false);
    for(std::list<StrawHitCluster>::iterator ic=clusters.begin();ic!=clusters.end();++ic){
      std::list<StrawHitCluster>::iterator jc=ic;jc++;
      double mindist(FLT_MAX);
      std::list<StrawHitCluster>::iterator imin=clusters.end();
      for(;jc!=clusters.end();++jc){
	double dist = distance(*ic,*jc);
	if(dist < mindist){
	  mindist = dist;
	  imin = jc;
	}
      }
      if(mindist < _dmerge){
      // merge the smaller cluster into the bigger
	ic->merge(*imin);
	clusters.erase(imin);
	merged = true;
      }
    }
    return merged;
  }

}

