//
// Object to perform helix fit to straw hits
//
#include "TrkReco/inc/TLTClusterer.hh"
// services
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
// utilities
#include "GeneralUtilities/inc/Angles.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_median.hpp>
// root, for diagnostics
#include "TTree.h"
#include "Math/VectorUtil.h"
// C++
#include <vector>
#include <set>
#include <algorithm>
using namespace boost::accumulators;
using namespace std;
using namespace ROOT::Math::VectorUtil;
namespace mu2e
{
 TLTClusterer::TLTClusterer(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _palg(static_cast<PosAlgorithm>(pset.get<int>("ClusterPositionAlgorithm",median))),
    _testflag(pset.get<bool>("TestFlag")),
    _bkgmask(pset.get<vector<string> >("BackgroundMask",vector<string>())),
    _sigmask(pset.get<vector<string> >("SignalMask",vector<string>())),
    _dseed(pset.get<float>("SeedDistance",100.0)), // minimum 'chisquared' to define a new cluster
    _dhit(pset.get<unsigned>("HitDistance",25.0)), // maximum 'chisquared' to add hit to a cluster
    _dd(pset.get<float>("ClusterDiameter",10.0)), // the natural cluster transverse size (mm)
    _dt(pset.get<float>("TimeDifference",30.0)), // the natural time spread (nsec)
    _maxdt(pset.get<float>("MaxTimeDifference",50.0)), // Maximum time difference (nsec)
    _maxdsum(pset.get<float>("MaxDistanceSum",100.0)), // iteration convergence
    _maxniter(pset.get<unsigned>("MaxNIterations",50)),
    _stereoinit(pset.get<bool>("StereoInit",false)), // initialize using stereo hits
    _minerr(pset.get<float>("MinHitError",5.0)), // corresponds to an error of the straw diameter
    _idiag(0)
  {
    // cache some values to avoid repeating FP operations on fixed numbers
    float trms(pset.get<float>("TimeRMS",2.0)); // effective individual hit time resolution (nsec)
    _trms2 = trms*trms;
    float maxdist(pset.get<float>("MaxDistance",50.0)); // Maximum transverse distance (mm)
    _md2 = maxdist*maxdist;
    _dd2 = _dd*_dd;
    _stereo = StrawHitFlag(StrawHitFlag::stereo);
    _stereo.merge(_sigmask);
    _maxwt = 1.0/_minerr;
 }

  TLTClusterer::~TLTClusterer()
  {}

  float TLTClusterer::distance(BkgCluster const& cluster, ComboHit const& ch) const {
    float retval = _dseed+1.0; // default return is above seed threshold
    // compute the simplest parts first to avoid expensive calculations on distant hits
    float dt = fabs(ch.time()-cluster.time());
    if( dt < _maxdt){
// compute spatial distance
      XYZVec psep = PerpVector(ch.pos()-cluster.pos(),Geom::ZDir());
      // work in squared magnitudes for efficiency
      float d2 = psep.mag2();
      if( d2 < _md2) {
	retval = 0.0; // if separation is inside natural size, count distance as 0
	// only count differences above the natural hit size (drift time, cluster size)
	if(dt > _dt){
	  float tdist = dt -_dt;
	  retval += tdist*tdist/_trms2;
	}
	if(d2 > _dd2){
	  // project separation along wire and transverse directions.  Only count distance beyond
	  // the natural cluster size
	  float dw = std::max(float(0.0),ch.wdir().Dot(psep)-_dd)/ch.posRes(ComboHit::wire);
	  XYZVec that(-ch.wdir().y(),ch.wdir().x(),0.0);
	  // minimum error is always larger than the transverse error
	  float dp = std::max(float(0.0),that.Dot(psep)-_dd)/_minerr;
	// add these contributions
	  retval += dw*dw + dp*dp;
	}
      }
    }
    return retval;
  }

  void TLTClusterer::init() {
   // setup diagnostics
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _idiag = tfs->make<TTree>("idiag","iteration diagnostics");
      _idiag->Branch("niter",&_niter,"niter/i");
      _idiag->Branch("nmerge",&_nmerge,"nmerge/i");
      _idiag->Branch("nclu",&_nclu,"nclu/i");
      _idiag->Branch("nchits",&_nchits,"nchits/i");
      _idiag->Branch("nhits",&_nhits,"nhits/i");
      _idiag->Branch("odist",&_odist,"odist/D");
      _idiag->Branch("tdist",&_tdist,"tdist/D");
    }
  }

  void TLTClusterer::findClusters(BkgClusterCollection& clusters,
      ComboHitCollection const& chcol) {
// initialize the clusters
    initClusters(clusters, chcol);
    _tdist = _odist = 0.0;
    _niter = 0;
    if(_diag > 0){
      countClusters(clusters, _nclu, _nchits, _tdist);
      _nhits = 0;
      for( auto const& ch : chcol )
	if((!_testflag) || (ch.flag().hasAllProperties(_sigmask) && !ch.flag().hasAnyProperty(_bkgmask)))++_nhits;
      _idiag->Fill();
    }
// iterate between assigning/creating clusters, updating and merging
    do {
      ++_niter;
      // assign hits to these clusters
      assignHits(clusters, chcol );
      // update the cluster positions for these new
      updateClusters(clusters, chcol);
      // check for merging
      _nmerge = mergeClusters(clusters, _dt, _dd2, chcol);
      // measure convergence
      _odist = _tdist;
      countClusters(clusters, _nclu, _nchits, _tdist);
      if(_diag > 0)
	_idiag->Fill();
    } while ( fabs(_odist - _tdist) > _maxdsum && _niter < _maxniter);
  // compress out the empty (merged) clusters from the list
    cleanClusters(clusters);
  }

  void TLTClusterer::initClusters( BkgClusterCollection& clusters,
	ComboHitCollection const& chcol) const {
// reset
    clusters.clear();
// dimension clusters proportional to the number of hits
    clusters.reserve( (unsigned)rint(chcol.size()/3.0) );
   // stereo init means take every stereo hit as a cluster seed, then merge
   // A more efficient implementation would loop over stereo hits themselves FIXME!
    if(_stereoinit){
      for(size_t ish = 0; ish < chcol.size(); ++ish) {
	ComboHit const& ch = chcol.at(ish);
	StrawHitFlag const& shf = ch.flag();
	if((!_testflag) || (shf.hasAllProperties(_stereo) && !shf.hasAnyProperty(_bkgmask))){
	  BkgCluster sclust(ch.pos(), ch.time());
	  sclust._hits.push_back(BkgClusterHit(0.0,ish,shf));
	  clusters.push_back(sclust);
	}
      }
      // merge with wide windows
      mergeClusters(clusters, _maxdt, _md2, chcol);
      // this process leaves lots of empties; flush them out to improve
      // efficiency later
      cleanClusters(clusters);
    } else {
      // otherwise, just take the 1st hit as a cluster and work from there
      for(size_t ish = 0; ish < chcol.size(); ++ish) {
	ComboHit const& ch = chcol.at(ish);
	StrawHitFlag const& shf = ch.flag();
	if((!_testflag) || (shf.hasAllProperties(_sigmask) && !shf.hasAnyProperty(_bkgmask))){
	  BkgCluster sclust(ch.pos(), ch.time());
	  sclust._hits.push_back(BkgClusterHit(0.0,ish,shf));
	  clusters.push_back(sclust);
	  break;
	}
      }
    }
  }

  void TLTClusterer::assignHits( BkgClusterCollection& clusters,
      ComboHitCollection const& chcol) const {
  // create a set of hits which are assigned to clusters
    set<uint16_t> assigned;
    // re-assign hits already assigned to the clusters.  This improves efficiency
    // loop over the clusters
    for(auto& clu : clusters) {
      if(clu.hits().size() > 0){
      // mark hits which are outside cuts for removal
	vector<size_t> toremove;
	toremove.reserve(clu.hits().size());
	// reverse-loop to facilitate removal
	for(int ihit = clu.hits().size()-1;ihit >= 0; --ihit){
	  if(clu.hits()[ihit].distance() > _dhit)
	    toremove.push_back(ihit);
	  else
	    assigned.insert(clu.hits()[ihit].index());
	}
	// actually remove the hits
	for(auto ihit : toremove){
	  std::swap(clu._hits[ihit],clu._hits.back());
	  clu._hits.pop_back();
	}
      }
    }
    // loop over the hits and find the closest cluster
    for(size_t ish = 0; ish < chcol.size(); ++ish) {
      ComboHit const& ch = chcol[ish];
      // select hits
      StrawHitFlag const& shf = ch.flag();
      if((!_testflag) || (shf.hasAllProperties(_sigmask) && !shf.hasAnyProperty(_bkgmask))){
	// skip hits already assigned
	auto ifind = assigned.find(ish);
	if(ifind == assigned.end()){
	  float mindist(FLT_MAX);
	  for(auto& clu : clusters) {
	    // loop over all non-empty clusters
	    if(clu.hits().size() > 0){
	      float dist = distance(clu,ch);
	      mindist = std::min(dist,mindist);
	      if(dist < _dhit){
		clu._hits.push_back(BkgClusterHit(dist,ish,shf));
		// assign to 1st cluster inside radius.  This doesn't correctly handle
		// overlapping clusters, but those should have been merged and this is more efficient
		break;
	      }
	    }
	  }
	  // hit isn't assigned and has minimal distance from other clusters, create a cluster from it
	  if(mindist > _dseed){
	    BkgCluster sclust(ch.pos(), ch.time());
	    sclust._hits.push_back(BkgClusterHit(0.0,ish,shf));
	    clusters.push_back(sclust);
	  }
	}
      }
    }
  }

  unsigned TLTClusterer::mergeClusters(BkgClusterCollection& clusters,
      float dt, float dd2,
      ComboHitCollection const& chcol) const {
// note: this algorithm merges pairwise iteratively. 
// More efficient could be to merge groups of clusters in a single iteration
    unsigned retval(0);
    unsigned niter(0);
    // struct to keep track of which clusters to merge
    vector<pair<size_t,size_t>> tomerge;
    do {
      tomerge.clear();
      for(size_t ic=0;ic < clusters.size(); ++ic){
	if(clusters[ic].hits().size()>0){
	  for(size_t jc = ic+1; jc < clusters.size(); ++jc){
	    if(clusters[jc].hits().size()>0){
	      if(fabs(clusters[ic].time() - clusters[jc].time()) < dt) {
		float d2 = (clusters[ic].pos() - clusters[jc].pos()).perp2();
		if(d2 < dd2){
		  tomerge.push_back(make_pair(ic,jc));
		}
	      }
	    }
	  }
	}
      }
      // acually merge
      for(auto& ipair : tomerge ){
	if(clusters[ipair.first].hits().size() > 0 &&
	    clusters[ipair.second].hits().size() > 0){
	  // merge smaller into bigger
	  if(clusters[ipair.first].hits().size() < clusters[ipair.second].hits().size()){
	    auto smaller = ipair.first;
	    ipair.first = ipair.second;
	    ipair.second = smaller;
	  }
	  // move the hits over from the smaller cluster
	  clusters[ipair.first]._hits.insert(clusters[ipair.first]._hits.end(),
	    clusters[ipair.second]._hits.begin(),clusters[ipair.second]._hits.end());
	  clusters[ipair.second]._hits.clear();
	  // re-compute the position for the merged cluster
	  updateCluster(clusters[ipair.first], chcol);
	}
      }
      retval += tomerge.size();
      ++niter;
    } while(tomerge.size() > 0 && niter < _maxniter);
    return retval;
  }

// update local cache (position, time)
  void TLTClusterer::updateClusters(BkgClusterCollection& clusters,
      ComboHitCollection const& chcol) const {
    for (auto& clu : clusters) {
      if(clu.hits().size() > 0)
	updateCluster(clu,chcol);
    }
  }

  void TLTClusterer::updateCluster(BkgCluster& cluster,
      ComboHitCollection const& chcol) const {
    // accumulators for position calculation.  For now just median,
    // should implement full 2-d finding including different errors FIXME!
    if(_palg == median) {
      accumulator_set<float, stats<tag::weighted_median(with_p_square_quantile) >, unsigned > racc, pacc, tacc;
      // work relative to the current cluster position
      // calculate cluster info
      float crho = sqrtf(cluster.pos().perp2());
      float cphi = cluster.pos().phi();
      // choose local cylindrical coordinates WRT the cluster; these map roughly onto the straw directions
//      XYZVec rdir = PerpVector(cluster.pos(),Geom::ZDir()).unit();
//      XYZVec pdir(-rdir.y(),rdir.x(),0.0);
      for(auto const& chit : cluster.hits()) {
	ComboHit const& ch = chcol.at(chit.index());
	float dt = ch.time() - cluster.time();
	float dr = sqrtf(ch.pos().perp2()) - crho;
	float phi = ch.pos().phi();
	float dp = Angles::deltaPhi(phi,cphi);
	// weight according to the # of hits
	// project weight along radial and azimuthal airections
	tacc(dt,weight=ch.nStrawHits());
	racc(dr,weight=ch.nStrawHits());
	pacc(dp,weight=ch.nStrawHits());
      }
      float dt = extract_result<tag::weighted_median>(tacc);
      float dr = extract_result<tag::weighted_median>(racc);
      float dp = extract_result<tag::weighted_median>(pacc);
      cluster._time += dt;
      crho += dr;
      cphi += dp;
      cluster._pos = XYZVec(crho*cos(cphi),crho*sin(cphi),0.0);
    } else {
      throw cet::exception("RECO")<<"mu2e::TLTClusterer: algorithm not implemented"<< endl;
    }
    // update hit distances
    for(auto& hit : cluster._hits)
      hit._dist = distance(cluster,chcol[hit.index()]);
  }

  void TLTClusterer::cleanClusters(BkgClusterCollection& clusters) const {
    auto iclu = clusters.begin();
    while(iclu < clusters.end()) {
      if(iclu->hits().size() ==0){
	// overwrite this element with the last and pop it off the vector
	// don't advance as we still need to check this element
	*iclu = clusters.back();
	clusters.pop_back();
      } else
	++iclu;
    }
  }

  void TLTClusterer::countClusters(BkgClusterCollection const& clusters,
      unsigned& nclu, unsigned& nchits, float& tdist) const {
    nclu = nchits = 0;
    tdist = 0.0;
    for (auto const& clu : clusters) {
      if(clu.hits().size() > 0){
	++nclu;
	nchits += clu.hits().size();
	// compute the total distance for all hits
	for(auto const& chit : clu.hits() )
	  tdist += chit.distance();
      }
    }
  }

} // mu2e namespace
