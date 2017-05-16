//
// Object to perform helix fit to straw hits
//
#include "TrkReco/inc/TLTClusterer.hh"
// services
#include "art/Framework/Services/Optional/TFileService.h"
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
// C++
#include <vector>
#include <set>
#include <algorithm>
using namespace boost::accumulators;
using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e
{
// struct for hit MVA
  struct BkgHitMVA {
    std::vector<Double_t> _pars;
    Double_t& _rho; // transverse radius of this hit WRT cluster center
    Double_t& _drho; // rho diff to cluster median radius
    Double_t& _crho; // normalized rho difference
    Double_t& _dt;  // time diff to cluster
    BkgHitMVA() : _pars(4,0.0),_rho(_pars[0]),_drho(_pars[1]),_crho(_pars[2]),_dt(_pars[3]){}
  };


 TLTClusterer::TLTClusterer(fhicl::ParameterSet const& pset) :
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _palg(static_cast<PosAlgorithm>(pset.get<int>("ClusterPositionAlgorith",median))),
    _bkgmask(pset.get<vector<string> >("BackgroundMask",vector<string>())),
    _sigmask(pset.get<vector<string> >("SignalMask",vector<string>())),
    _dseed(pset.get<double>("SeedDistance",100.0)), // minimum 'chisquared' to define a new cluster
    _dhit(pset.get<unsigned>("HitDistance",25.0)), // maximum 'chisquared' to add hit to a cluster
    _dmerge(pset.get<double>("MergeDistance",9.0)), // minimum 'chisquared' to merge clusters
    _dd(pset.get<double>("ClusterDiameter",10.0)), // the natural cluster transverse size (mm)
    _dt(pset.get<double>("TimeDifference",30.0)), // the natural time spread (nsec)
    _maxdt(pset.get<double>("MaxTimeDifference",50.0)), // Maximum time difference (nsec)
    _maxniter(pset.get<unsigned>("MaxNIterations",50)),
    _maxdist(pset.get<double>("MaxDistance",10.0)), // iteration convergence
    _refine(pset.get<bool>("RefineClusters",false)),
    _minnrefine(pset.get<unsigned>("MinNRefineHits",5)), // # of hits before refining the cluster
    _minerr(pset.get<double>("MinHitError",5.0)), // corresponds to an error of the straw diameter
    _minmva(pset.get<double>("MinHitMVA",0.2)), // used with MVA hit refining
    _hitMVA(pset.get<fhicl::ParameterSet>("HitMVA",fhicl::ParameterSet())),
    _idiag(0)
  {
    // cache some values to avoid repeating FP operations on fixed numbers
    double trms(pset.get<double>("TimeRMS",1.0)); // effective individual hit time resolution (nsec)
    _trms2 = trms*trms;
    double maxdist(pset.get<double>("MaxDistance",50.0)); // Maximum transverse distance (mm)
    _md2 = maxdist*maxdist;
    _dd2 = _dd*_dd;
    _stereo = StrawHitFlag(StrawHitFlag::stereo);
    _stereo.merge(_sigmask);
    _maxwt = 1.0/_minerr;
 }

  TLTClusterer::~TLTClusterer()
  {}

  double TLTClusterer::distance(BkgCluster const& cluster, StrawHit const& sh,
    StrawHitPosition const& shp) const {
    double retval = _dseed+1.0; // default return is above seed threshold
    // compute the simplest parts first to avoid expensive calculations on distant hits
    double dt = fabs(sh.time()-cluster.time());
    if( dt < _maxdt){
// compute spatial distance
      Hep3Vector psep = (shp.pos()-cluster.pos()).perpPart();
      // work in squared magnitudes for efficiency
      double d2 = psep.mag2();
      if( d2 < _md2) {
	retval = 0.0; // if separation is inside natural size, count distance as 0
	// only count differences above the natural hit size (drift time, cluster size)
	if(dt > _dt){
	  double tdist = dt -_dt;
	  retval += tdist*tdist/_trms2;
	}
	if(d2 > _dd2){
	  // project separation along wire and transverse directions
	  double dw = shp.wdir().dot(psep);
	  Hep3Vector that(-shp.wdir().y(),shp.wdir().x(),0.0);
	  double dp = that.dot(psep);
	// add these contributions
	  retval += (dw*dw)/(shp.posRes(StrawHitPosition::wire)*shp.posRes(StrawHitPosition::wire)) +
	    (dp*dp)/(shp.posRes(StrawHitPosition::trans)*shp.posRes(StrawHitPosition::trans));
	}
      }
    }
    return retval;
  }

  void TLTClusterer::init() {
    // initialize the MVAs
    _hitMVA.initMVA();
    if(_debug > 0){
      cout << "Hit MVA : " << endl;
      _hitMVA.showMVA();
    }
    // setup diagnostics
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _idiag = tfs->make<TTree>("idiag","iteration diagnostics");
      _idiag->Branch("niter",&_niter,"niter/U");
      _idiag->Branch("nmerge",&_nmerge,"nmerge/U");
      _idiag->Branch("nclu",&_nclu,"nclu/U");
      _idiag->Branch("nchits",&_nchits,"nchits/U");
      _idiag->Branch("nhits",&_nhits,"nhits/U");
      _idiag->Branch("odist",&_odist,"odist/D");
      _idiag->Branch("tdist",&_tdist,"tdist/D");
    }
  }

  void TLTClusterer::findClusters(BkgClusterCollection& clusters,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol ) const {
// require consistency
    if(shcol.size() != shpcol.size() || shcol.size() != shfcol.size()){
      ostringstream os;
      os <<  " TLTClusterer: inconsistent collection lengths ";
      throw out_of_range( os.str() );
    }
// initialize the clusters
    initClusters(clusters, shcol, shpcol, shfcol );
// compute the initial total distance
    _tdist = _odist = 0.0;
// iterate between assigning/creating clusters, updating and merging
    _niter = 0;
    do {
      ++_niter;
      // assign hits to these clusters
      assignHits(clusters, shcol, shpcol, shfcol );
      // update the cluster positions for these new
      updateClusters(clusters, shcol, shpcol);
      // check for merging
      _nmerge = mergeClusters(clusters, shcol, shpcol);
      // measure convergence
      _odist = _tdist;
      countClusters(clusters, _nclu, _nchits, _tdist);
      if(_diag > 0){
// record convergence information: nhits, distance, iteration
	_nhits = 0;
	for( auto const& shf : shfcol )
	  if(shf.hasAllProperties(_stereo) && !shf.hasAnyProperty(_bkgmask))++_nhits;
	_idiag->Fill();
      }
    } while ( fabs(_odist - _tdist) > _maxdist && _niter < _maxniter);
  // optionally refine the hit assignment using an MVA
    if(_refine)refineClusters(clusters,shcol, shpcol);
  // compress out the empty (merged) clusters from the list
    cleanClusters(clusters);
  }

  void TLTClusterer::initClusters( BkgClusterCollection& clusters,
	StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol) const {
// reset
    clusters.clear();
// dimension clusters proportional to the number of hits
    clusters.reserve( (unsigned)rint(shcol.size()/3.0) );
   // stereo init means take every stereo hit as a cluster seed, then merge
    if(_stereoinit){
      for(size_t ish = 0; ish < shcol.size(); ++ish) {
	StrawHitFlag const& shf = shfcol.at(ish);
	if(shf.hasAllProperties(_stereo) && !shf.hasAnyProperty(_bkgmask)){
	  StrawHit const& sh = shcol.at(ish);
	  StrawHitPosition const& shp = shpcol.at(ish);
	  BkgCluster sclust(shp.pos(), sh.time());
	  sclust._hits.push_back(BkgClusterHit(0.0,ish,shf));
	  clusters.push_back(sclust);
	}
      }
      mergeClusters(clusters, shcol, shpcol);
    } else {
      // otherwise, just take the 1st hit as a cluster and work from there
      for(size_t ish = 0; ish < shcol.size(); ++ish) {
	StrawHitFlag const& shf = shfcol.at(ish);
	if(shf.hasAllProperties(_sigmask) && !shf.hasAnyProperty(_bkgmask)){
	  StrawHit const& sh = shcol.at(ish);
	  StrawHitPosition const& shp = shpcol.at(ish);
	  BkgCluster sclust(shp.pos(), sh.time());
	  sclust._hits.push_back(BkgClusterHit(0.0,ish,shf));
	  clusters.push_back(sclust);
	  break;
	}
      }
    }
  }

  void TLTClusterer::assignHits( BkgClusterCollection& clusters,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol,
      StrawHitFlagCollection const& shfcol) const {
  // create a set of hits which are assigned to clusters
    set<StrawHitIndex> assigned;
    // re-assign hits already assigned to the clusters.  This improves efficiency
    // loop over the clusters
    for(auto& clu : clusters) {
      if(clu.hits().size() > 0){
	// remember which hits were in this cluster
	vector<BkgClusterHit> hits(clu.hits());
	// remove the hits
	clu._hits.clear();
	// loop over hits that were in the cluster
	for(auto& hit : hits) {
	  if(hit.distance() < _dhit){
	    assigned.insert(hit.index());
	    clu._hits.push_back(hit);
	  }
	}
      }
    }
    // loop over the hits and find the closest cluster
    for(size_t ish = 0; ish < shcol.size(); ++ish) {
      // select hits
      StrawHitFlag const& shf = shfcol.at(ish);
      if(shf.hasAllProperties(_sigmask) && !shf.hasAnyProperty(_bkgmask)){
	// skip hits already assigned
	auto ifind = assigned.find(ish);
	if(ifind == assigned.end()){
	  double mindist(FLT_MAX);
	  for(auto& clu : clusters) {
	    // loop over all non-empty clusters
	    if(clu.hits().size() > 0){
	      double dist = distance(clu,shcol.at(ish), shpcol.at(ish));
	      if(dist < _dhit){
		clu._hits.push_back(BkgClusterHit(dist,ish,shf));
		// assign to 1st cluster inside radius.  This doesn't correctly handle
		// overlapping clusters, but those should have been merged and this is more efficient
		break;
	      }
	      mindist = std::min(dist,mindist);
	    }
	  }
	  // hit isn't assigned and has minimal distance from other clusters, create a cluster from it
	  if(mindist > _dseed){
	    BkgCluster sclust(shpcol.at(ish).pos(), shcol.at(ish).time());
	    sclust._hits.push_back(BkgClusterHit(0.0,ish,shf));
	    clusters.push_back(sclust);
	  }
	}
      }
    }
  }

  unsigned TLTClusterer::mergeClusters(BkgClusterCollection& clusters,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol) const {
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
	      if(fabs(clusters[ic].time() - clusters[jc].time()) < _dt) {
		double d2 = (clusters[ic].pos() - clusters[jc].pos()).perp2();
		if(d2 < _dd2){
		  tomerge.push_back(make_pair(ic,jc));
		  break;
		}
	      }
	    }
	  }
	}
      }
      // acually merge
      for(auto const& ipair : tomerge )
	mergeClusters(clusters[ipair.first],clusters[ipair.second], shcol, shpcol);
      retval += tomerge.size();
      ++niter;
    } while(tomerge.size() > 0 && niter < _maxniter);
    return retval;
  }

  // merge 2 clusters by absorbing the smaller into the larger.  This leaves
  // a 'ghost' cluster which is ignored during processing and is compressed
  // out after all clustering is done
  void TLTClusterer::mergeClusters(BkgCluster& c1, BkgCluster& c2,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol) const {
    if(c1.hits().size() > c2.hits().size()){
      c1._hits.insert(c1._hits.end(),c2._hits.begin(),c2._hits.end());
      c2._hits.clear();
      // re-compute the position for the merged cluster
      if(c1.hits().size() > 1)updateCluster(c1, shcol, shpcol);
    } else {
      c2._hits.insert(c2._hits.end(),c1._hits.begin(),c1._hits.end());
      c1._hits.clear();
      if(c2.hits().size() > 1)updateCluster(c2, shcol, shpcol);
    }
  }

// update local cache (position, time)
  void TLTClusterer::updateClusters(BkgClusterCollection& clusters,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol) const {
    for (auto& clu : clusters) {
      if(clu.hits().size() > 0)
	updateCluster(clu,shcol, shpcol);
    }
  }

  void TLTClusterer::updateCluster(BkgCluster& cluster,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol) const {
    // accumulators for position calculation.  For now just median,
    // should implement full 2-d finding including different errors FIXME!
    if(_palg == median) {
      accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > racc, pacc, tacc;
      // work relative to the current cluster position
      // calculate cluster info
      double crho = cluster.pos().perp();
      double cphi = cluster.pos().phi();
      // choose local cylindrical coordinates WRT the cluster; these map reasonably onto the straw directions
      Hep3Vector rdir = cluster.pos().perpPart().unit();
      Hep3Vector pdir(-rdir.y(),rdir.x(),0.0);
      for(auto const& chit : cluster.hits()) {
	StrawHit const& sh = shcol.at(chit.index());
	StrawHitPosition const& shp = shpcol.at(chit.index());
	double dt = sh.time();
	double dr = shp.pos().perp() - crho;
	double phi = shp.pos().phi();
	double dp = Angles::deltaPhi(phi,cphi);
	// only the error along the straw contributes significantly
	double dw = shp.posRes(StrawHitPosition::wire);
	Hep3Vector const& wdir = shp.wdir();
	// weight according to the wire direction; linearly for now
	double twt = 1.0/dw;
	double rwt = fabs(rdir.dot(wdir))*twt;
	double pwt = fabs(pdir.dot(wdir))*twt;
	tacc(dt,weight=twt);
	racc(dr,weight=rwt);
	pacc(dp,weight=pwt);
      }
      double dt = extract_result<tag::weighted_median>(tacc);
      double dr = extract_result<tag::weighted_median>(racc);
      double dp = extract_result<tag::weighted_median>(pacc);
      cluster._time += dt;
      crho += dr;
      cphi += dp;
      cluster._pos = Hep3Vector(crho*cos(cphi),crho*sin(cphi),0.0);
    } else {
      throw cet::exception("RECO")<<"mu2e::TLTClusterer: algorithm not implemented"<< endl;
    }
    // update hit distances
    for(auto& hit : cluster._hits)
      hit._dist = distance(cluster,shcol.at(hit.index()), shpcol.at(hit.index()));
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
      unsigned& nclu, unsigned& nchits, double& tdist) const {
    nclu = nchits = 0;
    tdist = 0.0;
    for (auto const& clu : clusters) {
      if(clu.hits().size() > 0){
	++nclu;
	nchits += clu.hits().size();
	for(auto const& chit : clu.hits() )
	  tdist += chit.distance();
      }
    }
  }

  void TLTClusterer::refineClusters(BkgClusterCollection& clusters,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol) const {
    for( auto& cluster: clusters) {
      if(cluster.hits().size() >= _minnrefine)
	refineCluster(cluster,shcol, shpcol);
    }
  }

  // use an MVA to refine the cluster
  void TLTClusterer::refineCluster(BkgCluster& cluster,
      StrawHitCollection const& shcol,
      StrawHitPositionCollection const& shpcol) const {
      // determine a radius for the hits around the center of this cluster
    accumulator_set<double, stats<tag::weighted_median(with_p_square_quantile) >, double > racc;
    for(auto const& chit : cluster.hits()){
      StrawHitPosition const& shp = shpcol.at(chit.index());

      Hep3Vector psep = (shp.pos()-cluster.pos()).perpPart();
      double rho = psep.mag();
      double rwt = _maxwt;
      Hep3Vector pdir = psep.unit();
      double dr = shp.posRes(StrawHitPosition::wire)*pdir.dot(shp.wdir());
      if(dr > _minerr)rwt = 1.0/dr;
      racc(rho,weight=rwt);
    }
    cluster._rho = extract_result<tag::weighted_median>(racc);
    // now compute the MVA for each hit
    for(auto& chit : cluster._hits){
      StrawHit const& sh = shcol.at(chit.index());
      StrawHitPosition const& shp = shpcol.at(chit.index());

      Hep3Vector psep = (shp.pos()-cluster.pos()).perpPart();
      double rho = psep.mag();
      Hep3Vector pdir = psep.unit();
      double rerr = std::max(_minerr,shp.posRes(StrawHitPosition::wire)*pdir.dot(shp.wdir()));
      double dr = rho - cluster.rho();
      // fill the struct
      BkgHitMVA hmva;
      hmva._rho = rho;
      hmva._drho = dr;
      hmva._crho = dr/rerr;
      hmva._dt = sh.time() - cluster.time();
      // MVA value OVERWRITES the geometric + time distance
      chit._dist = _hitMVA.evalMVA(hmva._pars);
      // if the hit MVA value is too small, flag off the hit
      if(chit._dist < _minmva){ 
	chit._flag.clear(StrawHitFlag::active);
      }
    }
    // set cluster bit
    cluster._flag.merge(BkgClusterFlag::refined);
  }
} // mu2e namespace
