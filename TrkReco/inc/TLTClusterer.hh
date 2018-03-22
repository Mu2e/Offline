//
// Object to cluster straw hits, used in background removal and track fitting.  This
// implementation uses the 2-level distance scheme, where hits are associated if their
// distance to an existing cluster is below the first threshold, and are used to seed
// a new cluster if their distance to every other existing cluster is above the 2nd threshold
//
//  David Brown (LBNL) 2013
//
#ifndef TLTClusterer_HH
#define TLTClusterer_HH

// base class
#include "TrkReco/inc/BkgClusterer.hh"
// framework
#include "fhiclcpp/ParameterSet.h"
// root
#include "TTree.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/mean.hpp>

namespace mu2e 
{

  class TLTClusterer : public BkgClusterer
  {
  public:
  enum PosAlgorithm { median=0, mean}; // algorithm for comuting cluster position
// parameter set should be passed in on construction
  explicit TLTClusterer(fhicl::ParameterSet const&);
  virtual ~TLTClusterer();
  virtual void init();
  //  base class function: given the straw hits and associated data, find the clusters.
  virtual void findClusters(BkgClusterCollection& clusters,
      ComboHitCollection const& chcol);
  private:
  // functions to define the distance between a hit and cluster
  // this is the variable against which the clustering thresholds are applied
  float distance(BkgCluster const&, ComboHit const& ch) const; 
  // merge overlapping clusters
  unsigned mergeClusters(BkgClusterCollection& bkgcol,
      float dt, float drho2,
      ComboHitCollection const& chcol ) const;
 // compress out empties
  void cleanClusters(BkgClusterCollection& bkgcol ) const;
  // initialize clusters.  Different strategeis can be tried
  void initClusters( BkgClusterCollection& clusters, ComboHitCollection const& chcol) const;
  // assign hits to the given clusters.  Hits already assigned
  // are tested first, for efficiency
  void assignHits(BkgClusterCollection& bkgcol, ComboHitCollection const& chcol) const;
  // compute the total distance between hits and clusters
  float hitDistance(BkgClusterCollection const& bkgcol, ComboHitCollection const& chcol) const;
  // update the cluster state (position, distances)
  void updateClusters(BkgClusterCollection& clusters, ComboHitCollection const& chcol) const;
  //
  void updateCluster(BkgCluster& cluster, ComboHitCollection const& chcol) const;
  //
  void countClusters(BkgClusterCollection const& clusters, unsigned& nclu, unsigned& nchits, float& tdist) const;
  void refineClusters(BkgClusterCollection& clusters, ComboHitCollection const& chcol) const;
  void refineCluster(BkgCluster& cluster, ComboHitCollection const& chcol ) const;

  // configuration parameters
  int _diag,_debug;
  PosAlgorithm _palg;
  bool _testflag;
  StrawHitFlag _bkgmask; // mask for hits to ignore
  StrawHitFlag _sigmask; // mask for hits to use
  StrawHitFlag _stereo; // mask for selecting stereo hits
  float _dseed; // Minimum separation to seed a new cluster
  float _dhit; // Maximum separation to include a hit in a cluster
  float _dd; // cluster diameter
  float _dd2; // cluster diameter squared, cached for efficiency
  float _dt; // natural time spread
  float _maxdt; // maximum time difference
  float _trms2; // time RMS
  float _md2; // // cached square of maximum distance
  float _maxdsum; // maximum total distance change to consider 'converged'
  unsigned _maxniter; // maximum number of iterations
  bool _stereoinit; // initialize using only stereo hits
  float _maxwt; // maximum weight to give a single hit
  float _minerr; // minimum error to give a single hit

  // diagnostics
  mutable TTree* _idiag; // iteration diagnostics
  mutable unsigned _niter, _nmerge, _nclu, _nhits, _nchits;
  mutable float _odist, _tdist;
  };
}
#endif
