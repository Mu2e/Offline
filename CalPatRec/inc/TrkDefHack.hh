//
// Define a track; this provides the transfer between pat. rec. and fitting
//
// $Id: TrkDefHack.hh,v 1.18 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2014/08/30 12:19:38 $
//
// Original author David Brown, LBNL
//
#ifndef TrkDefHack_HH
#define TrkDefHack_HH
// Mu2e
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
// BTrk includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkT0.hh"
// CLHEP
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

class TrkDifPieceTraj;

namespace mu2e 
{
  //  typedef size_t hitIndex;
  class TrkDefHack {
  public:
    TrkDefHack(TimeCluster const& tcluster, HelixTraj const& helix,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    // TrkDefHack(const StrawHitCollection* shcol, std::vector<mu2e::hitIndex> const& hits, 
    //   HelixTraj const& htraj,
    //   TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    // TrkDefHack(const StrawHitCollection* shcol, std::vector<mu2e::hitIndex> const& hits, 
    //   TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDefHack(const StrawHitCollection* shcol, std::vector<StrawHitIndex> const& hits, 
      HelixTraj const& htraj,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDefHack(const StrawHitCollection* shcol, std::vector<StrawHitIndex> const& hits, 
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDefHack(TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDefHack(const TrkDefHack&);
    TrkDefHack& operator = (const TrkDefHack&);
    ~TrkDefHack();
  // append a straw hit to this track definition
    void appendHit(size_t index) { _timeCluster._strawHitIdxs.push_back(index); }
  // accessors
    //    const std::vector<hitIndex>& strawHitIndices() const { return _timeCluster._strawHitIdxs;}
    const std::vector<StrawHitIndex>& strawHitIndices() const { return _timeCluster._strawHitIdxs;}
    const HelixTraj& helix() const { return _h0; }
    const CLHEP::HepSymMatrix &helixCovMatr() const { return _h0.parameters()->covariance(); }
    const CLHEP::HepSymMatrix &covMatr() const { return _dcov; }
    TrkParticle const& particle() const { return _tpart; }
    TrkFitDirection const& fitdir() const { return _fdir; }
    void setHelix(HelixTraj const& helix) { _h0 = helix; }
    //    void setIndices(std::vector<hitIndex> const& indices ) { _timeCluster._strawHitIdxs = indices; }
    void setIndices(std::vector<StrawHitIndex> const& indices ) { _timeCluster._strawHitIdxs = indices; }
    TrkT0 const& t0() const { return _timeCluster._t0; }
    void setT0( TrkT0 const& t0) { _timeCluster._t0 = t0; }
    const StrawHitCollection* strawHitCollection() const { return _shcol ; }
  protected:
    const StrawHitCollection* _shcol;
    TimeCluster _timeCluster; // t0 and hit indices
    HelixTraj _h0; // helix estimate, valid in the region around z=0
// particle type.  Note this defines both the charge and the mass
    TrkParticle _tpart;
// fit direction
    TrkFitDirection _fdir;
// dummy variables
    static CLHEP::HepVector _dpar;
    static CLHEP::HepSymMatrix _dcov;
    static TrkParticle _eminus;
    static TrkFitDirection _downstream;
  };
}

#endif
