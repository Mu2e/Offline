//
// Define a track; this provides the transfer between pat. rec. and fitting
//
// $Id: TrkDef.hh,v 1.18 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2014/08/30 12:19:38 $
//
// Original author David Brown, LBNL
//
#ifndef TrkDef_HH
#define TrkDef_HH
// Mu2e
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
// BaBar
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

  class TrkDef {
  public:
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      const HelixTraj& helix, TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      const CLHEP::HepVector& parvec, const CLHEP::HepSymMatrix& covar,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const StrawHitCollection* strawcollection,TrackSeed const& seed,
    TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const StrawHitCollection* strawcollection,
    TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const TrkDef&);
    TrkDef& operator = (const TrkDef&);
    ~TrkDef();
  // append a straw hit to this track definition
    void appendHit(size_t index,int ambig=0) { _seed._timeCluster._strawHitIdxs.push_back(hitIndex(index,ambig)); }
  // accessors
    const StrawHitCollection* strawHitCollection() const { return _straws; }
    const std::vector<hitIndex>& strawHitIndices() const { return _seed._timeCluster._strawHitIdxs;}
    const HelixTraj& helix() const { return _h0; }
    const CLHEP::HepSymMatrix &helixCovMatr() const { return _h0.parameters()->covariance(); }
    const CLHEP::HepSymMatrix &covMatr() const { return _dcov; }
    TrkParticle const& particle() const { return _tpart; }
    TrkFitDirection const& fitdir() const { return _fdir; }
    void setHelix(HelixTraj const& helix) { _h0 = helix; }
    void setIndices(std::vector<hitIndex> const& indices ) { _seed._timeCluster._strawHitIdxs = indices; }
    TrkT0 const& t0() const { return _seed._timeCluster._t0; }
    void setT0( TrkT0 const& t0) { _seed._timeCluster._t0 = t0; }
  protected:
    const StrawHitCollection* _straws; // straw hit collection
    TrackSeed _seed; // track seed
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
