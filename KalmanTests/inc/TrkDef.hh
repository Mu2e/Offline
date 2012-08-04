//
// Define a track; this provides the transfer between pat. rec. and fitting
//
// $Id: TrkDef.hh,v 1.13 2012/08/04 00:38:06 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/08/04 00:38:06 $
//
// Original author David Brown, LBNL
//
#ifndef TrkDef_HH
#define TrkDef_HH
// Mu2e
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
// BaBar
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkT0.hh"
#include "TrkBase/TrkParticle.hh"
// CLHEP
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

class TrkDifPieceTraj;

namespace mu2e 
{

  struct hitIndex {
    size_t _index; // index into the straw hit container
    int _ambig; // hit ambiguity.  0 means compute from track
    hitIndex() : _index(0),_ambig(0) {}
    hitIndex(size_t index,int ambig=0) : _index(index),_ambig(ambig) {}
    hitIndex& operator = (size_t index) { _index = index; return *this; }
  };
  
  class TrkDef {
  public:
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      const HelixTraj& helix, TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream,
      double t0=0.0, double t0err=-1.0);
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      const HepVector& parvec, const HepSymMatrix& covar,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream,
      double t0=0.0, double t0err=-1.0 );
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const StrawHitCollection* strawcollection,
    TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(TrkParticle const& tpart=_eminus, TrkFitDirection const& fdir=_downstream);
    TrkDef(const TrkDef&);
    TrkDef& operator = (const TrkDef&);
    ~TrkDef();
  // append a straw hit to this track definition
    void appendHit(size_t index,int ambig=0) { _indices.push_back(hitIndex(index,ambig)); }
  // accessors
    const StrawHitCollection* strawHitCollection() const { return _straws; }
    const std::vector<hitIndex>& strawHitIndices() const { return _indices;}
    const HelixTraj& helix() const { return _h0; }
    const TrkDifPieceTraj* traj() const { return _ptraj; }
    const TrkT0& trkT0() const { return _t0; }
    TrkParticle const& particle() const { return _tpart; }
    TrkFitDirection const& fitdir() const { return _fdir; }
    void setHelix(HelixTraj const& helix) { _h0 = helix; }
    void setTraj(const TrkDifPieceTraj* ptraj) { _ptraj = ptraj; }
    void setTrkT0(double t0, double t0err) { _t0.setT0(t0,t0err); }
    void setTrkT0(const TrkT0& t0) { _t0 = t0; }
    void setIndices(std::vector<hitIndex> const& indices ) { _indices = indices; }
    unsigned eventId() const { return _eventid; }
    unsigned trackId() const { return _trkid; }
    void setEventId(unsigned eventid) { _eventid = eventid; }
    void setTrackId(unsigned trkid) { _trkid = trkid; }
  private:
    unsigned _eventid;
    unsigned _trkid;
    const StrawHitCollection* _straws; // straw hit collection
    std::vector<hitIndex> _indices; // indices to straw hits in the collection to use for this track
    HelixTraj _h0; // helix estimate, valid in the region around z=0
    const TrkDifPieceTraj* _ptraj; // optional initial estimate of the trajectory
// particle type.  Note this defines both the charge and the mass
    TrkParticle _tpart;
// fit direction
    TrkFitDirection _fdir;
    TrkT0 _t0; // t0 estimate
    // dummy variables
    static HepVector _dpar;
    static HepSymMatrix _dcov;
    static TrkParticle _eminus;
    static TrkFitDirection _downstream;
  };
}

#endif
