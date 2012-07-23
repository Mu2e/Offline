//
// Define a track; this provides the transfer between pat. rec. and fitting
//
// $Id: TrkDef.hh,v 1.11 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/23 17:52:27 $
//
// Original author David Brown, LBNL
//
#ifndef TrkDef_HH
#define TrkDef_HH
// Mu2e
#include "RecoDataProducts/inc/StrawHitCollection.hh"
// BaBar
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkT0.hh"
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
// define the track direction: either downstream (dZ/dt > 0) or upstream (dZ/dt < 0);
    enum trkdir {downstream=0,upstream};
    
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      const HelixTraj& helix, double t0=0.0, double t0err=-1.0);
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits,
      const HepVector& parvec, const HepSymMatrix& covar, double t0=0.0, double t0err=-1.0 );
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<hitIndex>& strawhits);
    TrkDef(const StrawHitCollection* strawcollection);
    TrkDef();
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
    void setHelix(HelixTraj const& helix) { _h0 = helix; }
    void setTraj(const TrkDifPieceTraj* ptraj) { _ptraj = ptraj; }
    void setTrkT0(double t0, double t0err) { _t0.setT0(t0,t0err); }
    void setTrkT0(const TrkT0& t0) { _t0 = t0; }
    void setIndices(std::vector<hitIndex> const& indices ) { _indices = indices; }
  private:
    const StrawHitCollection* _straws; // straw hit collection
    std::vector<hitIndex> _indices; // indices to straw hits in the collection to use for this track
    HelixTraj _h0; // helix estimate, valid in the region around z=0
    const TrkDifPieceTraj* _ptraj; // optional initial estimate of the trajectory
    TrkT0 _t0; // t0 estimate
// fit direction
    trkdir _fitdir;
// particle type.  Note this defines both the charge and the mass
  
    // dummy variables
    static HepVector _dpar;
    static HepSymMatrix _dcov;
  };
}

#endif
