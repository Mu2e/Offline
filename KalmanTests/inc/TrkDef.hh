//
// Define a track; this provides the transfer between pat. rec. and fitting
//
// $Id: TrkDef.hh,v 1.1 2011/06/02 00:00:05 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/02 00:00:05 $
//
// Original author David Brown, LBNL
//
#ifndef TrkDef_HH
#define TrkDef_HH
// Mu2e
#include "ToyDP/inc/StrawHitCollection.hh"
// BaBar
#include "TrkBase/HelixTraj.hh"
// CLHEP
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

namespace mu2e 
{
  class TrkDef {
  public:
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<size_t>& strawhits,
      const HelixTraj& helix, double t0, double t0err);
    TrkDef(const StrawHitCollection* strawcollection, const std::vector<size_t>& strawhits,
      const HepVector& parvec, const HepSymMatrix& covar, double t0, double t0err );
    TrkDef();
    TrkDef(const TrkDef&);
    TrkDef& operator = (const TrkDef&);
    ~TrkDef();
  // append a straw hit to this track definition
    void appendHit(size_t index) { _indices.push_back(index); }
  // accessors
    const StrawHitCollection* strawHitCollection() const { return _straws; }
    const std::vector<size_t>& strawHitIndices() const { return _indices;}
    const HelixTraj& helix() const { return _h0; }
    double t0() const { return _t0; }
    double t0Err() const { return _t0err; }
    void setT0(double t0, double t0err) { _t0 = t0; _t0err = t0err; }
  private:
    const StrawHitCollection* _straws; // straw hit collection
    std::vector<size_t> _indices; // indices to straw hits in the collection to use for this track
    HelixTraj _h0; // helix estimate, valid in the region around z=0
    double _t0; // initial estimate of track t0, defined when the track reaches z=0;
    double _t0err; // initial estimate of t0 error
    // dummy variables
    static HepVector _dpar;
    static HepSymMatrix _dcov;
  };
}

#endif
