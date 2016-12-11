
// BaBar hit object corresponding to a single straw hit
//
// Original author David Brown, LBNL
//
#ifndef TrkStrawHit_HH
#define TrkStrawHit_HH
// BaBar
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "RecoDataProducts/inc/HitT0.hh"
// Mu2e
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/HitT0.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// c++
#include <vector>
#include <functional>
// forward refs
class TrkDifTraj;

namespace mu2e
{
  class TrkStrawHit : public TrkHit {
  public:
  // enum for hit flags
    enum TrkStrawHitFlag {weededHit=-5, driftFail=-3, updateFail=-1,addedHit=3,unweededHit=4};
    TrkStrawHit(const StrawHit& strawhit, const Straw& straw,unsigned istraw,
    const HitT0& trkt0, double fltlen, double exterr, double maxdriftpull);
    virtual ~TrkStrawHit();
//  implementation of TrkHit interface
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }
    virtual int ambig() const { return _iamb; }
//    virtual void invert();
    virtual void setAmbig(int newambig);
    void setAmbigUpdate(bool update) { _ambigupdate = update; }
    unsigned index() const { return _istraw; } // index into StrawHit vector
    double hitRMS() const { return _t2d._rdrifterr;}
// strawhit specific interface
    const StrawHit& strawHit() const { return _strawhit; }
    const Straw& straw() const { return _straw; }
// correct the hit time for the wire propagation
    double time() const;
    double driftRadius() const { return _t2d._rdrift;}
    double driftRadiusErr() const { return _t2d._rdrifterr;}
    double driftVelocity() const { return _t2d._vdrift; }
    double timeDiffDist() const { return _tddist; }
    double timeDiffDistErr() const { return _tddist_err; }
    const CLHEP::Hep3Vector& wirePosition() const { return _wpos; }
    const CLHEP::Hep3Vector& wirePositionError() const { return _wpos_err; }
    void hitPosition(CLHEP::Hep3Vector& hpos) const;
    HitT0 const& hitT0() const { return _hitt0;}
    void updateHitT0(HitT0 const& t0) { _hitt0 = t0; }
    double signalTime() const { return _stime; } // time for signal to reach the end of the wire
// external hit error (mm); the intrinsic error comes from the t2d calibration object
    double extErr() const { return _exterr; }
// error to penalize mis-assigned ambiguity
    double penaltyErr() const { return _penerr; }
// error ON RDrift and residual coming from hit t0 error
    double t0Err() const { return _hitt0._t0err*_t2d._vdrift; }
// total error
    double totalErr() const { return _toterr; }
// intrinsic hit error (mm)
    double hitErr() const { return _t2d._rdrifterr; }
// changing the extneral hit error invalidates the cache, it should invalidate the fit, FIXME!!!!
    void setExtErr(double exterr) { _exterr = exterr; }
// set the penalty error: this is set when we can't resolve the ambiguity
    void setPenalty(double penerr) { _penerr = penerr; }
    bool physicalDrift(double maxchi) const;
// logical operators to allow searching for StrawHits
    bool operator == (StrawHit const& sh) const { return _strawhit == sh; }
    bool operator != (StrawHit const& sh) const { return !operator==(sh); }
    void print(std::ostream& ) const;
  protected:
    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj);
    virtual void updateDrift();
    virtual void updateSignalTime();
  //private:
    const StrawHit& _strawhit;
    const Straw& _straw;
    unsigned _istraw;
    TrkLineTraj* _hittraj;
    CLHEP::Hep3Vector _wpos;
    CLHEP::Hep3Vector _wpos_err;
    HitT0 _hitt0;
    double _stime;
    double _exterr,_penerr,_toterr;
    int _iamb;
    bool _ambigupdate;
    T2D _t2d; // current values of t2d
    double _tddist;
    double _tddist_err;
    double _maxdriftpull;
  };

// binary functor to sort TrkStrawHits by StrawHit index
  struct indexcomp : public std::binary_function<TrkStrawHit,TrkStrawHit, bool> {
    bool operator()(const TrkStrawHit* x,const TrkStrawHit* y) { return x->index() < y->index(); }
  };
 
// unary functor to select TrkStrawHit from a given hit
  struct FindTrkStrawHit {
    FindTrkStrawHit(StrawHit const& strawhit) : _strawhit(strawhit) {}
    bool operator () (TrkStrawHit* const& tsh ) { return tsh->strawHit() == _strawhit; }
    StrawHit const& _strawhit;
  };
// define TrkStrawHitVector, to allow explicit conversion and construction
  typedef std::vector<TrkStrawHit*> TrkStrawHitVector;
// utility function to convert vector of TrkHits into TrkStrawHits
  void convert(TrkHitVector const& thv, TrkStrawHitVector& tshv);
}

#endif
