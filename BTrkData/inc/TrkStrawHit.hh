
// BaBar hit object corresponding to a single straw hit
//
// Original author David Brown, LBNL
//
#ifndef TrkStrawHit_HH
#define TrkStrawHit_HH
// BTrk
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "Offline/BTrkLegacy/inc/TrkT0.hh"
// Mu2e
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// c++
#include <vector>
#include <functional>
// forward refs
class TrkDifTraj;
class TrkDifPieceTraj;

namespace mu2e
{
  class TrkStrawHit : public TrkHit {
  public:
  // enum for hit flags
    TrkStrawHit(StrawResponse::cptr_t strawResponse,
                const ComboHit& strawhit, const Tracker& tracker,StrawHitIndex index,
                const TrkT0& trkt0, double fltlen, double maxdriftpull,
                double timeWeight);
    virtual ~TrkStrawHit();
//  implementation of TrkHit interface
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }
    int ambig() const { return _iamb; }
    virtual void setAmbig(int newambig);
    StrawHitIndex index() const { return _index; } // index into StrawHit vector
    double hitRMS() const { return _rdrifterr;}
// strawhit specific interface
    const ComboHit& comboHit() const { return _combohit; }
    const Straw& straw() const { return _straw; }
    virtual double time() const { return _combohit.time(); }
    //    virtual bool time(HitT0& t0) ;// 2019-04-22: GIANIPEZ change
    StrawEnd const& earlyEnd() const { return _combohit.earlyEnd(); }
    double driftTime() const; // drift time for the current end strategy
    double driftPhi() const { return _phi;}
    double driftRadius() const { return _rdrift;}
    double driftRadiusErr() const { return _rdrifterr;}
    double driftVelocity() const { return _vdriftinst; }
    double timeDiffDist() const { return _combohit.wireDist(); }
    double timeDiffDistErr() const { return _combohit.wireRes(); }
    const CLHEP::Hep3Vector& wirePosition() const { return _wpos; }
    void hitPosition(CLHEP::Hep3Vector& hpos) const;
    virtual bool signalPropagationTime( TrkT0& t0 ); // time and error
    virtual void trackT0Time(double &htime, double t0flt, const TrkDifPieceTraj* ptraj, double vflt);

    double signalTime() const { return _stime; } // time for signal to reach the end of the wire
// error to penalize mis-assigned ambiguity
    double penaltyErr() const { return _penerr; }
// error ON RDrift and residual coming from hit t0 error
    double t0Err() const { return hitT0()._t0err*_vdriftinst; }
// total error
    double totalErr() const { return _toterr; }
// intrinsic hit error (mm)
    virtual double hitErr() const { return _rdrifterr; }
  // test the consistincy of this hit with 'physical' limts, with a given # of sigma
    virtual bool isPhysical(double maxchi) const;
    virtual void print(std::ostream& ) const;

    //**************************************************
    // SET VALUES
    //**************************************************
    // set the penalty error: this is set when we can't resolve the ambiguity
    void setPenalty(double penerr) { _penerr = penerr; }

  protected:

    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj);
    virtual void updateDrift();
    virtual void updateSignalTime();
  //private:

    StrawResponse::cptr_t _strawResponse;

    const ComboHit&   _combohit;
    const Straw&      _straw;
    double              _rstraw; // straw radius; cached from Tracker
    StrawHitIndex     _index;
    TrkLineTraj*      _hittraj;
    CLHEP::Hep3Vector _wpos;
    CLHEP::Hep3Vector _wpos_err;
    double            _penerr,_toterr;
    int               _iamb;
    double            _rdrift;
    double            _rdrifterr;
    double            _phi;
    double            _vdriftinst;
    double            _vprop; // effective signal propagation velocity
    double              _stime; // signal propagation time
    double            _maxdriftpull;

  };

// binary functor to sort TrkStrawHits by StrawHit index
  struct indexcomp {
    bool operator()(const TrkStrawHit* x,const TrkStrawHit* y) { return x->index() < y->index(); }
  };

// unary functor to select TrkStrawHit from a given hit
  struct FindTrkStrawHit {
    FindTrkStrawHit(ComboHit const& strawhit) : _combohit(strawhit) {}
    bool operator () (TrkStrawHit* const& tsh ) { return &tsh->comboHit() == &_combohit; }
    ComboHit const& _combohit;
  };
// define TrkStrawHitVector, to allow explicit conversion and construction
  typedef std::vector<TrkStrawHit*> TrkStrawHitVector;
// utility function to convert vector of TrkHits into TrkStrawHits
  void convert(TrkHitVector const& thv, TrkStrawHitVector& tshv);
}

#endif
