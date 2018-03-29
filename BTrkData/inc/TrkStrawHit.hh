
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
#include "BTrk/TrkBase/TrkT0.hh"
// Mu2e
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
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
    enum enduse { cal=TrkTypes::cal, hv = TrkTypes::hv, earliest, both};
    TrkStrawHit(const StrawHit& strawhit, const Straw& straw,StrawHitIndex index,
		const TrkT0& trkt0, double fltlen, double exterr, double maxdriftpull, 
		double timeWeight, double mint0doca);
    virtual ~TrkStrawHit();
//  implementation of TrkHit interface
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }
    int ambig() const { return _iamb; }
    enduse driftEnd() const { return _enduse; }
//    virtual void invert();
    virtual void setAmbig(int newambig);
    void setAmbigUpdate(bool update) { _ambigupdate = update; }
    StrawHitIndex index() const { return _index; } // index into StrawHit vector
    double hitRMS() const { return _rdrifterr;}
// strawhit specific interface
    const StrawHit& strawHit() const { return _strawhit; }
    const Straw& straw() const { return _straw; }
    virtual double time() const { return _strawhit.time(); }
// the following function is DEPRECATED as the underlying function is now end specific
    double driftTime(StrawEnd end) const; // drift time for a specific end
    double driftTime() const; // drift time for the current end strategy

    double driftRadius() const { return _rdrift;}
    double driftRadiusErr() const { return _rdrifterr;}
    double driftVelocity() const { return _vdriftinst; }
    double timeDiffDist() const { return _tddist; }
    double timeDiffDistErr() const { return _tddist_err; }
    const CLHEP::Hep3Vector& wirePosition() const { return _wpos; }
    void hitPosition(CLHEP::Hep3Vector& hpos) const;
    virtual bool signalPropagationTime(double &propTime, double&Doca, 
			       double resid, double &residErr, 
			       CLHEP::Hep3Vector trajDirection);//propagation time
    virtual void trackT0Time(double &htime, double t0flt, const TrkDifPieceTraj* ptraj, double vflt);

    double signalTime(StrawEnd end=TrkTypes::cal) const { return _stime[end]; } // time for signal to reach the end of the wire
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
    
// logical operators to allow searching for StrawHits
    bool operator == (StrawHit const& sh) const { return _strawhit == sh; }
    bool operator != (StrawHit const& sh) const { return !operator==(sh); }
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
    const StrawHit&   _strawhit;
    const Straw&      _straw;
    StrawHitIndex     _index;
    TrkLineTraj*      _hittraj;
    CLHEP::Hep3Vector _wpos;
    CLHEP::Hep3Vector _wpos_err;
    double _stime[2]; // time for the signal to get from the POCA to each wire end
    double            _penerr,_toterr;
    int               _iamb;
    enduse _enduse; // which ends are used in the drift measurement
    bool              _ambigupdate;
    double            _rdrift;
    double            _rdrifterr;
    double            _tddist;
    double            _tddist_err;
    double            _vdriftinst;
    double            _maxdriftpull;
    double            _mint0doca;	    // minimum doca for t0 calculation.  Note this is a SIGNED QUANTITITY
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
