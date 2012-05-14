//
// BaBar hit object corresponding to a single straw hit
//
// $Id: TrkStrawHit.hh,v 1.14 2012/05/14 19:20:02 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/05/14 19:20:02 $
//
// Original author David Brown, LBNL
//
#ifndef TrkStrawHit_HH
#define TrkStrawHit_HH
// BaBar
#include "KalmanTests/inc/DetStrawGasElem.hh"
#include "KalmanTests/inc/DetStrawWallElem.hh"
#include "KalmanTests/inc/DetStrawHitType.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/TrkEnums.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkDetElemId.hh"
#include "MatEnv/MatDBInfo.hh"
// Mu2e
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// forward refs
class TrkRep;
class TrkDifTraj;
class DchHitOnTrack;
class SvtHitOnTrack;

namespace mu2e 
{
// dummy hit to satisfy BaBar interface; FIXME!!
  class TrkDummyHit : public TrkFundHit{
  public:
    TrkDummyHit(TrkEnums::TrkViewInfo view, int id, TrkDetElemId::systemIndex sys);
    virtual ~TrkDummyHit();
    TrkDummyHit* clone() const;
    // obsolete inteface
    virtual TrkEnums::TrkViewInfo whatView() const      { return _view; }
    virtual TrkDetElemId elemId() const                 { return _eid; }
    virtual const GTrack* getGTrack() const             { return NULL; }
    // fake the layer number
    int layerNumber() const { return _eid.systemElemId();}
    virtual void print(std::ostream&) const {}
  protected:
    TrkDummyHit(const TrkDummyHit& other);
    TrkEnums::TrkViewInfo _view;
    TrkDetElemId _eid;
  };

  class TrkStrawHit : public TrkHitOnTrk {
  public:
    TrkStrawHit(const StrawHit& strawhit, const Straw& straw,unsigned istraw,
    const TrkT0& trkt0, double flt0, double fltlen, double exterr, double maxdriftpull);
    virtual ~TrkStrawHit();
//  Simplistic implementation of TrkHitOnTrk interface.  Lie where necessary
    virtual TrkStrawHit* clone(TrkRep* parentRep, const TrkDifTraj* trkTraj = 0) const;
    virtual TrkEnums::TrkViewInfo whatView() const              { return TrkEnums::xyView; }
    virtual unsigned layerNumber() const                        { return _straw.id().getLayer(); }
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }
    virtual bool timeResid(double& t, double& error) const      { return false; } // implement these:: FIXME!!!
    virtual bool timeAbsolute(double& t, double& error) const   { return false; }
    virtual int ambig() const { return _iamb; }
    virtual void setAmbig(int newambig);
    void setAmbigUpdate(bool update) { _ambigupdate = update; }
    unsigned index() const { return _istraw; } // index into StrawHit vector
// Implement the detector-specific casting as a logical-only return value:
// the pointers returned by these functions WILL BE INVALID, but non-null as appropriate, to allow
// hit-counting logic to function
    virtual const DchHitOnTrack* dchHitOnTrack() const { return (DchHitOnTrack*)1;}
    virtual const SvtHitOnTrack* svtHitOnTrack() const { return (SvtHitOnTrack*)0;}
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
    double hitT0() const { return _hitt0;}
    double hitT0Err() const { return _hitt0_err;}
    void updateT0(const TrkT0&, double t0flt);
    double wallPath(double pdist,Hep3Vector const& tdir) const; // track pathlength through one wall of the straw
    double gasPath(double pdist,Hep3Vector const& tdir) const; // track pathlength through 1/2 the gas of the straw
// external hit error (mm); the intrinsic error comes from the t2d calibration object
    double extErr() const { return _exterr; }
// error to penalize mis-assigned ambiguity
    double penaltyErr() const { return _penerr; }
// error ON RDrift and residual coming from hit t0 error
    double t0Err() const { return _hitt0_err*_t2d._vdrift; }
// total error
    double totalErr() const { return _toterr; }
// intrinsic hit error (mm)
    double hitErr() const { return _t2d._rdrifterr; }
// changing the extneral hit error invalidates the cache, it should invalidate the fit, FIXME!!!! 
    void setExtErr(double exterr) { _exterr = exterr; }
// set the penalty error: this is set when we can't resolve the ambiguity
    void setPenalty(double penerr) { _penerr = penerr; }
// access to associated detector elements
    DetStrawGasElem const& gasElem() const { return _gelem; }
    DetStrawWallElem const& wallElem() const { return _welem; }
// logical operators to allow searching for StrawHits
    bool operator == (StrawHit const& sh) const { return _strawhit == sh; }
    bool operator != (StrawHit const& sh) const { return !operator==(sh); }
  protected:
    TrkStrawHit(const TrkStrawHit& other, TrkRep* rep);
    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj);
    void updateDrift();
  private:
    const StrawHit& _strawhit;
    const Straw& _straw;
    unsigned _istraw;
    TrkLineTraj* _hittraj;
    CLHEP::Hep3Vector _wpos;
    CLHEP::Hep3Vector _wpos_err;
    double _hitt0, _hitt0_err;
    double _exterr,_penerr,_t0err,_toterr;
    int _iamb;
    bool _ambigupdate;
    T2D _t2d; // current values of t2d
    double _tddist;
    double _tddist_err;
    double _wtime;
    double _wtime_err;
    double _maxdriftpull;
// DetModel stuff
    static DetStrawHitType _wtype;
    static DetStrawHitType _gtype;
    DetStrawWallElem _welem;
    DetStrawGasElem _gelem;
// parameters that should come from some service: FIXME!!!
    static const double _vlight;
    static MatDBInfo* _matdbinfo;
  };
// unary functor to select TrkStrawHit from a given hit
  struct FindTrkStrawHit {
    FindTrkStrawHit(StrawHit const& strawhit) : _strawhit(strawhit) {}
    bool operator () (TrkStrawHit* const& tsh ) { return tsh->strawHit() == _strawhit; }
    StrawHit const& _strawhit;
  };

  
}

#endif
