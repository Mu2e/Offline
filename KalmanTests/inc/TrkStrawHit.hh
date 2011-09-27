//
// BaBar hit object corresponding to a single straw hit
//
// $Id: TrkStrawHit.hh,v 1.9 2011/09/27 21:49:09 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/09/27 21:49:09 $
//
// Original author David Brown, LBNL
//
#ifndef TrkStrawHit_HH
#define TrkStrawHit_HH
// BaBar
#include "KalmanTests/inc/DetStrawGasElem.hh"
#include "KalmanTests/inc/DetStrawWallElem.hh"
#include "KalmanTests/inc/DetStrawHitType.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/TrkEnums.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkDetElemId.hh"
#include "MatEnv/MatDBInfo.hh"
// Mu2e
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"
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
    double t0, double t0err, double herr);
    virtual ~TrkStrawHit();
//  Simplistic implementation of TrkHitOnTrk interface.  Lie where necessary
    virtual TrkStrawHit* clone(TrkRep* parentRep, const TrkDifTraj* trkTraj = 0) const;
    virtual TrkEnums::TrkViewInfo whatView() const              { return TrkEnums::xyView; }
    virtual unsigned layerNumber() const                        { return _straw.id().getLayer(); }
    virtual const TrkLineTraj* hitTraj() const                   { return _hittraj; }
    virtual bool timeResid(double& t, double& error) const      { return false; }
    virtual bool timeAbsolute(double& t, double& error) const   { return false; }
    virtual int ambig() const { return _iamb; }
    virtual void setAmbig(int newambig);
    unsigned index() const { return _istraw; } // index into StrawHit vector
// Implement the detector-specific casting as a logical-only return value:
// the pointers returned by these functions WILL BE INVALID, but non-null as appropriate, to allow
// hit-counting logic to function
    virtual const DchHitOnTrack* dchHitOnTrack() const { return (DchHitOnTrack*)1;}
    virtual const SvtHitOnTrack* svtHitOnTrack() const { return (SvtHitOnTrack*)0;}
    double hitRMS() const { return _rdrift_err;}
// strawhit specific interface
    const StrawHit& strawHit() const { return _strawhit; }
    const Straw& straw() const { return _straw; }
// correct the hit time for the wire propagation
    double time() const;
    double driftRadius() const { return _rdrift;}
    double driftRadiusErr() const { return _rdrift_err;}
    double timeDiffDist() const { return _tddist; }
    double timeDiffDistErr() const { return _tddist_err; }
    const CLHEP::Hep3Vector& wirePosition() const { return _wpos; }
    const CLHEP::Hep3Vector& wirePositionError() const { return _wpos_err; }
    void hitPosition(CLHEP::Hep3Vector& hpos) const;
    double hitT0() const { return _hitt0;}
    double hitT0Err() const { return _hitt0_err;}
    void updateT0(double hitt0,double hitt0_err);
    double wallPath(Hep3Vector const& tdir) const; // track pathlength through one wall of the straw
    double gasPath(Hep3Vector const& tdir) const; // track pathlength through 1/2 the gas of the straw
// intrinsic hit error
    double hitErr() const { return _herr; }
// allow configuring
    static void setMaxDriftPull(double maxdriftpull) { _maxdriftpull = maxdriftpull; }
// access to associated detector elements
    DetStrawGasElem const& gasElem() const { return _gelem; }
    DetStrawWallElem const& wallElem() const { return _welem; }
// logical operators to allow searching for StrawHits
    bool operator == (StrawHit const& sh) const { return _strawhit == sh; }
    bool operator != (StrawHit const& sh) const { return !operator==(sh); }
  protected:
    TrkStrawHit(const TrkStrawHit& other, TrkRep* rep);
    virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj, bool maintainAmbiguity);
    void updateDrift();
  private:
    const StrawHit& _strawhit;
    const Straw& _straw;
    unsigned _istraw;
    TrkLineTraj* _hittraj;
    CLHEP::Hep3Vector _wpos;
    CLHEP::Hep3Vector _wpos_err;
    double _hitt0, _hitt0_err;
    double _herr;
    int _iamb;
    double _rdrift;
    double _rdrift_err;
    double _tddist;
    double _tddist_err;
    double _wtime;
    double _wtime_err;
// DetModel stuff
    static DetStrawHitType _wtype;
    static DetStrawHitType _gtype;
    DetStrawWallElem _welem;
    DetStrawGasElem _gelem;
// parameters that should come from some service: FIXME!!!
    static double _vdrift;
    static double _maxdriftpull;
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
