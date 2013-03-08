#ifndef RecoDataProducts_StrawHitPosition_hh
#define RecoDataProducts_StrawHitPosition_hh
//
// Class to describe derived information from a StrawHit, in particular pos().
// 
// $Id: StrawHitPosition.hh,v 1.1 2013/03/08 04:29:49 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 04:29:49 $
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
namespace mu2e {

  class StrawHitPosition {
    public:
// construct from a straw hit.  Calibration and other information must be provided from outside
      StrawHitPosition(StrawHit const& hit, Tracker const& tracker, ConditionsHandle<TrackerCalibrations>& tcal, StrawHitFlag const& flag=_nullflag );
// construct from another hit, optionally with additional flag bits (these will be ORed with the existing bits
      StrawHitPosition(StrawHitPosition const& pos,StrawHitFlag const& orflag=_nullflag);
// construct from a stereo hit
      StrawHitPosition(StereoHit const& sthit);
// null constructor for root
      StrawHitPosition();
      virtual ~StrawHitPosition();
      StrawHitPosition& operator =(StrawHitPosition const& other);
      enum edir {phi=0,rho=1};
// accessors
      CLHEP::Hep3Vector const& pos() const { return _pos; }
      float posRes(edir dir) const { return dir==phi ? _pres : _rres; }
      float chisq() const { return _chisq; }
      StrawHitFlag const& flag() const { return _flag; }
    private:
      CLHEP::Hep3Vector _pos; // pos() of the hit
      float _pres; // pos resolution along phi
      float _rres; // pos resolution perpendicular to the Z axis
      float _chisq; // self-consistency chisquared
      StrawHitFlag _flag; // bit flags for this hit
      static StrawHitFlag _nullflag; // null flag
      static double _invsqrt12; // 1/sqrt(12)
  };
}
#endif


