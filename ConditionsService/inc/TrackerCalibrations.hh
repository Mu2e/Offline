#ifndef ConditionsService_TrackerCalibrations_hh
#define ConditionsService_TrackerCalibrations_hh
//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.hh,v 1.18 2014/03/01 11:13:30 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/01 11:13:30 $
//
// Original author Vadim Rusu
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "ConditionsService/inc/TrackerCalibrationStructs.hh"
#include "DataProducts/inc/StrawIndex.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e
{
  class SimpleConfig;
  class StrawHit;
  class Straw;
  struct TrackerCalibrations: virtual public ConditionsEntity{


    TrackerCalibrations ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

    double TimeDivisionResolution(StrawIndex strawIndex, double znorm) const;
    double SignalVelocity(StrawIndex strawIndex) const;
    // Drift time calibration for reconstruction.  This returns the transverse radius and the error on that, NOT INCLUDING the
    // effective error from uncertainty in the ambiguity assignment or t0 determination.
    // The track direction is used to compute Lorentz effects
    void TimeToDistance(StrawIndex strawIndex, double tdrift, CLHEP::Hep3Vector const& tdir, T2D& t2d) const;
    // Inverse drift time calibration for use in simulation.  takes the transverse radius and track direction as input.
    void DistanceToTime(StrawIndex strawIndex,double rdrift, CLHEP::Hep3Vector const& tdir, D2T& d2t) const;

    // time difference calibration
    virtual double TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const;
    // information about a hit's position and time.  This uses time difference to compute
    // the position along the wire
    virtual void StrawHitInfo(Straw const& straw, StrawHit const& strawhit, SHInfo& shinfo) const;
    void EnergyToAmplitude(StrawIndex strawIndex, double edep, E2A& e2a) const;
    void AmplitudeToEnergy(StrawIndex strawIndex, double ampl, A2E& a2e) const;

    double CrossTalk(StrawIndex strawIndex0, StrawIndex strawIndexN) const;

  protected:

    // We want to discourage multi-phase construction.
    TrackerCalibrations ();

    // time-division base resolution and length-dependent quadratic term
    double _tdresopar0;
    double _tdresopar1;
    // temoprary constant drift velocity and time resolution.  Replace these with a more physical model	FIXME!!!
    double _vdrift;
    double _rres_min, _rres_max, _rres_rad;

    //velocity of signal propagation in wire mm/ns
    double _distvsdeltat;

    double _edepToAmpl; // MeV/mV
    double _amplRes;    // relative
    double _crossTalk;  // amount of crosstalk

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const TrackerCalibrations&  ){
    ost << "( TrackerCalibrations: to be implemnted "
        << " )";

    return ost;
  }
}

#endif /* ConditionsService_TrackerCalibrations_hh */
