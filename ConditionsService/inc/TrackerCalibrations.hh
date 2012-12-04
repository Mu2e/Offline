#ifndef ConditionsService_TrackerCalibrations_hh
#define ConditionsService_TrackerCalibrations_hh
//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.hh,v 1.15 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
//
// Original author Vadim Rusu
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "DataProducts/inc/StrawIndex.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e
{
  class SimpleConfig;
  class StrawHit;

// simple struct to hold output of timeToDistance function
  struct T2D {
    double _rdrift;
    double _rdrifterr;
    double _vdrift; // local drift velocity at this radius, ie dr/dt(rdrift).
    T2D() : _rdrift(0.0), _rdrifterr(1.0), _vdrift(0.0) {}
  };

// simple struct to hold output of distanceToTime function
  struct D2T {
    double _tdrift;
    double _tdrifterr;
    double _vdrift; // local drift velocity at this radius, ie dr/dt(rdrift).
    D2T() : _tdrift(0.0),_tdrifterr(1.0), _vdrift(0.0) {}
  };

// simple struct to hold output of energyToAmplitude function
  struct E2A {
    double _ampl;
    double _amplerr;
    E2A() : _ampl(0.0), _amplerr(1.0) {}
  };

// simple struct to hold output of amplitudeToEnergy function
  struct A2E {
    double _edep;
    double _edeperr;
    A2E() : _edep(0.0),_edeperr(1.0) {}
  };

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
    virtual void StrawHitInfo(StrawHit const& strawhit,
      CLHEP::Hep3Vector& pos, double& time,double& tdres, double& timeres) const;

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
    double _rres;

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
