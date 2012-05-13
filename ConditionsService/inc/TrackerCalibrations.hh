#ifndef ConditionsService_TrackerCalibrations_hh
#define ConditionsService_TrackerCalibrations_hh
//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.hh,v 1.11 2012/05/13 21:26:32 ignatov Exp $
// $Author: ignatov $
// $Date: 2012/05/13 21:26:32 $
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
  };

// simple struct to hold output of distanceToTime function
  struct D2T {
    double _tdrift;
    double _tdrifterr;
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
    double TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const;
    // information about a hit's position and time.  This uses time difference to compute
    // the position along the wire
    void StrawHitInfo(StrawHit const& strawhit,
      CLHEP::Hep3Vector& pos, double& time,double& tdres, double& timeres) const;

  private:

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
