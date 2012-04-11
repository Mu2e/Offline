#ifndef ConditionsService_TrackerCalibrations_hh
#define ConditionsService_TrackerCalibrations_hh
//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.hh,v 1.9 2012/04/11 19:48:08 brownd Exp $
// $Author: brownd $
// $Date: 2012/04/11 19:48:08 $
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


  struct TrackerCalibrations: virtual public ConditionsEntity{


    TrackerCalibrations ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

    double TimeDivisionResolution(StrawIndex strawIndex, double znorm) const;
    double SignalVelocity(StrawIndex strawIndex) const;
    // Drift time calibration for reconstruction.  This returns the transverse radius and the error on that, NOT INCLUDING the
    // effective error from uncertainty in the ambiguity assignment.
    // The track direction is used to compute Lorentz effects
    void TimeToDistance(StrawIndex strawIndex, double tdrift, CLHEP::Hep3Vector const& tdir,
      double& rdrift, double& rdrifterr) const;
    // Inverse drift time calibration for use in simulation.  takes the transverse radius and track direction as input.
    void DistanceToTime(StrawIndex strawIndex,double rdrift, CLHEP::Hep3Vector const& tdir,
      double& tdrift, double& tdrifterr) const;

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
