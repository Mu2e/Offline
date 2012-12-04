#ifndef ConditionsService_TrackerCalibrationsI_hh
#define ConditionsService_TrackerCalibrationsI_hh
//
// Parameters for I-tracker calibrations.
//
// $Id: TrackerCalibrationsI.hh,v 1.2 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "ConditionsService/inc/TrackerCalibrations.hh"

namespace mu2e
{
  class SimpleConfig;
  class StrawHit;
  struct TrackerCalibrationsI: public TrackerCalibrations {


    TrackerCalibrationsI ( SimpleConfig const& config );
    // time difference calibration
    virtual double TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const;
    // information about a hit's position and time.  This uses time difference to compute
    // the position along the wire
    virtual void StrawHitInfo(StrawHit const& strawhit,
      CLHEP::Hep3Vector& pos, double& time,double& tdres, double& timeres) const;

  private:

    // We want to discourage multi-phase construction.
    TrackerCalibrationsI ();
  };
}

#endif /* ConditionsService_TrackerCalibrationsI_hh */
