#ifndef ConditionsService_TrackerCalibrationsI_hh
#define ConditionsService_TrackerCalibrationsI_hh
//
// Parameters for I-tracker calibrations.
//
// $Id: TrackerCalibrationsI.hh,v 1.3 2013/04/03 22:08:21 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/04/03 22:08:21 $
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
    virtual void StrawHitInfo(Straw const& straw, StrawHit const& strawhit, SHInfo& shinfo) const;

  private:

    // We want to discourage multi-phase construction.
    TrackerCalibrationsI ();
  };
}

#endif /* ConditionsService_TrackerCalibrationsI_hh */
