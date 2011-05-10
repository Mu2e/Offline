#ifndef TrackerCalibrations_H
#define TrackerCalibrations_H
//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.hh,v 1.1 2011/05/10 16:44:02 vrusu Exp $
// $Author: vrusu $
// $Date: 2011/05/10 16:44:02 $
//
// Original author Vadim Rusu
//

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "ConditionsService/inc/ConditionsEntity.hh"
#include "TrackerGeom/inc/StrawIndex.hh"


namespace mu2e
{
  class SimpleConfig;


  struct TrackerCalibrations: public ConditionsEntity{


    TrackerCalibrations ( SimpleConfig const& config );

    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

    const double TimeDivisionResolution(StrawIndex strawIndex, double znorm) const;
    const double SignalVelocity(StrawIndex strawIndex) const;

    //this shoule be called by the patt rec
    const double TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const;

  private:

    // We want to discourage multi-phase construction.
    TrackerCalibrations ();

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const TrackerCalibrations& daqpar ){
    ost << "( "
        << " )";

    return ost;
  }
}

#endif
