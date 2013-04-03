//
// Parameters for I-tracker calibrations.
//
// $Id: TrackerCalibrationsI.cc,v 1.3 2013/04/03 22:08:21 tassiell Exp $
// $Author: tassiell $
// $Date: 2013/04/03 22:08:21 $
//

// Mu2e include files
#include "ConditionsService/inc/TrackerCalibrationsI.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include <iostream>

#define TwoInvSqrt12 0.5773502691896258

namespace mu2e {

  TrackerCalibrationsI::TrackerCalibrationsI( SimpleConfig const& config ) : TrackerCalibrations(config)  {}
  
  double TrackerCalibrationsI::TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const{
    if (deltaT>-9999) {
            return 0.5 * SignalVelocity(strawIndex) * deltaT;
    } else {
            return 0.0;
    }
  }

  void TrackerCalibrationsI::StrawHitInfo(Straw const& straw, StrawHit const& strawhit, SHInfo& shinfo) const {
          const Tracker& tracker = getTrackerOrThrow();
          const ITracker& itr     = static_cast<const ITracker&>(tracker);
          CellGeometryHandle *itwp = itr.getCellGeometryHandle();
          itwp->SelectCellDet(strawhit.strawIndex().asUint());
          // compute the position as being on the wire the distance specified by time division.  Note this can
          // be beyond the physical wire!
          double vwire = SignalVelocity(strawhit.strawIndex());
          shinfo._tddist = TimeDiffToDistance(strawhit.strawIndex(),strawhit.dt());
          shinfo._pos = itwp->GetWireCenter() + shinfo._tddist*itwp->GetWireDirection();
          // this time represents when the particle passed by the wire.
          double shlen = itwp->GetCellHalfLength();
          shinfo._time = strawhit.time() - (shlen-shinfo._tddist)/vwire;
          // the error matrix is defined with x along the wire, y in the mu2e Z direction, and y perpendicular
          // error along the wire is given by the time division
          if (strawhit.dt()>-9999) {
                  shinfo._tdres = TimeDivisionResolution(strawhit.strawIndex(),0.5*(shlen-shinfo._tddist)/shlen);
          } else {
                  shinfo._tdres = shlen*TwoInvSqrt12;
          }
          // time resolution is due to intrinsic timing resolution and time difference resolution
          shinfo._timeres = shinfo._tdres/vwire;
  }

}
