//
// Parameters for I-tracker calibrations.
//
// $Id: TrackerCalibrationsI.cc,v 1.2 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
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

  void TrackerCalibrationsI::StrawHitInfo(StrawHit const& strawhit,
    CLHEP::Hep3Vector& pos, double& time,double& tdres, double& timeres) const {
    const Tracker& tracker = getTrackerOrThrow();
    const ITracker& itr     = static_cast<const ITracker&>(tracker);
    CellGeometryHandle *itwp = itr.getCellGeometryHandle();
    itwp->SelectCellDet(strawhit.strawIndex().asUint());
// compute the position as being on the wire the distance specified by time division.  Note this can
// be beyond the physical wire!
    double vwire = SignalVelocity(strawhit.strawIndex());
    double tddist = TimeDiffToDistance(strawhit.strawIndex(),strawhit.dt());
    pos = itwp->GetWireCenter() + tddist*itwp->GetWireDirection();
// this time represents when the particle passed by the wire.
    double shlen = itwp->GetCellHalfLength();
    time = strawhit.time() - (shlen-tddist)/vwire;
// the error matrix is defined with x along the wire, y in the mu2e Z direction, and y perpendicular
// error along the wire is given by the time division
    if (strawhit.dt()>-9999) {
            tdres = TimeDivisionResolution(strawhit.strawIndex(),0.5*(shlen-tddist)/shlen);
    } else {
            tdres = shlen*TwoInvSqrt12;
    }
    // time resolution is due to intrinsic timing resolution and time difference resolution
    timeres = tdres/vwire;
  }

}
