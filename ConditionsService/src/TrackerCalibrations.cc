//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.cc,v 1.3 2011/07/09 05:03:05 mu2ecvs Exp $
// $Author: mu2ecvs $
// $Date: 2011/07/09 05:03:05 $
//

// Mu2e include files
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"

namespace mu2e {

  TrackerCalibrations::TrackerCalibrations( SimpleConfig const& config ){

    // Here we should eventually interface to some database


  }

  const double TrackerCalibrations::TimeDivisionResolution(StrawIndex strawIndex, double znorm) const {
    double resopar0 = 64.2;
    double resopar1 = 60.7;

    double reso  = resopar0 + resopar1 * (znorm - 0.5) * (znorm - 0.5); //resolution in mm
    return reso;

  }

  const double TrackerCalibrations::SignalVelocity(StrawIndex strawIndex) const {
    double distvsdeltat = 231.;
    return distvsdeltat; //mm/ns
  }

  const double TrackerCalibrations::TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const{
    return 0.5 * SignalVelocity(strawIndex) * deltaT;
  }
  
  void TrackerCalibrations::StrawHitInfo(StrawHit const& strawhit,
    CLHEP::Hep3Vector& pos, double& time,double& tdres, double& timeres) const {
    const Tracker& tracker = getTrackerOrThrow();
    const Straw& straw = tracker.getStraw(strawhit.strawIndex());
// compute the position as being on the wire the distance specified by time division.  Note this can
// be beyond the physical wire! 
    double vwire = SignalVelocity(strawhit.strawIndex());
    double tddist = TimeDiffToDistance(strawhit.strawIndex(),strawhit.dt());
    pos = straw.getMidPoint() + tddist*straw.getDirection();    
// this time represents when the particle passed by the wire.
    double shlen = straw.getHalfLength();
    time = strawhit.time() - (shlen-tddist)/vwire;
// the error matrix is defined with x along the wire, y in the mu2e Z direction, and y perpendicular
// error along the wire is given by the time division
    tdres = TimeDivisionResolution(straw.index(),0.5*(shlen-tddist)/shlen);
// time resolution is due to intrinsic timing resolution and time difference resolution
    timeres = tdres/vwire;
  }
}
