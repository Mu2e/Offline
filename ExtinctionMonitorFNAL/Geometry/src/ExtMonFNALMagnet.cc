#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

//================================================================
mu2e::ExtMonFNALMagnet::ExtMonFNALMagnet()
  : outerHalfSize_()
  , apertureWidth_()
  , apertureHeight_()
  , magneticLength_()
  , nominalMomentum_()
{}

//================================================================
double mu2e::ExtMonFNALMagnet::trackBendRadius(double momentum) const {
  // In the bend plane: compute the gyroradius
  // The constant factor is 1/c_light scaled such as
  // to get rTrack in millimeters
  const double rTrack = 3335.64095198 * (momentum/CLHEP::GeV) / (bfield_.mag()/CLHEP::tesla);
  return rTrack;
}

//================================================================
double mu2e::ExtMonFNALMagnet::trackBendHalfAngle(double momentum) const {

  const double rTrack = trackBendRadius(momentum);

  // Can't do momenta that are too low.  For simplicity we just
  // check for the "absolutely impossible" requests here.  The real
  // momentum constraint is tighter because of other pieces of
  // geometry.

  if(rTrack < 0.5*magneticLength_) {
    throw cet::exception("GEOM")<<"ExtMonFNALBuilding::FilterMagnetExtMonFNAL::trackBendHalfAngle(): "
                                <<"requested momentum p="<<momentum/CLHEP::GeV<<" GeV is too low ";
  }

  return asin(0.5*magneticLength_/rTrack);
}

//================================================================
