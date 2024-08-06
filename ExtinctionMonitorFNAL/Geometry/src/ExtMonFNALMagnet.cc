#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib_except/exception.h"

//================================================================
mu2e::ExtMonFNALMagnet::ExtMonFNALMagnet()
  : outerHalfSize_()
  , apertureWidth_()
  , apertureHeight_()
  , magneticLength_()
  , nominalMomentum_()
{}

//================================================================
namespace {
  // 1/c_light in units that give rTrack in millimeters
  // for input momentum in MeV/c in trackBendRadius()
  static const double c_light_inv = 3335.64095198;
}

double mu2e::ExtMonFNALMagnet::trackBendRadius(double momentum) const {
  // In the bend plane: compute the gyroradius
  const double rTrack = c_light_inv * (momentum/CLHEP::GeV) / (bfield_.mag()/CLHEP::tesla);
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
CLHEP::Hep2Vector mu2e::ExtMonFNALMagnet::dxdzdydz() const {
  const auto v = magnetRotationInMu2e_ * CLHEP::Hep3Vector(0,0, -1);
  return CLHEP::Hep2Vector( v.x()/v.z(), v.y()/v.z());
}

//================================================================
