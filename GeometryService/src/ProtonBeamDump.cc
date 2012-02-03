#include "GeometryService/inc/ProtonBeamDump.hh"

#include <cmath>

#include "cetlib/exception.h"

namespace mu2e {

  double ProtonBeamDump::FilterMagnetExtMonFNAL::trackBendHalfAngle(double momentum) const {

    // In the bend plane: compute the gyroradius
    // The constant factor is 1/c_light scaled such as
    // to get rTrack in millimeters
    const double rTrack = 3335.64095198 * (momentum/CLHEP::GeV) / (_fieldStrength/CLHEP::tesla);
    
//    std::cerr<<"AG: got rTrack = "<<rTrack<<" mm for p = "
//	     <<(momentum/CLHEP::GeV)<<" GeV and  B = "
//	     <<(_fieldStrength()/CLHEP::tesla)<<" tesla"<<std::endl;

    // Can't do momenta that are too low.  For simplicity we just
    // check for the "absolutely impossible" requests here.  The real
    // momentum constraint is tighter because of other pieces of
    // geometry.

    if(_outerHalfSize[2] < rTrack) {
      return asin(_outerHalfSize[2]/rTrack);
    }
    else {
      throw cet::exception("GEOM")<<"ProtonBeamDump::FilterMagnetExtMonFNAL::trackBendHalfAngle(): "
				  <<"requested momentum p="<<momentum/CLHEP::GeV<<" GeV is too low ";
    }
  }


  CLHEP::Hep3Vector ProtonBeamDump::filterEntranceInMu2e() const {
    CLHEP::Hep3Vector filterEntranceInEnclosure(_coreCenterInEnclosure[0] + _filterEntranceOffsetX,
                                                _coreCenterInEnclosure[1] + _filterEntranceOffsetY,
                                                _enclosureHalfSize[2]);
    
    return _enclosureCenterInMu2e + _enclosureRotationInMu2e.inverse() * filterEntranceInEnclosure;
  }
}
