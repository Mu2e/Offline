#ifndef GeneralUtilities_OrientationResolver_hh
#define GeneralUtilities_OrientationResolver_hh

// ********************************************************
// This is a utility class used for parsing geometry
// orientation strings for Mu2e and generating the needed
// HEPRotations.

// Original Author:  David Norvil Brown
// University of Louisville, October 2014.  Moved to
// separate file in May 2016.

// ********************************************************

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <string>

namespace mu2e {

  class OrientationResolver
  {
  public:
    
    OrientationResolver() {}
    ~OrientationResolver() {}

    void getRotationFromOrientation ( CLHEP::HepRotation& aRotation,
				    std::string orient );

  private:
    
  };
}  // end of namespace

#endif // ExternalShieldingGeom/inc/ExtShieldDownstream.hh
