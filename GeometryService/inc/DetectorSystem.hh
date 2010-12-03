#ifndef DetectorSystem_HH
#define DetectorSystem_HH
//
// Transformations between the Mu2e coordinate system and the detector coordinate system.
// This is a dumb data class that does not know how to build itself.
//
// $Id: DetectorSystem.hh,v 1.1 2010/12/03 00:52:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/03 00:52:40 $
//
// Original author Rob Kutschke
//

// Mu2e includes.
#include "GeometryService/inc/Detector.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class DetectorSystem: public Detector{

  public:

    DetectorSystem( CLHEP::Hep3Vector const& origin ):
      _origin(origin){
    }

    // Accept the compiler generator d'tor, copy c'tor and copy assigment.

    // Take a 3-vector in the detector system; return it in the Mu2e system.
    CLHEP::Hep3Vector toMu2e( CLHEP::Hep3Vector const& v ){
      return v+_origin;
    }

    // Take a 3-vector in the mu2e system; return it in the detector system.
    CLHEP::Hep3Vector toDetector( CLHEP::Hep3Vector const& v ){
      return v-_origin;
    }

  private:

    // Location of the origin of the Detector coordinate system, as measured
    // in the Mu2e coordinate system.
    CLHEP::Hep3Vector _origin;

  };

} //namespace mu2e

#endif
