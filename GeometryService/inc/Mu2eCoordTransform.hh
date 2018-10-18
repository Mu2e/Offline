#ifndef GeometryService_Mu2eCoordTransform_HH
#define GeometryService_Mu2eCoordTransform_HH
// **********************************************************************
// Mu2eCoordTransform.hh
// This is the header file for the simple Mu2eCoordTransform class.

// Mu2eCoordTransform class
// This could be a simple function, except that some configurable
// values have to be read in at run-time and stored for later use.
// This class allows one to transform positions from one coordinate
// system to another.  Specifically, it works with the Mu2e coordinate
// system, the g4beamline (g4bl) system, and the Tracker system.
// As a reminder, the systems are defined as follows:

// Mu2e coordinate system:  the origin is at the center of TS3 with
// positive x in the building north direction, positive y vertically upward,
// and positive z in the building east direction.

// G4beamline coordinate system:  the axis directions are the same, but
// the origin is at the center of the production target.

// Tracker coordinate system:  the axis directions are the same, but the
// origin is at a convenient point in the Detector Solenoid.  Currently,
// that convenient point is the center of the Tracker active volume.

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class SimpleConfig;

  class Mu2eCoordTransform {
  public:

    Mu2eCoordTransform( const SimpleConfig& config );
    ~Mu2eCoordTransform( ) {};

    CLHEP::Hep3Vector Mu2eToDetec  ( const CLHEP::Hep3Vector& ptInMu2e );
    CLHEP::Hep3Vector DetecToMu2e  ( const CLHEP::Hep3Vector& ptInDetec  );

    CLHEP::Hep3Vector Mu2eToG4bl ( const CLHEP::Hep3Vector& ptInMu2e );
    CLHEP::Hep3Vector G4blToMu2e ( const CLHEP::Hep3Vector& ptInG4bl );

    CLHEP::Hep3Vector DetecToG4bl  ( const CLHEP::Hep3Vector& ptInDetec  );
    CLHEP::Hep3Vector G4blToDetec  ( const CLHEP::Hep3Vector& ptInG4bl );
    
  private:

    double DetecZ0InMu2e;
    double DetecX0InMu2e;

    double G4blZ0InMu2e;
    double G4blX0InMu2e;

    
  };  // end of class def

} // end of namespace


#endif
