// Mu2eCoordTransform.cc
// This is the definitions file for a simple class to perform coordinate
// transforms among the coordinate systems used in Mu2e.
// It is owned by the GeometryService.

// David Norvil Brown (UofL) August 2017

#include "GeometryService/inc/Mu2eCoordTransform.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  Mu2eCoordTransform::Mu2eCoordTransform( const SimpleConfig& config )
    : DetecZ0InMu2e ( config.getDouble("mu2e.detectorSystemZ0") ),
      G4blX0InMu2e  ( config.getDouble("mu2e.solenoidOffset") )
  {
    DetecX0InMu2e = -G4blX0InMu2e;
    G4blZ0InMu2e  = -7929.0; // mm, for now assumed this is fixed

  } // end of constructor

  CLHEP::Hep3Vector Mu2eCoordTransform::Mu2eToDetec( const CLHEP::Hep3Vector& ptInMu2e )
  {
    return CLHEP::Hep3Vector( ptInMu2e.x() - DetecX0InMu2e, ptInMu2e.y(),
			      ptInMu2e.z() - DetecZ0InMu2e );

  }// end of Mu2eToDetec function

  CLHEP::Hep3Vector Mu2eCoordTransform::DetecToMu2e( const CLHEP::Hep3Vector& ptInDetec )
  {
    return CLHEP::Hep3Vector( ptInDetec.x() + DetecX0InMu2e, ptInDetec.y(),
			      ptInDetec.z() + DetecZ0InMu2e );

  } // end of DetecToMu2e function


  CLHEP::Hep3Vector Mu2eCoordTransform::Mu2eToG4bl( const CLHEP::Hep3Vector& ptInMu2e )
  {
    return CLHEP::Hep3Vector( ptInMu2e.x() - G4blX0InMu2e, ptInMu2e.y(),
			      ptInMu2e.z() - G4blZ0InMu2e );

  } // end of Mu2eToG4bl

  CLHEP::Hep3Vector Mu2eCoordTransform::G4blToMu2e( const CLHEP::Hep3Vector& ptInG4bl )
  {
    return CLHEP::Hep3Vector( ptInG4bl.x() + G4blX0InMu2e, ptInG4bl.y(),
			      ptInG4bl.z() + G4blZ0InMu2e );

  } // end of G4blToMu2e

  CLHEP::Hep3Vector Mu2eCoordTransform::DetecToG4bl( const CLHEP::Hep3Vector& ptInDetec )
  {
    return CLHEP::Hep3Vector( ptInDetec.x() + DetecX0InMu2e - G4blX0InMu2e,
			      ptInDetec.y(),
			      ptInDetec.z() + DetecZ0InMu2e - G4blZ0InMu2e );
  } // end of DetecToG4bl

  CLHEP::Hep3Vector Mu2eCoordTransform::G4blToDetec( const CLHEP::Hep3Vector& ptInG4bl )
  {
    return CLHEP::Hep3Vector( ptInG4bl.x() + G4blX0InMu2e - DetecX0InMu2e,
			      ptInG4bl.y(),
			      ptInG4bl.z() + G4blZ0InMu2e - DetecZ0InMu2e );

  } // end of G4blToDetec


}
