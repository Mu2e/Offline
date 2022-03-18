#ifndef Mu2eG4_validGeometryOrThrow_hh
#define Mu2eG4_validGeometryOrThrow_hh
//
// Called after the G4 geometry has been instantiated.
// Inspect the stores of the geometry objects to check for a valid geometry.
// If the geometry is invalid, throw an exception.
//
// Original author Rob Kutschke
//

namespace mu2e {
  void validGeometryOrThrow ( int printLevel );
}

#endif /* Mu2eG4_validGeometryOrThrow_hh */
