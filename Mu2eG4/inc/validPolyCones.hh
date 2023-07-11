#ifndef Mu2eG4_validPolyCones_hh
#define Mu2eG4_validPolyCones_hh
//
// Called after the G4 geometry has been instantiated.
// Check all PolyCones in the G4SolidStore.
//   - Return true if the Polycone geometry is valid.
//   - Return false if any PolyCone has an invalid geometry
//
// Conditions checked for:
//  1) The geometry is invalid if any two planes have identical ( z, rmin, rmax).
//
// Original author Rob Kutschke
//

namespace mu2e {
  bool validPolyCones ( int printLevel );
}

#endif /* Mu2eG4_validPolyCones_hh */
