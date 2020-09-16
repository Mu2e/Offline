#ifndef Mu2eG4_Mu2eGlobalField_hh
#define Mu2eG4_Mu2eGlobalField_hh
//
// G4 interface to the Detector Solenoid full magnetic field.
//
//
// Original author Julie Managan and Bob Bernstein
// Major rewrite Rob Kutschke at version 1.4
//

#include <string>

#include "BFieldGeom/inc/BFCacheManager.hh"

#include "G4MagneticField.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

namespace mu2e {

  // Forward references
  class BFieldManager;

  class Mu2eGlobalField: public G4MagneticField {

  public:

    explicit Mu2eGlobalField(const G4ThreeVector& mapOrigin);
    virtual ~Mu2eGlobalField(){}

    // This is called by G4.
    virtual void GetFieldValue(const G4double Point[4],
                               G4double *Bfield) const;

    // Update the map and its origin.  Must be called whenever
    // the map or the offset changes (begin run probably).
    void update( const G4ThreeVector& mapOrigin );

  private:
    // The map is stored in the Mu2e coordinate system.
    // This is the location of the origin the Mu2e system, measured in the G4 world system.
    G4ThreeVector _mapOrigin;

    // Non-owning pointer to the field map object (it is owned by the geometry service).
    const BFieldManager* _map;

    // A copy of the bfield cache manager - must be thread local.
    BFCacheManager _cm;

  };
}
#endif /* Mu2eG4_Mu2eGlobalField_hh */
