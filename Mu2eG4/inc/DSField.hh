#ifndef Mu2eG4_DSField_hh
#define Mu2eG4_DSField_hh
//
// G4 interface to the Detector Solenoid full magnetic field.
//
// $Id: DSField.hh,v 1.6 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Julie Managan and Bob Bernstein
// Major rewrite Rob Kutschke at version 1.4
//

#include <string>

#include "G4MagneticField.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

namespace mu2e {

  // Forward references
  class BFMapBase;

  class DSField: public G4MagneticField {

  public:

    DSField( std::string name, G4ThreeVector mapOrigin );
    virtual ~DSField(){};

    // This is called by G4.
    virtual void GetFieldValue(const G4double Point[4],
                               G4double *Bfield) const;

    const std::string& name() const { return _name; }

    // Update the map and its origin.  Must be called whenever
    // the map or the offset changes (begin run probably).
    void update( const G4ThreeVector& mapOrigin );

  private:

    // Name of this map.
    std::string _name;

    // The map is stored in the Mu2e coordinate system.
    // This is the location of the origin the Mu2e system, measured in the G4 world system.
    G4ThreeVector _mapOrigin;

    // Non-owning pointer to the field map object (it is owned by the geometry service).
    const BFMapBase* _map;

  };
}
#endif /* Mu2eG4_DSField_hh */
