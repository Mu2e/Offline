#ifndef Mu2eG4_Mu2eG4DSGradientMagneticField_hh
#define Mu2eG4_Mu2eG4DSGradientMagneticField_hh
//
// G4 interface to the Detector Solenoid gradient magnetic field.
// Right now the field direction is assumed to be along Z.
// If neccessary, the class can easily be extended to allow any direction,
// but it is unlikely to be ever needed.
//
// Author Ivan Logashenko
//

#include <string>

#include "Geant4/G4MagneticField.hh"
#include "Geant4/G4Types.hh"
#include "Geant4/G4ThreeVector.hh"

namespace mu2e {

  class Mu2eG4DSGradientMagneticField: public G4MagneticField {

  public:

    Mu2eG4DSGradientMagneticField( std::string name, G4ThreeVector mapOrigin, double gradient, double B0 );
    virtual ~Mu2eG4DSGradientMagneticField(){}

    // This is called by G4.
    virtual void GetFieldValue(const G4double Point[4],
                               G4double *Bfield) const;

    const std::string& name() const { return _name; }

  private:

    // Name of this map.
    std::string _name;

    // Location of the field "origin" - at this point the field is Bz=B0, Br=0
    G4ThreeVector _mapOrigin;

    // Field parameters
    double _gradient;
    double _B0;

  };
}
#endif /* Mu2eG4_Mu2eG4DSGradientMagneticField_hh */
