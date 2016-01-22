#ifndef STMGeom_GeDetector_hh
#define STMGeom_GeDetector_hh

// Germanium Detector Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class GeDetector {
  public:
    GeDetector(bool build, 
               std::string const & crystalMaterial,
               double crystalRadiusIn, double crystalRadiusOut, double crystalHalfLength,
               std::string const & canMaterial,
               double canRadiusIn, double canRadiusOut, double canHalfLength,
               std::string const & canUpStrWindowMaterial, double canUpStrWindowHalfLength,
               std::string const & canGasMaterial,
               CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
               CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
              ) :
      _build( build ),          
      _crystalMaterial( crystalMaterial ),
      _crystalRadiusIn( crystalRadiusIn ),
      _crystalRadiusOut( crystalRadiusOut ),
      _crystalHalfLength( crystalHalfLength ),
      _canMaterial( canMaterial ),
      _canRadiusIn( canRadiusIn ),
      _canRadiusOut( canRadiusOut ),
      _canHalfLength( canHalfLength ),
      _canUpStrWindowMaterial( canUpStrWindowMaterial ),
      _canUpStrWindowHalfLength( canUpStrWindowHalfLength ),
      _canGasMaterial( canGasMaterial ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     )
    {}

    bool   build()                               const { return _build;  }
    std::string const & crystalMaterial()        const { return _crystalMaterial; }    
    double crystalRadiusIn()                     const { return _crystalRadiusIn;  }
    double crystalRadiusOut()                    const { return _crystalRadiusOut;  }
    double crystalHalfLength()                   const { return _crystalHalfLength;  }
    std::string const & canMaterial()            const { return _canMaterial; }    
    double canRadiusIn()                         const { return _canRadiusIn;  }
    double canRadiusOut()                        const { return _canRadiusOut;  }
    double canHalfLength()                       const { return _canHalfLength;  }
    std::string const & canUpStrWindowMaterial() const { return _canUpStrWindowMaterial; }    
    double canUpStrWindowHalfLength()            const { return _canUpStrWindowHalfLength; }    
    std::string const & canGasMaterial()         const { return _canGasMaterial; }    
    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    //double zBegin()          const { return _originInMu2e.z() - zHalfLength(); }
    //double zEnd()            const { return _originInMu2e.z() + zHalfLength(); }


    // Genreflex can't do persistency of vector<GeDetector> without a default constructor
    GeDetector() {}
  private:

    bool               _build;
    std::string        _crystalMaterial;
    double             _crystalRadiusIn;
    double             _crystalRadiusOut;
    double             _crystalHalfLength;
    std::string        _canMaterial;
    double             _canRadiusIn;
    double             _canRadiusOut;
    double             _canHalfLength;
    double             _canUpStrSpace;
    std::string        _canUpStrWindowMaterial;
    double             _canUpStrWindowHalfLength;
    std::string        _canGasMaterial;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

}

#endif/*STMGeom_GeDetector_hh*/
