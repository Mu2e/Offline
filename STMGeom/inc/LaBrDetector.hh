#ifndef STMGeom_LaBrDetector_hh
#define STMGeom_LaBrDetector_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LaBrDetector {
  public:
    LaBrDetector(bool build,
               std::string const & crystalMaterial,
               std::string const & windowMaterial,
               std::string const & wallMaterial,
               double EndcapR, double EndcapL,
               double CrystalR, double CrystalL,
               double Z_LaBr,
               double WindowD, double EndcapD, double AirD, double offset_LaBr,
               CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
               CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
              ):
      _build(build),
      _crystalMaterial(crystalMaterial),
      _windowMaterial(windowMaterial),
      _wallMaterial(wallMaterial),
      _EndcapR(EndcapR), _EndcapL(EndcapL),
      _CrystalR(CrystalR), _CrystalL(CrystalL),
      _Z_LaBr(Z_LaBr),
      _WindowD(WindowD), _EndcapD(EndcapD), _AirD(AirD),
      _originInMu2e(originInMu2e),
      _rotation(rotation)
    {}

    bool   build()                               const {return _build;}
    std::string const & crystalMaterial()        const {return _crystalMaterial;}
    std::string const & windowMaterial()         const {return _windowMaterial;}
    std::string const & wallMaterial()           const {return _wallMaterial;}

    double EndcapR()                    const {return _EndcapR;}
    double EndcapL()                    const {return _EndcapL;}
    double CrystalR()                   const {return _CrystalR;}
    double CrystalL()                   const {return _CrystalL;}
    double Z_LaBr()                          const {return _Z_LaBr;}
    double WindowD()                         const {return _WindowD;}
    double EndcapD()                         const {return _EndcapD;}
    double AirD()                            const {return _AirD;}
    double offset_LaBr()                     const {return _offset_LaBr;}


    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    //double zBegin()          const { return _originInMu2e.z() - zHalfLength(); }
    //double zEnd()            const { return _originInMu2e.z() + zHalfLength(); }


    // Genreflex can't do persistency of vector<LaBrDetector> without a default constructor
    LaBrDetector() {}
  private:

    bool               _build;
    std::string        _crystalMaterial;
    std::string        _windowMaterial;
    std::string        _wallMaterial;

    double _EndcapR;
    double _EndcapL;
    double _CrystalR;
    double _CrystalL;
    double _Z_LaBr;
    double _WindowD;
    double _EndcapD;
    double _AirD;
    double _offset_LaBr;
    
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

}

#endif/*STMGeom_LaBrDetector_hh*/
