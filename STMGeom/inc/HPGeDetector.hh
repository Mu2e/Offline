#ifndef STMGeom_HPGeDetector_hh
#define STMGeom_HPGeDetector_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class HPGeDetector {
  public:
    HPGeDetector(bool build,
               std::string const & crystalMaterial,
               std::string const & holeMaterial,
               std::string const & windowMaterial,
               std::string const & wallMaterial,
               std::string const & capsuleMaterial,
               double EndcapR, double EndcapL,
               double CrystalR, double CrystalL,
               double Z_HPGe, double HoleR, double HoleL,
               double Capsule_Wallthick, double Capsule_Windowthick,
               double Capsule_Endthick, double Capsule_Walllength,
               double WindowD, double EndcapD, double AirD, double offset_HPGe,
               CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
               CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
              ):
      _build(build),
      _crystalMaterial(crystalMaterial),
      _holeMaterial(holeMaterial),
      _windowMaterial(windowMaterial),
      _wallMaterial(wallMaterial),
      _capsuleMaterial(capsuleMaterial),
      _EndcapR(EndcapR), _EndcapL(EndcapL),
      _CrystalR(CrystalR), _CrystalL(CrystalL),
      _Z_HPGe(Z_HPGe), _HoleR(HoleR), _HoleL(HoleL),
      _Capsule_Wallthick(Capsule_Wallthick), _Capsule_Windowthick(Capsule_Windowthick),
      _Capsule_Endthick(Capsule_Endthick), _Capsule_Walllength(Capsule_Walllength),
      _WindowD(WindowD), _EndcapD(EndcapD), _AirD(AirD), _offset_HPGe(offset_HPGe),
      _originInMu2e(originInMu2e),
      _rotation(rotation)
    {}

    bool   build()                               const {return _build;}
    std::string const & crystalMaterial()        const {return _crystalMaterial;}
    std::string const & holeMaterial()           const {return _holeMaterial;}
    std::string const & windowMaterial()         const {return _windowMaterial;}
    std::string const & wallMaterial()           const {return _wallMaterial;}
    std::string const & capsuleMaterial()        const {return _capsuleMaterial;}

    double EndcapR()                    const {return _EndcapR;}
    double EndcapL()                    const {return _EndcapL;}
    double CrystalR()                   const {return _CrystalR;}
    double CrystalL()                   const {return _CrystalL;}
    double Z_HPGe()                          const {return _Z_HPGe;}
    double HoleR()                           const {return _HoleR;}
    double HoleL()                           const {return _HoleL;}
    double Capsule_Wallthick()               const {return _Capsule_Wallthick;}
    double Capsule_Windowthick()             const {return _Capsule_Windowthick;}
    double Capsule_Endthick()                const {return _Capsule_Endthick;}
    double Capsule_Walllength()              const {return _Capsule_Walllength;}
    double WindowD()                         const {return _WindowD;}
    double EndcapD()                         const {return _EndcapD;}
    double AirD()                            const {return _AirD;}
    double offset_HPGe()                     const {return _offset_HPGe;}


    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    //double zBegin()          const { return _originInMu2e.z() - zHalfLength(); }
    //double zEnd()            const { return _originInMu2e.z() + zHalfLength(); }


    // Genreflex can't do persistency of vector<HPGeDetector> without a default constructor
    HPGeDetector() {}
  private:

    bool               _build;
    std::string        _crystalMaterial;
    std::string        _holeMaterial;
    std::string        _windowMaterial;
    std::string        _wallMaterial;
    std::string        _capsuleMaterial;
    double _EndcapR;
    double _EndcapL;
    double _CrystalR;
    double _CrystalL;
    double _Z_HPGe;
    double _HoleR;
    double _HoleL;
    double _Capsule_Wallthick;
    double _Capsule_Windowthick;
    double _Capsule_Endthick;
    double _Capsule_Walllength;
    double _WindowD;
    double _EndcapD;
    double _AirD;
    double _offset_HPGe;
    
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

}

#endif/*STMGeom_HPGeDetector_hh*/
