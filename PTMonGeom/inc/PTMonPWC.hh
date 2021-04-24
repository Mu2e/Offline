#ifndef PTMonGeom_PTMonPWC_hh
#define PTMonGeom_PTMonPWC_hh

// C++ includes
#include <vector>
#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "GeomPrimitives/inc/Box.hh"

// Proportional Wire Chamber object, part of the 
// production target monitor (PTMon)
//
// Author: Helenka Casler
//

namespace mu2e {

  class PTMonPWC {
  public:
    PTMonPWC(std::string nameSuffix,
             double frameHeight, 
             double frameWidth, 
             double frameThick, 
             double outerPlateThick,
             std::string frameMaterialName, 
             double windowHeight,
             double windowWidth,
             double windowThick,
             std::string windowMaterialName,
             std::string gasMaterialName,
             int numVertWires,
             int numHorizWires,
             CLHEP::Hep3Vector const & originInParent = CLHEP::Hep3Vector(),
             CLHEP::HepRotation const & rotationInParent = CLHEP::HepRotation(),
             int wireNumStart = 0);
    // default ctor
    PTMonPWC() {}

    CLHEP::Hep3Vector const &  originInParent()   const { return _originInParent; }
    CLHEP::HepRotation const & rotationInParent() const { return _rotationInParent; }
    std::vector<double>        vertWireYPos()     const { return _vertWireYpos; }
    std::vector<double>        horizWireXPos()    const { return _horizWireXpos; }

    const Box*                 pwcWindow()           const { return _pwcWindow.get(); }
    const Box*                 vertWireGasSection()  const { return _vertWireGasSection.get(); }
    const Box*                 horizWireGasSection() const { return _horizWireGasSection.get(); }
    const Box*                 gasSection1()         const { return _gasSection1.get(); }
    const Box*                 gasSection4()         const { return _gasSection4.get(); }

    std::string nameSuffix()         const { return _nameSuffix; }
    double      frameHeight()        const { return _frameHeight; }
    double      frameWidth()         const { return _frameWidth; }
    double      frameThick()         const { return _frameThick; }
    double      outerPlateThick()    const { return _outerPlateThick; }
    std::string frameMaterialName()  const { return _frameMaterialName; }
    std::string windowMaterialName() const { return _windowMaterialName; }
    std::string gasMaterialName()    const { return _gasMaterialName; }
    int         wireNumStart()       const { return _wireNumStart; }
    int         numVertWires()       const { return _numVertWires; }
    double      vertWireZ()          const { return _vertWireZ; }
    int         numHorizWires()      const { return _numHorizWires; }
    double      horizWireZ()         const { return _horizWireZ; }
    double      ground1Z()           const { return _ground1Z; }
    double      gasInZ()             const { return _gasInZ; }
    double      hv1Z()               const { return _hv1Z; }
    double      hv2Z()               const { return _hv2Z; }
    double      hv3Z()               const { return _hv3Z; }
    double      gasOutZ()            const { return _gasOutZ; }
    double      ground2Z()           const { return _ground2Z; }
    double      totalThick()         const {return _totalThick; }
    


  private:
    std::string _nameSuffix;
    CLHEP::Hep3Vector _originInParent;
    CLHEP::HepRotation _rotationInParent;
    double _frameHeight;
    double _frameWidth;
    double _frameThick;
    double _outerPlateThick;
    std::string _frameMaterialName;

    std::unique_ptr<Box> _pwcWindow;
    std::string _windowMaterialName;

    std::string _gasMaterialName;
    int _wireNumStart;

    // "Vert" wires have their long dimension in the horizontal.
    int _numVertWires; // MEASURES the vertical profile!
    std::unique_ptr<Box> _vertWireGasSection;
    std::vector<double> _vertWireYpos;
    double _vertWireZ;
    // "Horiz" wires have their long dimension in the vertical.
    int _numHorizWires; // MEASURES the horizontal profile!
    std::unique_ptr<Box> _horizWireGasSection;
    std::vector<double> _horizWireXpos;
    double _horizWireZ;

    // numbering below is becase gas sections 2 and 3 contain wires
    std::unique_ptr<Box> _gasSection1;
    std::unique_ptr<Box> _gasSection4;

    // positions of windows inside detector
    double _ground1Z;
    double _hv1Z;
    double _hv2Z;
    double _hv3Z;
    double _ground2Z;
    double _gasInZ;
    double _gasOutZ;

    double _totalThick;

  }; // class PTMonPWC
} // namespace mu2e


#endif