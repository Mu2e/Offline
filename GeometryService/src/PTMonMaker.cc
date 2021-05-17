#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "PTMonGeom/inc/PTMonPWC.hh"
#include "PTMonGeom/inc/PTMon.hh"
#include "GeometryService/inc/PTMonMaker.hh"
//
// construct and return a PTMon
//
// original author Helenka Casler
//

using namespace std;
namespace mu2e {
  std::unique_ptr<PTMon> PTMonMaker::make(SimpleConfig const& _config) {
    double xPosInMu2e = _config.getDouble("pTargetMon_positionX");
    double yPosInMu2e = _config.getDouble("pTargetMon_positionY");
    double zPosInMu2e = _config.getDouble("pTargetMon_positionZ");
    CLHEP::Hep3Vector originInMu2e = CLHEP::Hep3Vector(xPosInMu2e, yPosInMu2e, zPosInMu2e);

    double yRotInMu2e = _config.getDouble("pTargetMon_rotY");
    double xRotInMu2e = _config.getDouble("pTargetMon_rotX");
    CLHEP::HepRotation rotationInMu2e = CLHEP::HepRotation();
    rotationInMu2e.rotateX(xRotInMu2e*CLHEP::deg);
    rotationInMu2e.rotateY(yRotInMu2e*CLHEP::deg);

    double pwcSeparation = _config.getDouble("pTargetMon_pwcSeparation");
    double frameHeight = _config.getDouble("pTargetMon_frameHeight");
    double frameWidth = _config.getDouble("pTargetMon_frameWidth");
    double frameThick = _config.getDouble("pTargetMon_frameThick");
    double outerPlateThick = _config.getDouble("pTargetMon_outerPlateLength");
    std::string frameMaterialName = _config.getString("pTargetMon_frameMaterial");

    double windowHeight = _config.getDouble("pTargetMon_windowHeight");
    double windowWidth = _config.getDouble("pTargetMon_windowWidth");
    double windowThick = _config.getDouble("pTargetMon_windowThick");
    std::string windowMaterialName = _config.getString("pTargetMon_windowMaterial");

    std::string gasMaterialName = _config.getString("pTargetMon_innerGas");
    int numVertWires = _config.getInt("pTargetMon_vertWiresPerPlane");
    int numHorizWires = _config.getInt("pTargetMon_horizWiresPerPlane");

    // mother volume contains both detectors, container contains one detector
    // each is slightly larger than all its contents
    double motherMargin = _config.getDouble("pTargetMon_motherMargin");
    double containerMargin = _config.getDouble("pTargetMon_containerMargin");

    // the PTMon is made of two identical PWC's, placed such that their center
    // lines match up.

    // "Near" PWC -- the more upstream of the two.
    CLHEP::Hep3Vector nearPWCPos = CLHEP::Hep3Vector(0.0, 0.0, -0.5*pwcSeparation);
    std::shared_ptr<PTMonPWC> nearPWC( new PTMonPWC("_1", frameHeight, frameWidth, frameThick, outerPlateThick, frameMaterialName, 
                                     windowHeight, windowWidth, windowThick, windowMaterialName,
                                     gasMaterialName, numVertWires, numHorizWires, nearPWCPos, 0, containerMargin) );
    // "Far" PWC
    CLHEP::Hep3Vector farPWCPos = CLHEP::Hep3Vector(0.0, 0.0, 0.5*pwcSeparation);
    int farWireNumStart = numHorizWires + numVertWires;
    std::shared_ptr<PTMonPWC> farPWC( new PTMonPWC("_2", frameHeight, frameWidth, frameThick, outerPlateThick, frameMaterialName, 
                                     windowHeight, windowWidth, windowThick, windowMaterialName,
                                     gasMaterialName, numVertWires, numHorizWires, farPWCPos, farWireNumStart, containerMargin) );

    std::unique_ptr<PTMon> ptmon(new PTMon(originInMu2e, rotationInMu2e, nearPWC, farPWC, pwcSeparation, motherMargin));
    return ptmon;
  } // PTMonMaker::make()

} // namespace mu2e