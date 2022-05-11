#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/PTMGeom/inc/PTM.hh"
#include "Offline/GeometryService/inc/PTMMaker.hh"
//
// construct and return a PTM
//
// original author Helenka Casler
//

using namespace std;
namespace mu2e {
  std::unique_ptr<PTM> PTMMaker::make(SimpleConfig const& _config) {
    double xPosInMu2e = _config.getDouble("PTM.positionX");
    double yPosInMu2e = _config.getDouble("PTM.positionY");
    double zPosInMu2e = _config.getDouble("PTM.positionZ");
    CLHEP::Hep3Vector originInMu2e = CLHEP::Hep3Vector(xPosInMu2e, yPosInMu2e, zPosInMu2e);

    double yRotInMu2e = _config.getDouble("PTM.rotY");
    double xRotInMu2e = _config.getDouble("PTM.rotX");
    CLHEP::HepRotation rotationInMu2e = CLHEP::HepRotation();
    rotationInMu2e.rotateX(xRotInMu2e*CLHEP::deg);
    rotationInMu2e.rotateY(yRotInMu2e*CLHEP::deg);

    double pwcSeparation = _config.getDouble("PTM.pwcSeparation");
    double frameHeight = _config.getDouble("PTM.frameHeight");
    double frameWidth = _config.getDouble("PTM.frameWidth");
    double frameThick = _config.getDouble("PTM.frameThick");
    double outerPlateThick = _config.getDouble("PTM.outerPlateLength");
    std::string frameMaterialName = _config.getString("PTM.frameMaterial");

    int framesInDetector = _config.getInt("PTM.framesInDetector");
    int outerPlatesInDetector = _config.getInt("PTM.outerPlatesInDetector");
    double ground1Zframes = _config.getDouble("PTM.ground1Zframes");
    double hv1Zframes = _config.getDouble("PTM.hv1Zframes");
    double hv2Zframes = _config.getDouble("PTM.hv2Zframes");
    double hv3Zframes = _config.getDouble("PTM.hv3Zframes");
    double ground2Zframes = _config.getDouble("PTM.ground2Zframes");

    double windowHeight = _config.getDouble("PTM.windowHeight");
    double windowWidth = _config.getDouble("PTM.windowWidth");
    double windowThick = _config.getDouble("PTM.windowThick");
    std::string windowMaterialName = _config.getString("PTM.windowMaterial");

    std::string gasMaterialName = _config.getString("PTM.innerGas");
    int numVertWires = _config.getInt("PTM.vertWiresPerPlane");
    int numHorizWires = _config.getInt("PTM.horizWiresPerPlane");

    // mother volume contains both detectors, container contains one detector
    // each is slightly larger than all its contents
    double motherMargin = _config.getDouble("PTM.motherMargin");
    double containerMargin = _config.getDouble("PTM.containerMargin");

    // the PTM is made of two identical PWC's, placed such that their center
    // lines match up.

    // "Near" PWC -- the more upstream of the two.
    CLHEP::Hep3Vector nearPWCPos = CLHEP::Hep3Vector(0.0, 0.0, -0.5*pwcSeparation);
    std::shared_ptr<PTMPWC> nearPWC( new PTMPWC("_1", frameHeight, frameWidth, frameThick, outerPlateThick, frameMaterialName,
                                     windowHeight, windowWidth, windowThick, windowMaterialName,
                                     gasMaterialName, numVertWires, numHorizWires, nearPWCPos, 0, containerMargin,
                                     framesInDetector, outerPlatesInDetector,
                                     ground1Zframes, hv1Zframes, hv2Zframes, hv3Zframes, ground2Zframes) );
    // "Far" PWC
    CLHEP::Hep3Vector farPWCPos = CLHEP::Hep3Vector(0.0, 0.0, 0.5*pwcSeparation);
    int farWireNumStart = numHorizWires + numVertWires;
    std::shared_ptr<PTMPWC> farPWC( new PTMPWC("_2", frameHeight, frameWidth, frameThick, outerPlateThick, frameMaterialName,
                                     windowHeight, windowWidth, windowThick, windowMaterialName,
                                     gasMaterialName, numVertWires, numHorizWires, farPWCPos, farWireNumStart, containerMargin,
                                     framesInDetector, outerPlatesInDetector,
                                     ground1Zframes, hv1Zframes, hv2Zframes, hv3Zframes, ground2Zframes) );

    std::unique_ptr<PTM> ptmon(new PTM(originInMu2e, rotationInMu2e, nearPWC, farPWC, pwcSeparation, motherMargin));
    return ptmon;
  } // PTMMaker::make()

} // namespace mu2e
