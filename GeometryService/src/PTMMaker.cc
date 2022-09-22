#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "cetlib_except/exception.h"

#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/PTMGeom/inc/PTM.hh"
#include "Offline/PTMGeom/inc/PTMHead.hh"
#include "Offline/PTMGeom/inc/PTMStand.hh"
#include "Offline/GeometryService/inc/PTMMaker.hh"
#include "Offline/GeomPrimitives/inc/Box.hh"
#include "Offline/GeomPrimitives/inc/ExtrudedSolid.hh"
//
// construct and return a PTM
//
// original author Helenka Casler
//

using namespace std;
namespace mu2e {

  std::unique_ptr<PTM> PTMMaker::make(SimpleConfig const& _config) {
    if (_config.getInt("PTM.version") == 1) {
        return makeFloatingPWCs(_config);
    } else if (_config.getInt("PTM.version") == 2) {
        return makeWithBasicStand(_config);
    } else {
        throw cet::exception("GEOM") << " illegal PTM version specified = " << _config.getInt("PTM.version")  << std::endl;
    }
    return 0;
  } // PTMMaker::make()

  std::unique_ptr<PTM> PTMMaker::makeWithBasicStand(SimpleConfig const& _config) {
    int version = _config.getInt("PTM.version");

    double xPosInMu2e = _config.getDouble("PTM.positionX");
    double yPosInMu2e = _config.getDouble("PTM.positionY");
    double zPosInMu2e = _config.getDouble("PTM.positionZ");
    CLHEP::Hep3Vector originInMu2e = CLHEP::Hep3Vector(xPosInMu2e, yPosInMu2e, zPosInMu2e);
    double ptmTotalLength = _config.getDouble("PTM.totalLength");
    double ptmTotalWidth = _config.getDouble("PTM.totalWidth");
    double ptmTotalHeight = _config.getDouble("PTM.totalHeight");

    double yRotInMu2e = _config.getDouble("PTM.rotY");
    double xRotInMu2e = _config.getDouble("PTM.rotX");
    CLHEP::HepRotation rotationInMu2e = CLHEP::HepRotation();
    rotationInMu2e.rotateX(xRotInMu2e*CLHEP::deg);
    rotationInMu2e.rotateY(yRotInMu2e*CLHEP::deg);

    // These values represent the center between the two detectors, not the
    // center of the head volume (which includes the handle). We will add the
    // handle's contribution later.
    double xPosInMu2eHead = _config.getDouble("PTM.head.positionX");
    double yPosInMu2eHead = _config.getDouble("PTM.head.positionY");
    double zPosInMu2eHead = _config.getDouble("PTM.head.positionZ");
    CLHEP::Hep3Vector headVolumeOriginInMu2e = CLHEP::Hep3Vector(xPosInMu2eHead, yPosInMu2eHead, zPosInMu2eHead);

    double yRotInMu2eHead = _config.getDouble("PTM.head.rotY");
    double xRotInMu2eHead = _config.getDouble("PTM.head.rotX");
    CLHEP::HepRotation headRotationInMu2e = CLHEP::HepRotation();
    headRotationInMu2e.rotateX(xRotInMu2eHead*CLHEP::deg);
    headRotationInMu2e.rotateY(yRotInMu2eHead*CLHEP::deg);

    // the PWCs:

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

    // The holder:

    std::string holderMaterialName = _config.getString("PTM.holder.extrusionMaterial");
    double holderExtrusionWidth = _config.getDouble("PTM.holder.extrusionWidth");
    double holderShortExtrusionPos = _config.getDouble("PTM.holder.supportDistFromCenter");
    double holderLongExtrusionSep = _config.getDouble("PTM.holder.separatorGap");

    // handle for RHS on top of the holder
    double handleHeight = _config.getDouble("PTM.handleHeight");
    double handleWidth = _config.getDouble("PTM.handleWidth");
    double handleThick = _config.getDouble("PTM.handleThick");
    double handleCornerCutoffs = _config.getDouble("PTM.handleCornerCutoffs");
    double handleCutoutHeight = _config.getDouble("PTM.handleCutoutHeight");
    double handleCutoutWidth = _config.getDouble("PTM.handleCutoutWidth");
    double handleCutoutCenterY = _config.getDouble("PTM.handleCutoutPositionY");
    std::string handleMaterial = _config.getString("PTM.handleMaterial");

    // mother volume contains both detectors, container contains one detector
    // each is slightly larger than all its contents
    double motherMargin = _config.getDouble("PTM.motherMargin");
    double containerMargin = _config.getDouble("PTM.containerMargin");

    // The stand for the PTM

    // top wedge that the PWCs are placed on
    double wedgeLength = _config.getDouble("PTM.wedge.length");
    double wedgeWidth = _config.getDouble("PTM.wedge.width");
    double wedgeMinHeight = _config.getDouble("PTM.wedge.minHeight");
    double wedgeMaxHeight = _config.getDouble("PTM.wedge.maxHeight");
    std::string wedgeMaterialName = _config.getString("PTM.wedge.materialName");
    // cutout for cables
    double wedgeCutoutLength = _config.getDouble("PTM.wedge.cutoutLength");
    double wedgeCutoutWidth = _config.getDouble("PTM.wedge.cutoutWidth");
    double wedgeCutoutRelX = _config.getDouble("PTM.wedge.cutoutPositionX");
    double wedgeCutoutRelZ = _config.getDouble("PTM.wedge.cutoutPositionZ");
    double wedgeCutoutRelRotY = _config.getDouble("PTM.wedge.cutoutRelRotY");
    CLHEP::Hep3Vector wedgeCutoutRelPosition = CLHEP::Hep3Vector(wedgeCutoutRelX, 0.0, wedgeCutoutRelZ);
    CLHEP::HepRotation wedgeCutoutRelRotation = CLHEP::HepRotation();
    wedgeCutoutRelRotation.rotateY(wedgeCutoutRelRotY*CLHEP::deg);
    double wedgeShiftDown = _config.getDouble("PTM.wedge.shiftDown");
    // aluminum extrusions that form a column, on which rests the top wedge
    double columnHeight = _config.getDouble("PTM.stand.mainColumnHeight");
    double columnExtrusionWidth = _config.getDouble("PTM.stand.extrusionWidth");
    std::string columnMaterialName = _config.getString("PTM.stand.extrusionMaterial");
    double columnExtrusionPositionAngle = _config.getDouble("PTM.stand.extrusionPositionAngle");
    double columnExtrusionRotation = _config.getDouble("PTM.stand.extrusionRotation");


    // First, the "head" portion: the PWC's and their holder

    // the PTM is made of two identical PWC's, placed such that their center
    // lines match up.
    // Y-position of the PWC's inside the mother volume, which also includes the holder and handle
    double headTotalHeight = holderLongExtrusionSep + handleHeight + motherMargin;
    double bottomOfHeadVolume = -0.5*headTotalHeight;
    double holderHeight = holderLongExtrusionSep + holderExtrusionWidth;
    double pwcY = bottomOfHeadVolume + (0.5*holderHeight);

    // "Near" PWC -- the more upstream of the two.
    CLHEP::Hep3Vector nearPWCPos = CLHEP::Hep3Vector(0.0, pwcY, -0.5*pwcSeparation);
    std::shared_ptr<PTMPWC> nearPWC( new PTMPWC("_1", frameHeight, frameWidth, frameThick, outerPlateThick, frameMaterialName,
                                     windowHeight, windowWidth, windowThick, windowMaterialName,
                                     gasMaterialName, numVertWires, numHorizWires, nearPWCPos, 0, containerMargin,
                                     framesInDetector, outerPlatesInDetector,
                                     ground1Zframes, hv1Zframes, hv2Zframes, hv3Zframes, ground2Zframes) );
    // "Far" PWC
    CLHEP::Hep3Vector farPWCPos = CLHEP::Hep3Vector(0.0, pwcY, 0.5*pwcSeparation);
    int farWireNumStart = numHorizWires + numVertWires;
    std::shared_ptr<PTMPWC> farPWC( new PTMPWC("_2", frameHeight, frameWidth, frameThick, outerPlateThick, frameMaterialName,
                                     windowHeight, windowWidth, windowThick, windowMaterialName,
                                     gasMaterialName, numVertWires, numHorizWires, farPWCPos, farWireNumStart, containerMargin,
                                     framesInDetector, outerPlatesInDetector,
                                     ground1Zframes, hv1Zframes, hv2Zframes, hv3Zframes, ground2Zframes) );

    // the holder consists of Al extrusions that hold the PWC's at the correct separation -- one in each corner.
    // It is supported partway down by some connecting extrusions of the same type
    double holderLength = pwcSeparation - 0.5*nearPWC->totalThick() - 0.5*farPWC->totalThick(); // separation is measured center-to-center
    std::shared_ptr<Box> holderExtrusionLong(new Box(0.5*holderExtrusionWidth, 0.5*holderExtrusionWidth, 0.5*holderLength));
    double shortExtrusionLength = 0.5*(holderLongExtrusionSep - holderExtrusionWidth);
    std::shared_ptr<Box> holderExtrusionShort(new Box(0.5*holderExtrusionWidth, 0.5*holderExtrusionWidth, shortExtrusionLength));

    // The handle is a rectagle with the top corners cut off, and an oval-shaped cutout
    std::shared_ptr<Box> handleBase(new Box(0.5*handleWidth, 0.5*handleHeight, 0.5*handleThick));
    double handleHoleSemiMajor = 0.5*handleCutoutHeight;
    double handleHoleSemiMinor = 0.5*handleCutoutWidth;
    CLHEP::Hep3Vector handleCutoutCenter = CLHEP::Hep3Vector(0.0, handleCutoutCenterY, 0.0);

    std::shared_ptr<PTMHead> ptmHead(new PTMHead(headVolumeOriginInMu2e,
                                                 headRotationInMu2e,
                                                 nearPWC,
                                                 farPWC,
                                                 pwcSeparation,
                                                 holderExtrusionLong,
                                                 holderExtrusionShort,
                                                 holderMaterialName,
                                                 holderLongExtrusionSep,
                                                 holderShortExtrusionPos,
                                                 handleBase,
                                                 handleMaterial,
                                                 handleHoleSemiMajor,
                                                 handleHoleSemiMinor,
                                                 handleCutoutCenter,
                                                 handleCornerCutoffs,
                                                 motherMargin,
                                                 originInMu2e,
                                                 rotationInMu2e));

    // build the top wedge in the stand
    // G4ExtrudedSolid assumes the vertices are given in (x, y) with the Z distance between given as dz
    // but our GeomPrimitives class returns it via getYhalfThickness()
    std::vector<CLHEP::Hep2Vector> wedgeVertices;
    // "bottom left" vertex
    wedgeVertices.push_back(CLHEP::Hep2Vector(-0.5*wedgeLength, -0.5*wedgeMaxHeight));
    // bottom right
    wedgeVertices.push_back(CLHEP::Hep2Vector(0.5*wedgeLength, -0.5*wedgeMaxHeight));
    // top right
    wedgeVertices.push_back(CLHEP::Hep2Vector(0.5*wedgeLength, wedgeMinHeight-(0.5*wedgeMaxHeight)));
    // top left
    wedgeVertices.push_back(CLHEP::Hep2Vector(-0.5*wedgeLength, 0.5*wedgeMaxHeight));

    // position of wedge is right under the "head"
    double wedgeY = headVolumeOriginInMu2e.y() - (0.5 * ptmHead->totalHeight()) - (0.5*wedgeMaxHeight) - (ptmHead->totalLength() * tan(xRotInMu2eHead*CLHEP::deg)) - wedgeShiftDown;
    CLHEP::Hep3Vector wedgeOriginInMu2e = CLHEP::Hep3Vector(headVolumeOriginInMu2e.x(), wedgeY, headVolumeOriginInMu2e.z());
    std::shared_ptr<ExtrudedSolid> topWedge(new ExtrudedSolid("PTMStandWedge", wedgeMaterialName, wedgeOriginInMu2e, 0.5*wedgeWidth, wedgeVertices));

    // build the support column that holds the wedge at the right height
    std::shared_ptr<Box> supportColExtrusion(new Box(0.5*columnExtrusionWidth, 0.5*columnHeight, 0.5*columnExtrusionWidth));
    // the three extrusions are arranged around a triangular central space
    double columnY = wedgeOriginInMu2e.y() - 0.5*wedgeMaxHeight - 0.5*columnHeight;
    CLHEP::Hep3Vector columnLocalCenter = CLHEP::Hep3Vector(wedgeOriginInMu2e.x(), columnY, wedgeOriginInMu2e.z());
    double distFromLocalCenter = 0.5*columnExtrusionWidth*tan(30*CLHEP::deg) + 0.5*columnExtrusionWidth;
    CLHEP::Hep3Vector columnExtrusionPos1 = CLHEP::Hep3Vector(0.0, 0.0, distFromLocalCenter) + columnLocalCenter;
    CLHEP::Hep3Vector columnExtrusionPos2 = CLHEP::Hep3Vector(0.0, 0.0, distFromLocalCenter);
    columnExtrusionPos2.rotateY(columnExtrusionPositionAngle*CLHEP::deg);
    columnExtrusionPos2 += columnLocalCenter;
    CLHEP::Hep3Vector columnExtrusionPos3 = CLHEP::Hep3Vector(0.0, 0.0, distFromLocalCenter);
    columnExtrusionPos3.rotateY(2*columnExtrusionPositionAngle*CLHEP::deg);
    columnExtrusionPos3 += columnLocalCenter;
    std::vector<CLHEP::Hep3Vector> columnPositions;
    columnPositions.push_back(columnExtrusionPos1);
    columnPositions.push_back(columnExtrusionPos2);
    columnPositions.push_back(columnExtrusionPos3);
    // rotations of the columns, so they're all "facing" the center
    CLHEP::HepRotation columnRotation1 = CLHEP::HepRotation();
    CLHEP::HepRotation columnRotation2 = CLHEP::HepRotation();
    columnRotation2.rotateY(columnExtrusionRotation *CLHEP::deg);
    CLHEP::HepRotation columnRotation3 = CLHEP::HepRotation();
    columnRotation3.rotateY(-1*columnExtrusionRotation *CLHEP::deg);
    std::vector<CLHEP::HepRotation> columnRotations;
    columnRotations.push_back(columnRotation1);
    columnRotations.push_back(columnRotation2);
    columnRotations.push_back(columnRotation3);



    std::shared_ptr<PTMStand> ptmStand(new PTMStand(topWedge,
          wedgeCutoutLength,
          0.5*wedgeCutoutWidth,
          wedgeCutoutRelPosition,
          wedgeCutoutRelRotation,
          supportColExtrusion,
          columnPositions,
          columnRotations,
          columnMaterialName,
          wedgeMaterialName));
    std::unique_ptr<PTM> ptmon(new PTM(version, originInMu2e, rotationInMu2e, ptmStand, ptmHead, ptmTotalLength, ptmTotalWidth, ptmTotalHeight));
    return ptmon;
  } // PTMMaker::makeWithBasicStand()

std::unique_ptr<PTM> PTMMaker::makeFloatingPWCs(SimpleConfig const& _config) {
    int version = _config.getInt("PTM.version");

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

    std::unique_ptr<PTM> ptmon(new PTM(version, originInMu2e, rotationInMu2e, nearPWC, farPWC, pwcSeparation, motherMargin));
    return ptmon;
  } // PTMMaker::makeFloatingPWCs()

} // namespace mu2e
