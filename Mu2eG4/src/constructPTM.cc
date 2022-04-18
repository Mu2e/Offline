//
// Free function. Approach borrowed from constructPS
// Constructs the downstream production target scanning monitor.
// Parent volume is the air in the target hall.
//

// C++ includes
#include <iostream>
#include <vector>
#include <string>

// Mu2e includes
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"
#include "Offline/PTMGeom/inc/PTM.hh"
#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

// G4 inludes
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4VSolid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"


using namespace std;

namespace mu2e {

  void insertOuterFrame(VolumeInfo const& container,
                        const PTMPWC* pwc) {
    G4Box *outerBox = new G4Box("pwcFrameOuter",
                pwc->frameWidth()/2.,
                pwc->frameHeight()/2.,
                pwc->detectorThick()/2.);
    G4Box *innerBox = new G4Box("pwcFrameInner",
                pwc->pwcWindow()->getXhalfLength(),
                pwc->pwcWindow()->getYhalfLength(),
                pwc->totalThick()/2.);
    std::string frameName = "PTMFrame";
    frameName.append(pwc->nameSuffix());
    G4Material *frameMaterial = findMaterialOrThrow(pwc->frameMaterialName());

    VolumeInfo frameInfo;
    frameInfo.name = frameName;
    frameInfo.solid = new G4SubtractionSolid(frameName, outerBox, innerBox);
    finishNesting(frameInfo,
          frameMaterial,
          nullptr,
          G4ThreeVector(0.0, 0.0, 0.0),
          container.logical,
          0,
          G4Colour::Blue(),
          "PTM");

  } // insertOuterFrame

  void insertWindows(VolumeInfo const& container,
                     const PTMPWC* pwc,
                     int const vdNum,
                     SimpleConfig const& _config,
                     bool const visible,
                     bool const forceSolid,
                     bool const forceAuxEdgeVisible,
                     bool const placePV,
                     bool const doSurfaceCheck,
                     int const verbosity) {
    G4Material *windowMaterial = findMaterialOrThrow(pwc->windowMaterialName());
    // most-upstream ground plane is a virtual detector
    std::string ground1Name = "VirtualDetector_PTMGroundIn"+pwc->nameSuffix();
    double windowHalfDims[] = {pwc->pwcWindow()->getXhalfLength(),
                               pwc->pwcWindow()->getYhalfLength(),
                               pwc->pwcWindow()->getZhalfLength()};
    nestBox (ground1Name,
             windowHalfDims,
             windowMaterial,
             nullptr,
             G4ThreeVector(0.0, 0.0, pwc->ground1Z()),
             container,
             vdNum, // copyNo
             G4Colour::Blue(),
             "PTM"); // lookup token for doSurfaceCheck, etc
    // The remaining windows are all identical, so make one logical volume and
    // paste it repeatedly.
    std::string windowName = "PTMWindow"+pwc->nameSuffix();
    G4VSolid *windowBox = new G4Box(windowName,
                            pwc->pwcWindow()->getXhalfLength(),
                            pwc->pwcWindow()->getYhalfLength(),
                            pwc->pwcWindow()->getZhalfLength());
    G4LogicalVolume *windowLogical = new G4LogicalVolume(windowBox, windowMaterial, windowName);
    // 3 high-voltage planes
    VolumeInfo HV1Info;
    HV1Info.solid = windowBox;
    HV1Info.logical = windowLogical;
    HV1Info.name = "PTMHV1"+pwc->nameSuffix();
    finishNesting(HV1Info,
                  windowMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(0, 0, pwc->hv1Z()),
                  container.logical,
                  0,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    VolumeInfo HV2Info;
    HV2Info.solid = windowBox;
    HV2Info.logical = windowLogical;
    HV2Info.name = "PTMHV2"+pwc->nameSuffix();
    finishNesting(HV2Info,
                  windowMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(0, 0, pwc->hv2Z()),
                  container.logical,
                  0,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    VolumeInfo HV3Info;
    HV3Info.solid = windowBox;
    HV3Info.logical = windowLogical;
    HV3Info.name = "PTMHV3"+pwc->nameSuffix();
    finishNesting(HV3Info,
                  windowMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(0, 0, pwc->hv3Z()),
                  container.logical,
                  0,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    // downstream ground plane
    VolumeInfo ground2Info;
    ground2Info.solid = windowBox;
    ground2Info.logical = windowLogical;
    ground2Info.name = "PTMGround2"+pwc->nameSuffix();
    finishNesting(ground2Info,
                  windowMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(0, 0, pwc->ground2Z()),
                  container.logical,
                  0,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);

  } // insertWindows

  void insertOuterGasBlocks(VolumeInfo const& container,
                            const PTMPWC* pwc,
                            G4Material* gasMaterial) {
    // between ground plane 1 and HV plane 1
    std::string gasName1 = "PTMGas1";
    gasName1.append(pwc->nameSuffix());
    std::vector<double> gas1HalfDims;
    gas1HalfDims.push_back(pwc->gasSection1()->getXhalfLength());
    gas1HalfDims.push_back(pwc->gasSection1()->getYhalfLength());
    gas1HalfDims.push_back(pwc->gasSection1()->getZhalfLength());
    nestBox (gasName1,
             gas1HalfDims,
             gasMaterial,
             nullptr,
             G4ThreeVector(0.0, 0.0, pwc->gasInZ()),
             container,
             0, // copyNo
             G4Colour::Red(),
             "PTM"); // lookup token for doSurfaceCheck, etc


    // between HV plane 3 and ground plane 2
    std::string gasName4 = "PTMGas4";
    gasName4.append(pwc->nameSuffix());
    std::vector<double> gas4HalfDims;
    gas4HalfDims.push_back(pwc->gasSection4()->getXhalfLength());
    gas4HalfDims.push_back(pwc->gasSection4()->getYhalfLength());
    gas4HalfDims.push_back(pwc->gasSection4()->getZhalfLength());
    nestBox (gasName4,
             gas4HalfDims,
             gasMaterial,
             nullptr,
             G4ThreeVector(0.0, 0.0, pwc->gasOutZ()),
             container,
             0, // copyNo
             G4Colour::Red(),
             "PTM");

  } // insertOuterGasBlocks

  void insertVerticalProfileWires(VolumeInfo const& container,
                                  const PTMPWC* pwc,
                                  G4Material* gasMaterial,
                                  std::string const& wireNameSuffix,
                                  SimpleConfig const& _config,
                                  bool const visible,
                                  bool const forceSolid,
                                  bool const forceAuxEdgeVisible,
                                  bool const placePV,
                                  bool const doSurfaceCheck,
                                  int const verbosity) {
    // In the real detector wires run HORIZONTALLY, so as to measure the
    // VERTICAL profile.
    // G4 geometry currently contains no actual wires.
    // Contains gas broken into sections -- each section represents the region
    // of the gas that is closest to one wire in the real detector.
    std::string wireGasNameVert = "PTMWireVert";
    wireGasNameVert.append(wireNameSuffix);
    G4VSolid *vertWireBox = new G4Box(wireGasNameVert,
                                    pwc->vertWireGasSection()->getXhalfLength(),
                                    pwc->vertWireGasSection()->getYhalfLength(),
                                    pwc->vertWireGasSection()->getZhalfLength());
    G4LogicalVolume *vertWireLogical = new G4LogicalVolume(vertWireBox, gasMaterial, wireGasNameVert);
    for (int i=0; i < pwc->numVertWires(); ++i) {
      int wireNum = pwc->wireNumStart() + i;
      std::string wireGasName = wireGasNameVert;
      wireGasName.append(std::to_string(wireNum));
      // wire numbering such that the lowest-numered wire is
      // on the bottom
      double gasY2 = pwc->vertWireYPos()[i];
      VolumeInfo vertWireInfo;
      vertWireInfo.solid = vertWireBox;
      vertWireInfo.logical = vertWireLogical;
      vertWireInfo.name = wireGasName;
      finishNesting(vertWireInfo,
                    gasMaterial,
                    nullptr,
                    G4ThreeVector(0.0, gasY2, pwc->vertWireZ()),
                    container.logical,
                    wireNum,
                    visible, // visible
                    G4Colour::Blue(),
                    forceSolid, // forceSolid
                    forceAuxEdgeVisible, // forceAuxEdgeVisible
                    placePV,
                    doSurfaceCheck,
                    verbosity>0);
    }

  } // insertVerticalProfileWires

  void insertHorizontalProfileWires(VolumeInfo const& container,
                                    const PTMPWC* pwc,
                                    G4Material* gasMaterial,
                                    std::string const& wireNameSuffix,
                                    SimpleConfig const& _config,
                                    bool const visible,
                                    bool const forceSolid,
                                    bool const forceAuxEdgeVisible,
                                    bool const placePV,
                                    bool const doSurfaceCheck,
                                    int const verbosity) {
    // In the real detector wires run VERTICALLY so as to measure the
    // HORIZONTAL profile.
    // G4 geometry currently contains no actual wires.
    // Contains gas broken into sections -- each section represents the region
    // of the gas that is closest to one wire in the real detector.
    std::string wireGasNameHoriz = "PTMWireHoriz";
    wireGasNameHoriz.append(wireNameSuffix);
    G4VSolid *horizWireBox = new G4Box(wireGasNameHoriz,
                                     pwc->horizWireGasSection()->getXhalfLength(),
                                     pwc->horizWireGasSection()->getYhalfLength(),
                                     pwc->horizWireGasSection()->getZhalfLength());
    G4LogicalVolume *horizWireLogical = new G4LogicalVolume(horizWireBox, gasMaterial, wireGasNameHoriz);
    for (int i=0; i < pwc->numVertWires(); ++i) {
      int wireNum = pwc->wireNumStart() + pwc->numVertWires() + i;
      std::string wireGasName = wireGasNameHoriz;
      wireGasName.append(std::to_string(wireNum));
      // wire numbering such that the lowest-numered wire is
      // on the left, from the point of view of the oncoming beam.
      double gasX3 = pwc->horizWireXPos()[i];
      VolumeInfo horizWireInfo;
      horizWireInfo.solid = horizWireBox;
      horizWireInfo.logical = horizWireLogical;
      horizWireInfo.name = wireGasName;
      finishNesting(horizWireInfo,
                    gasMaterial,
                    nullptr,
                    G4ThreeVector(gasX3, 0.0, pwc->horizWireZ()),
                    container.logical,
                    wireNum,
                    visible, // visible
                    G4Colour::Blue(),
                    forceSolid, // forceSolid
                    forceAuxEdgeVisible, // forceAuxEdgeVisible
                    placePV,
                    doSurfaceCheck,
                    verbosity>0);
    }
  } // insertHorizontalProfileWires

  void constructTargetHallPWC(VolumeInfo const& motherVolume,
                              const PTMPWC* pwc,
                              int const vdNum,
                              SimpleConfig const& _config) {

    // collect geomOptions to be used in helper methods
    const int verbosity = _config.getInt("PTM.verbosityLevel",1);
    const auto& geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "PTM", "PTM" );
    const bool visible = geomOptions->isVisible("PTM");
    const bool forceSolid = geomOptions->isSolid("PTM");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("PTM");
    const bool doSurfaceCheck = geomOptions->doSurfaceCheck("PTM");
    bool placePV = geomOptions->placePV("PTM");
    // "container": box representing the location of the individual PWC
    G4Material* containerMaterial = motherVolume.logical->GetMaterial();

    std::vector<double> containerHalfDims;
    containerHalfDims.push_back(pwc->totalWidth()/2.);
    containerHalfDims.push_back(pwc->totalHeight()/2.);
    containerHalfDims.push_back(pwc->totalThick()/2.);

    std::string containerName = "PTMInnerContainer";
    containerName.append(pwc->nameSuffix());
    VolumeInfo pwcContainerInfo = nestBox(containerName,
              containerHalfDims,
              containerMaterial,
              nullptr,
              pwc->originInParent(),
              motherVolume,
              0,
              G4Colour::Green(),
              "PTM");

    // G10 frame of the PWC, made up of two endplates, plus
    // 13 interior boards that hold the windows and wires.
    // Represented as a single solid piece here.
    insertOuterFrame(pwcContainerInfo, pwc);

    insertWindows(pwcContainerInfo, pwc, vdNum, _config, visible, forceSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck, verbosity);

    // the sections of gas between the outer grounded planes and the bias planes
    // Going to use gasMaterial a few times; collect it out here so we don't do
    // a findMaterialOrThrow several times for the same material.
    G4Material *gasMaterial = findMaterialOrThrow(pwc->gasMaterialName());
    insertOuterGasBlocks(pwcContainerInfo, pwc, gasMaterial);

    // the wire planes
    // Going to use the wireNameSuffix in a couple of places, so just do the
    // string append once out here.
    std::string wireNameSuffix = pwc->nameSuffix();
    wireNameSuffix.append("_");
    insertVerticalProfileWires(pwcContainerInfo, pwc, gasMaterial, wireNameSuffix, _config, visible, forceSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck, verbosity);
    insertHorizontalProfileWires(pwcContainerInfo, pwc, gasMaterial, wireNameSuffix, _config, visible, forceSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck, verbosity);
  } // constructTargetHallPWC

  void constructPWCHolder(VolumeInfo const& motherVolume,
                          const PTMHead* ptmHead,
                          SimpleConfig const& _config ) {
    // collect geomOptions
    const int verbosity = _config.getInt("PTM.verbosityLevel",1);
    const auto& geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "PTM", "PTM" );
    const bool visible = geomOptions->isVisible("PTM");
    const bool forceSolid = geomOptions->isSolid("PTM");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("PTM");
    const bool doSurfaceCheck = geomOptions->doSurfaceCheck("PTM");
    bool placePV = geomOptions->placePV("PTM");

    double centerY = ptmHead->nearPWC()->originInParent().y();

    // pull out the pieces you need for the aluminum frame
    const Box* longExtrusion = ptmHead->holderExtrusionLong();
    const Box* shortExtrusion = ptmHead->holderExtrusionShort();
    std::string holderExtrusionMaterialName = ptmHead->holderExtrusionMaterialName();
    double holderExtrusionLongSep = ptmHead->holderExtrusionLongSep();
    double holderExtrusionShortPos = ptmHead->holderExtrusionShortPos();

    std::string longBarName = "PTMPWCHolderLong";
    G4Material *holderMaterial = findMaterialOrThrow(holderExtrusionMaterialName);
    G4VSolid *longBar = new G4Box("PTMPWCHolderLong",
                            longExtrusion->getXhalfLength(),
                            longExtrusion->getYhalfLength(),
                            longExtrusion->getZhalfLength());
    G4LogicalVolume *longBarLogical = new G4LogicalVolume(longBar, holderMaterial, longBarName);

    // 4 long bars holding the PWC at the right distance apart
    double bar1X = 0.5*holderExtrusionLongSep;
    double bar1Y = 0.5*holderExtrusionLongSep + centerY;
    VolumeInfo longBarInfo1;
    longBarInfo1.solid = longBar;
    longBarInfo1.logical = longBarLogical;
    longBarInfo1.name = longBarName+"01";
    finishNesting(longBarInfo1,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(bar1X, bar1Y, 0),
                  motherVolume.logical,
                  1,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    double bar2X = 0.5*holderExtrusionLongSep;
    double bar2Y = -0.5*holderExtrusionLongSep + centerY;
    VolumeInfo longBarInfo2;
    longBarInfo2.solid = longBar;
    longBarInfo2.logical = longBarLogical;
    longBarInfo2.name = longBarName+"02";
    finishNesting(longBarInfo2,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(bar2X, bar2Y, 0),
                  motherVolume.logical,
                  2,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    double bar3X = -0.5*holderExtrusionLongSep;
    double bar3Y = 0.5*holderExtrusionLongSep + centerY;
    VolumeInfo longBarInfo3;
    longBarInfo3.solid = longBar;
    longBarInfo3.logical = longBarLogical;
    longBarInfo3.name = longBarName+"03";
    finishNesting(longBarInfo3,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(bar3X, bar3Y, 0),
                  motherVolume.logical,
                  3,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    double bar4X = -0.5*holderExtrusionLongSep;
    double bar4Y = -0.5*holderExtrusionLongSep + centerY;
    VolumeInfo longBarInfo4;
    longBarInfo4.solid = longBar;
    longBarInfo4.logical = longBarLogical;
    longBarInfo4.name = longBarName+"04";
    finishNesting(longBarInfo4,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  nullptr, // no rotation
                  G4ThreeVector(bar4X, bar4Y, 0),
                  motherVolume.logical,
                  4,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);

    // four short bars, holding the long bars in the right positions relative to each other
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    std::string shortBarName = "PTMPWCHolderShort";
    G4VSolid *shortBar = new G4Box("PTMPWCHolderShort",
                            shortExtrusion->getXhalfLength(),
                            shortExtrusion->getYhalfLength(),
                            shortExtrusion->getZhalfLength());
    G4LogicalVolume *shortBarLogical = new G4LogicalVolume(shortBar, holderMaterial, shortBarName);
    double shortBarZ = holderExtrusionShortPos;

    G4RotationMatrix* vertBarRotation;
    vertBarRotation = reg.add(new G4RotationMatrix());
    vertBarRotation->rotateX(90.*CLHEP::deg);
    double short1X = 0.5*holderExtrusionLongSep;
    double short1Y = 0.0 + centerY;
    VolumeInfo shortBarInfo1;
    shortBarInfo1.solid = shortBar;
    shortBarInfo1.logical = shortBarLogical;
    shortBarInfo1.name = shortBarName+"01";
    finishNesting(shortBarInfo1,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  vertBarRotation,
                  G4ThreeVector(short1X, short1Y, shortBarZ),
                  motherVolume.logical,
                  1,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    double short2X = -0.5*holderExtrusionLongSep;
    double short2Y = 0.0 + centerY;
    VolumeInfo shortBarInfo2;
    shortBarInfo2.solid = shortBar;
    shortBarInfo2.logical = shortBarLogical;
    shortBarInfo2.name = shortBarName+"02";
    finishNesting(shortBarInfo2,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  vertBarRotation,
                  G4ThreeVector(short2X, short2Y, shortBarZ),
                  motherVolume.logical,
                  2,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);

    G4RotationMatrix* horizBarRotation;
    horizBarRotation = reg.add(new G4RotationMatrix());
    horizBarRotation->rotateY(90.*CLHEP::deg);
    double short3X = 0.0;
    double short3Y = 0.5*holderExtrusionLongSep + centerY;
    VolumeInfo shortBarInfo3;
    shortBarInfo3.solid = shortBar;
    shortBarInfo3.logical = shortBarLogical;
    shortBarInfo3.name = shortBarName+"03";
    finishNesting(shortBarInfo3,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  horizBarRotation,
                  G4ThreeVector(short3X, short3Y, shortBarZ),
                  motherVolume.logical,
                  3,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
    double short4X = 0.0;
    double short4Y = -0.5*holderExtrusionLongSep + centerY;
    VolumeInfo shortBarInfo4;
    shortBarInfo4.solid = shortBar;
    shortBarInfo4.logical = shortBarLogical;
    shortBarInfo4.name = shortBarName+"04";
    finishNesting(shortBarInfo4,
                  holderMaterial, // probably unnecessary, since we made the logical already
                  horizBarRotation,
                  G4ThreeVector(short4X, short4Y, shortBarZ),
                  motherVolume.logical,
                  4,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);

    // the handle at the top of the holder -- rectangle with top corners cut off, and an elliptical cutout
    double subtractionThickness = 1.1 * ptmHead->handleBase()->getZhalfLength();
    G4RotationMatrix* rotation45;
    rotation45 = reg.add(new G4RotationMatrix());
    rotation45->rotateZ(45.*CLHEP::deg);
    G4Box *handleOuter = new G4Box("ptmHandleOuter",
                ptmHead->handleBase()->getXhalfLength(),
                ptmHead->handleBase()->getYhalfLength(),
                ptmHead->handleBase()->getZhalfLength());
    G4EllipticalTube *handleHole = new G4EllipticalTube("ptmHandleHole",
                ptmHead->handleHoleSemiMinor(),
                ptmHead->handleHoleSemiMajor(),
                subtractionThickness);
    G4SubtractionSolid *handleLarge = new G4SubtractionSolid("ptmHandleLarge", handleOuter, handleHole, 0, ptmHead->handleHoleCenter());
    G4Box *cornerCutoff = new G4Box("ptmHandleCorner",
                0.5*ptmHead->handleCornerCutSide(),
                0.5*ptmHead->handleCornerCutSide(),
                subtractionThickness);
    G4SubtractionSolid *oneCorner = new G4SubtractionSolid("ptmHandleOneCorner",
                                           handleLarge,
                                           cornerCutoff,
                                           rotation45,
                                           G4ThreeVector(ptmHead->handleBase()->getXhalfLength(), ptmHead->handleBase()->getYhalfLength(), 0.0));
    G4SubtractionSolid *handleFinal = new G4SubtractionSolid("PTMHandle",
                                           oneCorner,
                                           cornerCutoff,
                                           rotation45,
                                           G4ThreeVector(-1*ptmHead->handleBase()->getXhalfLength(), ptmHead->handleBase()->getYhalfLength(), 0.0));
    std::string handleName = "PTMHandle";
    G4Material *handleMaterial = findMaterialOrThrow(ptmHead->handleMaterialName());
    G4ThreeVector handlePosition = G4ThreeVector(0.0, short3Y-shortExtrusion->getYhalfLength()+ptmHead->handleBase()->getYhalfLength(), shortBarZ+shortExtrusion->getXhalfLength()+ptmHead->handleBase()->getZhalfLength());


    VolumeInfo handleInfo;
    handleInfo.name = handleName;
    handleInfo.solid = handleFinal;
    finishNesting(handleInfo,
                  handleMaterial,
                  0,
                  handlePosition,
                  motherVolume.logical,
                  0,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);
  }

  void constructStand(VolumeInfo const& parent, const PTMStand* ptmStand, SimpleConfig const& _config) {
    // collect geomOptions
    const int verbosity = _config.getInt("PTM.verbosityLevel",1);
    const auto& geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "PTM", "PTM" );
    const bool visible = geomOptions->isVisible("PTM");
    const bool forceSolid = geomOptions->isSolid("PTM");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("PTM");
    const bool doSurfaceCheck = geomOptions->doSurfaceCheck("PTM");
    bool placePV = geomOptions->placePV("PTM");
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    // top wedge
    G4Material *wedgeMaterial = findMaterialOrThrow(ptmStand->wedgeMaterialName());
    G4ThreeVector wedgePositionInParent = ptmStand->topWedge()->getOffsetFromMu2eOrigin() - parent.centerInMu2e();
    G4ExtrudedSolid *outerWedge = new G4ExtrudedSolid("PTMWedgeOuter",
                     ptmStand->topWedge()->getVertices(),
                     ptmStand->topWedge()->getYhalfThickness(),
                     G4TwoVector(0.0, 0.0),
                     1,
                     G4TwoVector(0.0, 0.0),
                     1);
    //G4Box *cutout = new G4Box("PTMWedgeCutout", ptmStand->wedgeCutout()->getXhalfLength(), ptmStand->wedgeCutout()->getYhalfLength(), ptmStand->wedgeCutout()->getZhalfLength());
    G4EllipticalTube *cutout = new G4EllipticalTube("PTMWedgeCutout",
                ptmStand->wedgeCutoutSemiMinor(),
                ptmStand->wedgeCutoutSemiMajor(),
                ptmStand->topWedge()->getYhalfThickness()); // overkill, but will do the job
    G4RotationMatrix* wedgeRotation;
    wedgeRotation = reg.add(new G4RotationMatrix());
    wedgeRotation->rotateY(-90*CLHEP::deg); // so the long side points along z
    G4RotationMatrix* cutoutRotation;
    cutoutRotation = reg.add(new G4RotationMatrix());
    cutoutRotation->rotateX(90*CLHEP::deg); // so the long side points along z
    VolumeInfo wedgeInfo;
    wedgeInfo.name = "PTMStandTopWedge";
    wedgeInfo.solid = new G4SubtractionSolid("PTMStandTopWedge", outerWedge, cutout, cutoutRotation, G4ThreeVector(ptmStand->wedgeCutoutRelPosition()));
    finishNesting(wedgeInfo,
                  wedgeMaterial,
                  wedgeRotation,
                  wedgePositionInParent,
                  parent.logical,
                  0,
                  visible, // visible
                  G4Colour::Blue(),
                  forceSolid, // forceSolid
                  forceAuxEdgeVisible, // forceAuxEdgeVisible
                  placePV,
                  doSurfaceCheck,
                  verbosity>0);

    // support column
    std::string columnExtrusionName = "PTMBaseColumn";
    G4Material* columnMaterial = findMaterialOrThrow(ptmStand->columnMaterialName());
    G4Box* columnSolid = new G4Box(columnExtrusionName,
                                  ptmStand->columnExtrusion()->getXhalfLength(),
                                  ptmStand->columnExtrusion()->getYhalfLength(),
                                  ptmStand->columnExtrusion()->getZhalfLength());
    G4LogicalVolume* columnLogical = new G4LogicalVolume(columnSolid, columnMaterial, columnExtrusionName);
    int columnCopyNum = 0;
    for (auto columnPosition : ptmStand->columnOriginsInMu2e()) {
      G4RotationMatrix* columnRotation;
      columnRotation = reg.add(new G4RotationMatrix(ptmStand->columnRotations().at(columnCopyNum))); // TODO do this differently
      VolumeInfo columnInfo;
      columnInfo.name = columnExtrusionName + std::to_string(columnCopyNum);
      columnInfo.logical = columnLogical;
      columnInfo.solid = columnSolid;
      finishNesting(columnInfo,
                    columnMaterial,
                    columnRotation,
                    columnPosition - parent.centerInMu2e(),
                    parent.logical,
                    columnCopyNum,
                    visible, // visible
                    G4Colour::Blue(),
                    forceSolid, // forceSolid
                    forceAuxEdgeVisible, // forceAuxEdgeVisible
                    placePV,
                    doSurfaceCheck,
                    verbosity>0);
      columnCopyNum++;

    }
  }


  void constructPTM(VolumeInfo const& parent, SimpleConfig const& _config) {

    GeomHandle<PTM> ptmon;
    // collect geomOptions
    //const int verbosity = _config.getInt("PTM.verbosityLevel",1);
    const auto& geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "PTM", "PTM" );
    const bool visible = geomOptions->isVisible("PTM");
    const bool forceSolid = geomOptions->isSolid("PTM");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("PTM");
    const bool doSurfaceCheck = geomOptions->doSurfaceCheck("PTM");
    bool placePV = geomOptions->placePV("PTM");

    // create and place the main volume first
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    G4ThreeVector parentPosition = parent.centerInMu2e(); // HallAir
    G4RotationMatrix* parentRotation = parent.physical->GetObjectRotation();
    // The "main" volume that contains the stand, detectors, everything
    G4ThreeVector mainPTMPositionMu2e = ptmon->originInMu2e();
    G4ThreeVector mainPTMPositionInParent = mainPTMPositionMu2e - parentPosition;
    G4RotationMatrix* mainRotation;
    if (parentRotation->isIdentity()) {
      mainRotation = reg.add(new G4RotationMatrix(ptmon->rotationInMu2e()));
    } else {
      mainRotation = reg.add(new G4RotationMatrix(ptmon->rotationInMu2e()));
      mainRotation->transform(parentRotation->inverse());
    }
    G4Material* mainMaterial = parent.logical->GetMaterial();
    std::vector<double> mainHalfDims;
    mainHalfDims.push_back(ptmon->totalWidth()/2.);
    mainHalfDims.push_back(ptmon->totalHeight()/2.);
    mainHalfDims.push_back(ptmon->totalLength()/2.);
    VolumeInfo pTargetMonMain = nestBox("PTMMain",
               mainHalfDims,
               mainMaterial,
               mainRotation,
               mainPTMPositionInParent,
               parent,
               0,
               visible,
               G4Colour::Green(),
               forceSolid,
               forceAuxEdgeVisible,
               placePV,
               doSurfaceCheck
               );

    // the "mother" volume contains the PWCs, the PWC holder, and the handle for the RHS
    G4ThreeVector motherPosition = ptmon->ptmHead()->originInMu2e() - mainPTMPositionMu2e;
    // try making the rotation matrix object first, then doing the
    // main rotation, then un-rotating by whatever the parent's rotation is
    G4RotationMatrix* motherRotation;
    if (mainRotation->isIdentity()) {
      motherRotation = reg.add(new G4RotationMatrix(ptmon->ptmHead()->rotationInMu2e()));
    } else {
      motherRotation = reg.add(new G4RotationMatrix(ptmon->ptmHead()->rotationInMu2e()));
      *motherRotation = *motherRotation * ptmon->rotationInMu2e().inverse(); // transform() reverses the order of multiplication, gives the wrong result
    }

    std::vector<double> motherHalfDims;
    motherHalfDims.push_back(ptmon->ptmHead()->totalWidth()/2.);
    motherHalfDims.push_back(ptmon->ptmHead()->totalHeight()/2.);
    motherHalfDims.push_back(ptmon->ptmHead()->totalLength()/2.);

    VolumeInfo pTargetMonContainer = nestBox("PTMMother",
                motherHalfDims,
                mainMaterial, // same material as main volume
                motherRotation,
                motherPosition,
                pTargetMonMain,
                0,
                G4Colour::Green(),
                "PTM");

    // add the first PWC to the mother volume
    constructTargetHallPWC(pTargetMonContainer, ptmon->ptmHead()->nearPWC(), VirtualDetectorId::PTM_1_In, _config);
    // and the second PWC
    constructTargetHallPWC(pTargetMonContainer, ptmon->ptmHead()->farPWC(), VirtualDetectorId::PTM_2_In, _config);
    // the aluminum frame that holds the PWCs
    constructPWCHolder(pTargetMonContainer,
                          ptmon->ptmHead(),
                          _config);

    constructStand(pTargetMonMain, ptmon->ptmStand(), _config);

  } // constructPTM



} // namespace mu2e
