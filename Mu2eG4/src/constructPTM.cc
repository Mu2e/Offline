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


  void constructPTM(VolumeInfo const& parent, SimpleConfig const& _config) {

    GeomHandle<PTM> ptmon;

    // create and place the mother volume first
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    G4ThreeVector parentPosition = parent.centerInMu2e();
    G4RotationMatrix* parentRotation = parent.physical->GetObjectRotation();
    G4ThreeVector motherPosition = ptmon->originInMu2e() - parentPosition;
    // try making the rotation matrix object first, then doing the
    // main rotation, then un-rotating by whatever the parent's rotation is
    G4RotationMatrix* motherRotation;
    if (parentRotation->isIdentity()) {
      motherRotation = reg.add(new G4RotationMatrix(ptmon->rotationInMu2e()));
    } else {
      motherRotation = reg.add(new G4RotationMatrix(ptmon->rotationInMu2e()));
      motherRotation->transform(parentRotation->inverse());
    }

    G4Material* motherMaterial = parent.logical->GetMaterial();
    std::vector<double> motherHalfDims;
    motherHalfDims.push_back(ptmon->totalWidth()/2.);
    motherHalfDims.push_back(ptmon->totalHeight()/2.);
    motherHalfDims.push_back(ptmon->totalLength()/2.);

    VolumeInfo pTargetMonContainer = nestBox("PTMMother",
                motherHalfDims,
                motherMaterial,
                motherRotation,
                motherPosition,
                parent,
                0,
                G4Colour::Green(),
                "PTM");

    // add the first PWC to the mother volume
    constructTargetHallPWC(pTargetMonContainer, ptmon->nearPWC(), VirtualDetectorId::PTM_1_In, _config);
    // and the second PWC
    constructTargetHallPWC(pTargetMonContainer, ptmon->farPWC(), VirtualDetectorId::PTM_2_In, _config);

  } // constructPTM



} // namespace mu2e
