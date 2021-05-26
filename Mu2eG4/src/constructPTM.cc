//
// Free function. Approach borrowed from constructPS
// Constructs the downstream production target scanning monitor.
// Parent volume is the air in the target hall. Probably?
//

// C++ includes
#include <iostream>
#include <vector>
#include <string>

// Mu2e includes
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "PTMGeom/inc/PTM.hh"
#include "PTMGeom/inc/PTMPWC.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

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
                        const PTMPWC* pwc, 
                        SimpleConfig const& _config,
                        bool const doSurfaceCheck,
                        int const verbosity) {
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

    if (doSurfaceCheck) checkForOverlaps( frameInfo.physical, _config, verbosity>0);

  } // insertOuterFrame

  void insertWindows(VolumeInfo const& container, 
                     const PTMPWC* pwc, 
                     int const vdNum,
                     SimpleConfig const& _config,
                     bool const doSurfaceCheck,
                     int const verbosity) {
    // The windows are all identical, so make one logical volume and paste it repeatedly.
    // The one exception is the most-upstream window, which is a virtual detector.
    G4Material *windowMaterial = findMaterialOrThrow(pwc->windowMaterialName());
    G4VSolid* windowBox = new G4Box("PTMWindow",
                  pwc->pwcWindow()->getXhalfLength(),
                  pwc->pwcWindow()->getYhalfLength(),
                  pwc->pwcWindow()->getZhalfLength());
    G4LogicalVolume* windowLogical = new G4LogicalVolume(windowBox,
                              windowMaterial,
                              "PTMWindow");
    G4LogicalVolume* windowLogicalVD = new G4LogicalVolume(windowBox,
                              windowMaterial,
                              "VirtualDetector_PTMGroundIn"+pwc->nameSuffix());
    std::string ground1Name = "VirtualDetector_PTMGroundIn"+pwc->nameSuffix();
    G4VPhysicalVolume* ground1Phys =
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, pwc->ground1Z()),
              windowLogicalVD,
              ground1Name,
              container.logical,
              false,
              vdNum,
              false);
    std::string hv1Name = "PTMHV1"+pwc->nameSuffix();
    G4VPhysicalVolume* hv1Phys =
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, pwc->hv1Z()),
              windowLogical,
              ground1Name,
              container.logical,
              false,
              0,
              false);
    std::string hv2Name = "PTMHV2"+pwc->nameSuffix();
    G4VPhysicalVolume* hv2Phys =
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, pwc->hv2Z()),
              windowLogical,
              ground1Name,
              container.logical,
              false,
              0,
              false);
    std::string hv3Name = "PTMHV3"+pwc->nameSuffix();
    G4VPhysicalVolume* hv3Phys =
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, pwc->hv3Z()),
              windowLogical,
              ground1Name,
              container.logical,
              false,
              0,
              false);
    std::string ground2Name = "PTMGroundOut"+pwc->nameSuffix();
    G4VPhysicalVolume* ground2Phys =
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, pwc->ground2Z()),
              windowLogical,
              ground2Name,
              container.logical,
              false,
              0,
              false);

    if (doSurfaceCheck){
      checkForOverlaps( ground1Phys, _config, verbosity>0);
      checkForOverlaps( hv1Phys, _config, verbosity>0);
      checkForOverlaps( hv2Phys, _config, verbosity>0);
      checkForOverlaps( hv3Phys, _config, verbosity>0);
      checkForOverlaps( ground2Phys, _config, verbosity>0);
    }
    

  } // insertWindows

  void insertOuterGasBlocks(VolumeInfo const& container, 
                            const PTMPWC* pwc, 
                            G4Material* gasMaterial, 
                            SimpleConfig const& _config,
                            bool const doSurfaceCheck,
                            int const verbosity) {
    // between ground plane 1 and HV plane 1
    std::string gasName1 = "PTMGas1";
    gasName1.append(pwc->nameSuffix());
    std::vector<double> gas1HalfDims;
    gas1HalfDims.push_back(pwc->gasSection1()->getXhalfLength());
    gas1HalfDims.push_back(pwc->gasSection1()->getYhalfLength());
    gas1HalfDims.push_back(pwc->gasSection1()->getZhalfLength());
    VolumeInfo gas1Info =
    nestBox (gasName1,
             gas1HalfDims,
             gasMaterial,
             nullptr,
             G4ThreeVector(0.0, 0.0, pwc->gasInZ()),
             container,
             0, // copyNo
             false, // isVisible
             G4Colour::Red(),
             true, // forceSolid
             false, // forceAuxEdgeVisible
             true, // placePV
             false); // doSurfaceCheck


    // between HV plane 3 and ground plane 2
    std::string gasName4 = "PTMGas4";
    gasName4.append(pwc->nameSuffix());
    std::vector<double> gas4HalfDims;
    gas4HalfDims.push_back(pwc->gasSection4()->getXhalfLength());
    gas4HalfDims.push_back(pwc->gasSection4()->getYhalfLength());
    gas4HalfDims.push_back(pwc->gasSection4()->getZhalfLength());
    VolumeInfo gas4Info =
    nestBox (gasName4,
             gas4HalfDims,
             gasMaterial,
             nullptr,
             G4ThreeVector(0.0, 0.0, pwc->gasOutZ()),
             container,
             0, // copyNo
             false, // isVisible
             G4Colour::Red(),
             true, // forceSolid
             false, // forceAuxEdgeVisible
             true, // placePV
             false); // doSurfaceCheck

    if (doSurfaceCheck) {
      checkForOverlaps( gas1Info.physical, _config, verbosity>0);
      checkForOverlaps( gas4Info.physical, _config, verbosity>0);
    }

  } // insertOuterGasBlocks

  void insertVerticalProfileWires(VolumeInfo const& container, 
                                  const PTMPWC* pwc, 
                                  G4Material* gasMaterial, 
                                  std::string const& wireNameSuffix, 
                                  SimpleConfig const& _config,
                                  bool const doSurfaceCheck,
                                  int const verbosity) {
    // In the real detector wires run HORIZONTALLY, so as to measure the 
    // VERTICAL profile.
    // G4 geometry currently contains no actual wires.
    // Contains gas broken into sections -- each section represents the region 
    // of the gas that is closest to one wire in the real detector.
    std::string wireGasNameVert = "PTMWireVert";
    wireGasNameVert.append(wireNameSuffix);
    G4VSolid* vertWireBox = new G4Box(wireGasNameVert,
                      pwc->vertWireGasSection()->getXhalfLength(),
                      pwc->vertWireGasSection()->getYhalfLength(),
                      pwc->vertWireGasSection()->getZhalfLength());
    G4LogicalVolume* vertWireLogical = new G4LogicalVolume(vertWireBox,
                                gasMaterial,
                                wireGasNameVert);
    for (int i=0; i < pwc->numVertWires(); ++i) {
      int wireNum = pwc->wireNumStart() + i;
      std::string wireGasName = wireGasNameVert;
      wireGasName.append(std::to_string(wireNum));
      // wire numbering such that the lowest-numered wire is 
      // on the bottom
      double gasY2 = pwc->vertWireYPos()[i];
      G4VPhysicalVolume* wireGas =
      new G4PVPlacement(nullptr,
                G4ThreeVector(0.0, gasY2, pwc->vertWireZ()),
                vertWireLogical,
                wireGasName,
                container.logical,
                false,
                wireNum,
                false);
      if (doSurfaceCheck) checkForOverlaps( wireGas, _config, verbosity>0);
    }

  } // insertVerticalProfileWires

  void insertHorizontalProfileWires(VolumeInfo const& container, 
                                    const PTMPWC* pwc, 
                                    G4Material* gasMaterial, 
                                    std::string const& wireNameSuffix, 
                                    SimpleConfig const& _config,
                                    bool const doSurfaceCheck,
                                    int const verbosity) {
    // In the real detector wires run VERTICALLY so as to measure the 
    // HORIZONTAL profile.
    // G4 geometry currently contains no actual wires.
    // Contains gas broken into sections -- each section represents the region 
    // of the gas that is closest to one wire in the real detector.
    std::string wireGasNameHoriz = "PTMWireHoriz";
    wireGasNameHoriz.append(wireNameSuffix);
    G4VSolid* horizWireBox = new G4Box(wireGasNameHoriz,
                      pwc->horizWireGasSection()->getXhalfLength(),
                      pwc->horizWireGasSection()->getYhalfLength(),
                      pwc->horizWireGasSection()->getZhalfLength());
    G4LogicalVolume* horizWireLogical = new G4LogicalVolume(horizWireBox,
                                gasMaterial,
                                wireGasNameHoriz);
    for (int i=0; i < pwc->numHorizWires(); ++i) {
      int wireNum = pwc->wireNumStart() + pwc->numVertWires() + i;
      std::string wireGasName = wireGasNameHoriz;
      wireGasName.append(std::to_string(wireNum));
      // wire numbering such that the lowest-numered wire is 
      // on the bottom
      double gasX3 = pwc->horizWireXPos()[i];
      G4VPhysicalVolume* wireGas =
      new G4PVPlacement(nullptr,
                G4ThreeVector(gasX3, 0.0, pwc->horizWireZ()),
                horizWireLogical,
                wireGasName,
                container.logical,
                false,
                wireNum,
                false);
      if (doSurfaceCheck) checkForOverlaps( wireGas, _config, verbosity>0);
    }

  } // insertHorizontalProfileWires

  void constructTargetHallPWC(VolumeInfo const& motherVolume, 
                              const PTMPWC* pwc, 
                              int const vdNum,
                              SimpleConfig const& _config, 
                              bool const doSurfaceCheck, 
                              int const verbosity) {
    
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
    if (doSurfaceCheck) checkForOverlaps( pwcContainerInfo.physical, _config, verbosity>0);

    // G10 frame of the PWC, made up of two endplates, plus
    // 13 interior boards that hold the windows and wires.
    // Represented as a single solid piece here.
    insertOuterFrame(pwcContainerInfo, pwc, _config, doSurfaceCheck, verbosity);

    insertWindows(pwcContainerInfo, pwc, vdNum, _config, doSurfaceCheck, verbosity);

    // the sections of gas between the outer grounded planes and the bias planes
    // Going to use gasMaterial a few times; collect it out here so we don't do
    // a findMaterialOrThrow several times for the same material.
    G4Material *gasMaterial = findMaterialOrThrow(pwc->gasMaterialName());
    insertOuterGasBlocks(pwcContainerInfo, pwc, gasMaterial, _config, doSurfaceCheck, verbosity);

    // the wire planes
    // Going to use the wireNameSuffix in a couple of places, so just do the 
    // string append once out here.
    std::string wireNameSuffix = pwc->nameSuffix();
    wireNameSuffix.append("_");
    insertVerticalProfileWires(pwcContainerInfo, pwc, gasMaterial, wireNameSuffix, _config, doSurfaceCheck, verbosity);
    insertHorizontalProfileWires(pwcContainerInfo, pwc, gasMaterial, wireNameSuffix, _config, doSurfaceCheck, verbosity);
  } // constructTargetHallPWC


  void constructPTM(VolumeInfo const& parent, SimpleConfig const& _config) {
    const int verbosity = _config.getInt("PTM_verbosity",1);
    const auto& geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "PTM", "PTM" );
    const bool doSurfaceCheck = geomOptions->doSurfaceCheck("PTM");

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
      motherRotation = reg.add(new G4RotationMatrix(parentRotation->inverse()));
      motherRotation->transform(ptmon->rotationInMu2e());
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
    if (doSurfaceCheck) checkForOverlaps( pTargetMonContainer.physical, _config, verbosity>0);

    // add the first PWC to the mother volume
    constructTargetHallPWC(pTargetMonContainer, ptmon->nearPWC(), VirtualDetectorId::PTM_1_In, _config, doSurfaceCheck, verbosity);
    // and the second PWC
    constructTargetHallPWC(pTargetMonContainer, ptmon->farPWC(), VirtualDetectorId::PTM_2_In, _config, doSurfaceCheck, verbosity);

  } // constructPTM


  
} // namespace mu2e