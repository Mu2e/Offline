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
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

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

  // nameSuffix lets you give multiple PWC's unique names.
  // wireNumStart is the first number to be used when naming the sections of
  // gas corresponding to wires. This same wire numbering is used for the
  // copyNo argument when placing the wire gas in the geometry.
  void constructTargetHallPWC(VolumeInfo const & parent, SimpleConfig const & _config, std::string const & nameSuffix, G4ThreeVector const & positionInParent, int const wireNumStart) {
    double gasLength = _config.getDouble("pTargetMon_gasLength");
    double outerPlateLength = _config.getDouble("pTargetMon_outerPlateLength");
    double windowWidth = _config.getDouble("pTargetMon_windowWidth");
    double windowHeight = _config.getDouble("pTargetMon_windowHeight");
    double height = _config.getDouble("pTargetMon_height");
    double width = _config.getDouble("pTargetMon_width");
    double windowThick = _config.getDouble("pTargetMon_windowThick");
    double frameThick = _config.getDouble("pTargetMon_frameThick");
    int wiresPerPlane = _config.getInt("pTargetMon_wiresPerPlane");
    double detectorLength = gasLength + (2*outerPlateLength) + frameThick;
    std::string wireNameSuffix = nameSuffix;
    wireNameSuffix.append("_");

    std::vector<double> halfDims;
    halfDims.push_back(width/2.);
    halfDims.push_back(height/2.);
    halfDims.push_back(detectorLength/2.);

    G4Material* baseMaterial = parent.logical->GetMaterial();

    // "container": box representing the location of the individual PWC

    std::string containerName = "pTargetMonInnerContainer";
    containerName.append(nameSuffix);
    VolumeInfo PWCContainerInfo = nestBox(containerName,
              halfDims,
              baseMaterial,
              nullptr,
              positionInParent,
              parent,
              0,
              G4Colour::Green(),
              "PTM");

    
    // G10 frame of the PWC, represented as one piece here
    G4Box *outerBox = new G4Box("pwcFrameOuter", 
                height/2., 
                width/2., 
                detectorLength/2.);
    G4Box *innerBox = new G4Box("pwcFrameInner", 
                0.001+windowHeight/2., 
                0.001+windowWidth/2., 
                0.001+detectorLength/2.);
    std::string frameName = "pTargetMonFrame";
    frameName.append(nameSuffix);
    G4Material *frameMaterial = findMaterialOrThrow(_config.getString("pTargetMon_frameMaterial"));

    VolumeInfo frameInfo;
    frameInfo.name = frameName;
    frameInfo.solid = new G4SubtractionSolid(frameName, outerBox, innerBox);
    finishNesting(frameInfo, 
          frameMaterial, 
          nullptr, 
          G4ThreeVector(0.0, 0.0, 0.0), 
          PWCContainerInfo.logical, 
          0, 
          G4Colour::Blue(), 
          "PTM");

    // insert the windows
    G4Material *windowMaterial = findMaterialOrThrow(_config.getString("pTargetMon_windowMaterial"));
    std::vector<double> windowHalfDims;
    windowHalfDims.push_back(windowWidth/2.);
    windowHalfDims.push_back(windowHeight/2.);
    windowHalfDims.push_back(windowThick/2.);
    G4VSolid* windowBox = new G4Box("pTargetMonWindow",
                  windowWidth/2.,
                  windowHeight/2.,
                  windowThick/2.);
    G4LogicalVolume* windowLogical = new G4LogicalVolume(windowBox,
                              windowMaterial,
                              "pTargetMonWindow");
    // first ground plane
    std::string ground1Name = "pTargetMonGroundIn";
    ground1Name.append(nameSuffix);
    double ground1Z = -5.5*frameThick;
    //G4VPhysicalVolume* pv = 
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, ground1Z),
              windowLogical,
              ground1Name,
              PWCContainerInfo.logical,
              false,
              0,
              false);
    // first HV plane
    std::string hv1Name = "pTargetMonHV1";
    hv1Name.append(nameSuffix);
    double hv1Z = -3.5*frameThick;
    //G4VPhysicalVolume* pv = 
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, hv1Z),
              windowLogical,
              hv1Name,
              PWCContainerInfo.logical,
              false,
              0,
              false);
    // second HV plane
    std::string hv2Name = "pTargetMonHV2";
    hv2Name.append(nameSuffix);
    double hv2Z = 0.5*frameThick;
    //G4VPhysicalVolume* pv = 
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, hv2Z),
              windowLogical,
              hv2Name,
              PWCContainerInfo.logical,
              false,
              0,
              false);
    // third HV plane
    std::string hv3Name = "pTargetMonHV3";
    hv3Name.append(nameSuffix);
    double hv3Z = 4.5*frameThick;
    //G4VPhysicalVolume* pv = 
    new G4PVPlacement(nullptr
              G4ThreeVector(0.0, 0.0, hv3Z),
              windowLogical,
              hv3Name,
              PWCContainerInfo.logical,
              false,
              0,
              false);
    // last ground plane
    std::string ground2Name = "pTargetMonGroundOut";
    ground2Name.append(nameSuffix);
    double ground2Z = 6.5*frameThick;
    new G4PVPlacement(nullptr,
              G4ThreeVector(0.0, 0.0, ground2Z),
              windowLogical,
              ground2Name,
              PWCContainerInfo.logical,
              false,
              0,
              false);
    
    // gas inside PWC
    G4Material *gasMaterial = findMaterialOrThrow(_config.getString("pTargetMon_innerGas"));
    // between ground plane 1 and HV plane 1
    std::string gasName1 = "pTargetMonGas1";
    gasName1.append(nameSuffix);
    double gasLength1 = hv1Z - ground1Z - windowThick;
    double gasZ1 = 0.5*(hv1Z + ground1Z);
    std::vector<double> gasHalfDims1;
    gasHalfDims1.push_back(windowWidth/2.);
    gasHalfDims1.push_back(windowHeight/2.);
    gasHalfDims1.push_back(gasLength1/2.);
    nestBox(gasName1,
        gasHalfDims1,
        gasMaterial,
        nullptr,
        G4ThreeVector(0.0, 0.0, gasZ1),
        PWCContainerInfo,
        0,
        G4Colour::Red(),
        "PTM");
    // between HV plane 1 and HV plane 2
    // "Vert" here means it measures the vertical profile
    std::string wireGasNameVert = "pTargetMonWireVert";
    wireGasNameVert.append(wireNameSuffix);
    double gasLength2 = hv2Z - hv1Z - windowThick;
    double gasZ2 = 0.5*(hv2Z + hv1Z);
    double wireSpacing = windowHeight/wiresPerPlane;
    G4VSolid* vertWireBox = new G4Box(wireGasNameVert,
                      windowWidth/2.,
                      wireSpacing/2.,
                      gasLength2/2.);
    G4LogicalVolume* vertWireLogical = new G4LogicalVolume(vertWireBox,
                                gasMaterial,
                                wireGasNameVert);
    for (int i=0; i<wiresPerPlane; i++) {
      int wireNum = wireNumStart + i;
      std::string wireGasName = wireGasNameVert;
      wireGasName.append(std::to_string(wireNum));
      // wire numbering such that the lowest-numered wire is 
      // on the bottom
      double gasY2 = (-0.5*windowHeight) + ((i+0.5)*wireSpacing);
      //G4VPhysicalVolume* pv = 
      new G4PVPlacement(nullptr,
                G4ThreeVector(0.0, gasY2, gasZ2),
                vertWireLogical,
                wireGasName,
                PWCContainerInfo.logical,
                false,
                wireNum,
                false);
    }
    
    // "Horiz" here means it measures the horizontal profile
    std::string wireGasNameHoriz = "pTargetMonWireHoriz";
    wireGasNameHoriz.append(wireNameSuffix);
    double gasLength3 = hv3Z - hv2Z - windowThick;
    double gasZ3 = 0.5*(hv3Z + hv2Z);
    wireSpacing = windowWidth/wiresPerPlane;
    G4VSolid* horizWireBox = new G4Box(wireGasNameHoriz,
                       wireSpacing/2.,
                       windowHeight/2.,
                       gasLength3/2.);
    G4LogicalVolume* horizWireLogical = new G4LogicalVolume(horizWireBox,
                                gasMaterial,
                                wireGasNameHoriz);
    int horizWireStart = wireNumStart+wiresPerPlane;
    for (int i=0; i<wiresPerPlane; i++) {
      int wireNum = horizWireStart + i;
      std::string wireGasName = wireGasNameHoriz;
      wireGasName.append(std::to_string(wireNum));
      // Wire numbering is reversed, because we rotate the overall 
      // container by almost 180 degrees to be along the beam path.
      // This puts the lowest-numbered wire all the way on the left
      // from the oncoming beam's point of view.
      double gasX3 = (0.5*windowWidth) - ((i+0.5)*wireSpacing);
      //G4VPhysicalVolume* pv = 
      new G4PVPlacement(nullptr,
                G4ThreeVector(gasX3, 0.0, gasZ3),
                horizWireLogical,
                wireGasName,
                PWCContainerInfo.logical,
                false,
                wireNum,
                false);
    }

    // between HV plane 3 and ground plane 2
    std::string gasName4 = "pTargetMonGas4";
    gasName4.append(nameSuffix);
    double gasLength4 = ground2Z - hv3Z - windowThick;
    double gasZ4 = 0.5*(ground2Z + hv3Z);
    std::vector<double> gasHalfDims4;
    gasHalfDims4.push_back(windowWidth/2.);
    gasHalfDims4.push_back(windowHeight/2.);
    gasHalfDims4.push_back(gasLength4/2.);
    nestBox(gasName4,
        gasHalfDims4,
        gasMaterial,
        nullptr,
        G4ThreeVector(0.0, 0.0, gasZ4),
        PWCContainerInfo,
        0,
        G4Colour::Red(),
        "PTM");

    
  } //constructTargetHallPWC


  void constructProductionTargetMon(VolumeInfo const & parent, SimpleConfig const & _config) {
    int wiresPerPlane = _config.getInt("pTargetMon_wiresPerPlane");
    double height = _config.getDouble("pTargetMon_height");
    double width = _config.getDouble("pTargetMon_width");
    double length = _config.getDouble("pTargetMon_length");
    std::vector<double> halfDims;
    halfDims.push_back(0.1+width/2.);
    halfDims.push_back(0.1+height/2.);
    halfDims.push_back(0.1+length/2.);

    G4double xPosInMu2e = _config.getDouble("pTargetMon_positionX");
    G4double yPosInMu2e = _config.getDouble("pTargetMon_positionY");
    G4double zPosInMu2e = _config.getDouble("pTargetMon_positionZ");
    G4ThreeVector positionInMu2e = G4ThreeVector(xPosInMu2e, yPosInMu2e, zPosInMu2e);

    double yRotInMu2e = _config.getDouble("pTargetMon_rotY");
    double xRotInMu2e = _config.getDouble("pTargetMon_rotX");
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    G4RotationMatrix* rotation = reg.add(new G4RotationMatrix);
    //G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateY(yRotInMu2e*CLHEP::deg);
    rotation->rotateX(xRotInMu2e*CLHEP::deg);

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    G4Material* baseMaterial = parent.logical->GetMaterial();

    // container: holds the 2 actual detectors
    VolumeInfo pTargetMonContainer = nestBox("pTargetMonMother",
                halfDims,
                baseMaterial,
                rotation,
                positionInMu2e-_hallOriginInMu2e,
                parent,
                0,
                G4Colour::Green(),
                "PTM");
    checkForOverlaps( pTargetMonContainer.physical, _config, true);

    double gasLength = _config.getDouble("pTargetMon_gasLength");
    double outerPlateLength = _config.getDouble("pTargetMon_outerPlateLength");
    double detectorLength = gasLength + (2*outerPlateLength);

    // location of one PWC:
    double z1 = (-0.5*length)+(detectorLength/2.);
    G4ThreeVector position_1 = G4ThreeVector(0.0, 0.0, z1);
    constructTargetHallPWC(pTargetMonContainer, _config, "_1", position_1, 0);
    // second PWC:
    double z2 = (0.5*length)-(detectorLength/2.);
    G4ThreeVector position_2 = G4ThreeVector(0.0, 0.0, z2);
    constructTargetHallPWC(pTargetMonContainer, _config, "_2", position_2, 2*wiresPerPlane);

  } // constructProductionTargetMon


  
} // namespace mu2e