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
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 inludes
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4Colour.hh"
#include "G4SubtractionSolid.hh"

#include "CLHEP/Units/SystemOfUnits.h"


using namespace std;

namespace mu2e {

    // TODO: add more structure
    // TODO: run mu2e -c Mu2eG4/fcl/gdmldump.fcl

    void constructTargetHallPWC(VolumeInfo const & parent, SimpleConfig const & _config, std::string nameSuffix, G4ThreeVector positionInParent) {
        double gasLength = _config.getDouble("pTargetMon_gasLength");
        double outerPlateLength = _config.getDouble("pTargetMon_outerPlateLength");
        double detectorLength = gasLength + (2*outerPlateLength);
        double windowWidth = _config.getDouble("pTargetMon_windowWidth");
        double windowHeight = _config.getDouble("pTargetMon_windowHeight");
        double height = _config.getDouble("pTargetMon_height");
        double width = _config.getDouble("pTargetMon_width");

        std::vector<double> halfDims;
        halfDims.push_back(width/2.);
        halfDims.push_back(height/2.);
        halfDims.push_back(detectorLength/2.);

        AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
        G4RotationMatrix* noRotation = reg.add(new G4RotationMatrix);

        G4Material* baseMaterial = parent.logical->GetMaterial();

        // "container": box representing the location of the individual PWC

        std::string containerName = "PWCContainer";
        containerName.append(nameSuffix);
        VolumeInfo PWCContainerInfo = nestBox(containerName,
                            halfDims,
                            baseMaterial,
                            noRotation,
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
                    noRotation, 
                    G4ThreeVector(0.0, 0.0, 0.0), 
                    PWCContainerInfo.logical, 
                    0, 
                    G4Colour::Blue(), 
                    "PTM");
        
        // gas inside PWC
        G4Material *gasMaterial = findMaterialOrThrow(_config.getString("pTargetMon_innerGas"));
        std::string gasName = "pTargetMonGas";
        gasName.append(nameSuffix);
        std::vector<double> gasHalfDims;
        gasHalfDims.push_back(windowWidth/2.);
        gasHalfDims.push_back(windowHeight/2.);
        gasHalfDims.push_back(gasLength/2.);
        nestBox(gasName,
                gasHalfDims,
                gasMaterial,
                noRotation,
                G4ThreeVector(0.0, 0.0, 0.0),
                PWCContainerInfo,
                0,
                G4Colour::Red(),
                "PTM");
        
    } //constructTargetHallPWC


    void constructProductionTargetMon(VolumeInfo const & parent, SimpleConfig const & _config) {
        cout << endl << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Now ENTERING constructProductionTargetMon" << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
        // for now, do this in a simple way that hopefully works?
        double height = _config.getDouble("pTargetMon_height");
        double width = _config.getDouble("pTargetMon_width");
        double length = _config.getDouble("pTargetMon_length");
        //double halfDimsRaw[3] = {width/2., height/2., length/2.};
        std::vector<double> halfDims;
        halfDims.push_back(width/2.);
        halfDims.push_back(height/2.);
        halfDims.push_back(length/2.);

        G4double xPosInMu2e = _config.getDouble("pTargetMon_positionX");
        G4double yPosInMu2e = _config.getDouble("pTargetMon_positionY");
        G4double zPosInMu2e = _config.getDouble("pTargetMon_positionZ");
        G4ThreeVector positionInMu2e = G4ThreeVector(xPosInMu2e, yPosInMu2e, zPosInMu2e);

        double yRotInMu2e = _config.getDouble("pTargetMon_rotY")*-1.;
        double xRotInMu2e = _config.getDouble("pTargetMon_rotX")*-1.;
        AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
        G4RotationMatrix* rotation = reg.add(new G4RotationMatrix);
        //G4RotationMatrix* rotation = new G4RotationMatrix();
        rotation->rotateY(yRotInMu2e*CLHEP::deg);
        rotation->rotateX(xRotInMu2e*CLHEP::deg);

        G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

        G4Material* baseMaterial = parent.logical->GetMaterial();

        // container: holds the 2 actual detectors
        VolumeInfo pTargetMonContainer = nestBox("pTargetMonContainer",
                                halfDims,
                                baseMaterial,
                                rotation,
                                positionInMu2e-_hallOriginInMu2e,
                                parent,
                                0,
                                G4Colour::Green(),
                                "PTM");

        double gasLength = _config.getDouble("pTargetMon_gasLength");
        double outerPlateLength = _config.getDouble("pTargetMon_outerPlateLength");
        double detectorLength = gasLength + (2*outerPlateLength);

        // location of one PWC:
        double z1 = (-0.5*length)+(detectorLength/2.);
        G4ThreeVector position_1 = G4ThreeVector(0.0, 0.0, z1);
        constructTargetHallPWC(pTargetMonContainer, _config, "_1", position_1);
        // second PWC:
        double z2 = (0.5*length)-(detectorLength/2.);
        G4ThreeVector position_2 = G4ThreeVector(0.0, 0.0, z2);
        constructTargetHallPWC(pTargetMonContainer, _config, "_2", position_2);


        cout << endl << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Now LEAVING constructProductionTargetMon" << endl;
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;

    } // constructProductionTargetMon


    

    
} // namespace mu2e