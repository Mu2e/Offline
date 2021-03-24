//
// Free function to create MSTM.
// Muon Stopping Target Monitor
//
// constructMSTM.cc
// Author: A. Palladino
// Date: see git for version history
//
// Original free function author K.L.Genser 
//

// Mu2e includes.
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructMSTM.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4UniformMagField.hh"
#include "Geant4/G4Mag_UsualEqRhs.hh"
#include "Geant4/G4ExactHelixStepper.hh"
#include "Geant4/G4ChordFinder.hh"
#include "Geant4/G4FieldManager.hh"

#include "Geant4/G4UserLimits.hh"

#include "Geant4/G4SDManager.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
using namespace std;

namespace mu2e {

  void constructMSTM( const VolumeInfo& parent,
                      const SimpleConfig& _config
                      ){

    MaterialFinder materialFinder(_config);

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "mstm", "mstm" );
    geomOptions->loadEntry( _config, "mstmMagnetField", "mstm.magnet.field" );
    
    int const verbosityLevel            = _config.getInt("mstm.verbosityLevel",0);
    const bool ismstmVisible            = geomOptions->isVisible("mstm"); 
    const bool ismstmSolid              = geomOptions->isSolid("mstm"); 
    const bool ismstmMagnetFieldVisible = geomOptions->isVisible("mstmMagnetField"); 
    const bool forceAuxEdgeVisible      = geomOptions->forceAuxEdgeVisible("mstm"); 
    const bool doSurfaceCheck           = geomOptions->doSurfaceCheck("mstm");
    const bool placePV                  = geomOptions->placePV("mstm");
    
    if ( verbosityLevel > 0 ) {
      std::cout << __func__ << " Constructing MSTM..." << std::endl;
    }

    // Fetch parent (hall) position
    G4ThreeVector parentCenterInMu2e = parent.centerInMu2e();

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    CLHEP::Hep3Vector const & dsP ( ds->position() );

    G4ThreeVector zeroVector(0.,0.,0.);


//     GeomHandle<DetectorSolenoidShielding> dss;

    //Create a reference position (everything in MSTM will be defined w.r.t. this position)
    
    G4ThreeVector mstmReferencePositionInMu2e(dsP.x(), 
                                              0.0, 
                                              _config.getDouble("mstm.refz0InMu2e") );
    
    G4ThreeVector mstmReferencePositionInParent = mstmReferencePositionInMu2e - parentCenterInMu2e;
    
    //----- Create the Mother volume for everything in the MSTM area--------------------------------
    
    //We want the Mother and the Shielding Wall to go down to the floor, so get the necessary info:
    const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );    //std::cout << " yExtentLow = " << yExtentLow << std::endl;
    
    const double mstmMotherHalfHeight =  fabs(yExtentLow);
    const double mstmMotherHalfWidth  =  _config.getDouble("mstm.wallUpStr.halfWidth");
    const double mstmMotherHalfLength =   1.0 * _config.getDouble("mstm.pipe0.halfLength")
                                        + 0.5 * _config.getDouble("mstm.collimator1.UpStrSpace")
                                        + 1.0 * _config.getDouble("mstm.collimator1.halfLength")
                                        + 0.5 * _config.getDouble("mstm.collimator2.UpStrSpace")
                                        + 1.0 * _config.getDouble("mstm.collimator2.halfLength")
                                        + 0.5 * _config.getDouble("mstm.shutter.UpStrSpace")
                                        + 1.0 * 500.0 //over-estimate 1 meter for the shutter length
                                        + 1.0 * _config.getDouble("mstm.pipe1.halfLength")
                                        + 0.5 * _config.getDouble("mstm.collimator3.UpStrSpace")
                                        + 1.0 * _config.getDouble("mstm.collimator3.halfLength")
                                        + 0.5 * _config.getDouble("mstm.can.UpStrSpace")
                                        + 1.0 * _config.getDouble("mstm.can.halfLength")
                                        + 0.5 * 1500.0; //make the Mother another 1.5m longer
                                        
    const double mstmMotherHalfLengths[3] = {mstmMotherHalfWidth, mstmMotherHalfHeight, mstmMotherHalfLength};
    
    std::string hallAirMaterialName = _config.getString("hall.insideMaterialName");
    G4Material* hallAirMaterial = findMaterialOrThrow(hallAirMaterialName);
    
    
    
    G4ThreeVector mstmMotherPositionInMu2e      = mstmReferencePositionInMu2e + G4ThreeVector(0.0, 0.0, mstmMotherHalfLength);
    G4ThreeVector mstmMotherPositionInParent    = mstmReferencePositionInParent + G4ThreeVector(0.0, 0.0, mstmMotherHalfLength);
    
    //  Make the mother volume for the MSTM.
    VolumeInfo mstmMotherInfo = nestBox("MSTMMother",
                                        mstmMotherHalfLengths,
                                        hallAirMaterial,  // Hall Air
                                        0x0,
                                        mstmMotherPositionInParent, //mstmMotherPositionInMu2e,
                                        parent,
                                        0,
                                        ismstmVisible,
                                        G4Color::Gray(),
                                        ismstmSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                       );
    
    if ( verbosityLevel > 0){
       cout << __func__ << " MSTM mother center in Mu2e   : " << mstmMotherPositionInMu2e << endl;
       cout << __func__ << " MSTM mother center in Parent : " << mstmMotherPositionInParent << endl;
       cout << __func__ << " parent.centerInMu2e()="<<parent.centerInMu2e() << endl;
       cout << __func__ << " mstmMotherInfo.centerInMu2e()="<<parent.centerInMu2e() << endl;
    }
    
    
    
    
    //----- Upstream shielding wall of MSTM area (2 ft thick? concrete wall)-------

    G4Material*  mstmUpStreamWallMaterial   = materialFinder.get("mstm.wallUpStr.material");
    const double mstmUpStreamWallUpStrSpace =  _config.getDouble("mstm.wallUpStr.UpStrSpace");
    const double mstmUpStreamWallHalfLength =  _config.getDouble("mstm.wallUpStr.halfLength");
    const double mstmUpStreamWallHalfWidth  =  _config.getDouble("mstm.wallUpStr.halfWidth");
    const double mstmUpStreamWallHoleROut   =  _config.getDouble("mstm.wallUpStr.holeRadius");

    // This box has a window.  implemented as a
    // G4SubtractionSolid to allow for another volume placement through it

    // Make the box for the wall
    G4Box* boxWallUpStream = new G4Box("boxWallUpStream",mstmUpStreamWallHalfWidth,fabs(yExtentLow),mstmUpStreamWallHalfLength);

    //Make the tube for the hole
    const TubsParams windparams(0.0,                        //inner radius
                                mstmUpStreamWallHoleROut,   //outer radius
                                mstmUpStreamWallHalfLength, //half length
                                0.0,                        //start angle
                                CLHEP::twopi                //end angle
                               );

    G4Tubs* windowTub = new G4Tubs( "window", 
                                   windparams.data()[0], 
                                   windparams.data()[1], 
                                   windparams.data()[2],
                                   windparams.data()[3], 
                                   windparams.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxWithWindow;
    boxWithWindow.name = "boxWallUpStreamWithWindow";
          
    // we need to put the window on the z axis
    G4ThreeVector mstmUpStreamWallPositionInMu2e   = mstmReferencePositionInMu2e + G4ThreeVector(0.0,0.0,+mstmUpStreamWallUpStrSpace + mstmUpStreamWallHalfLength);
    G4ThreeVector mstmUpStreamWallPositionInMother = mstmUpStreamWallPositionInMu2e - mstmMotherPositionInMu2e;
    //G4ThreeVector mstmUpStreamWallPositionInMother = G4ThreeVector(0.0,0.0,+mstmUpStreamWallUpStrSpace + mstmUpStreamWallHalfLength);
    
    //boxWithWindow.solid = new G4SubtractionSolid(boxWithWindow.name,boxWallUpStream,windowTub,0,mstmUpStreamWallPositionInMu2e);
    boxWithWindow.solid = new G4SubtractionSolid(boxWithWindow.name,boxWallUpStream,windowTub,0,zeroVector);    

    finishNesting(boxWithWindow,
                  mstmUpStreamWallMaterial,
                  0,
                  mstmUpStreamWallPositionInMother, 
                  mstmMotherInfo.logical, 
                  0,
                  ismstmVisible,
                  G4Colour::Magenta(),
                  ismstmSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);


    //----- Downstream shielding wall of MSTM area (2 ft thick? concrete wall)-------

    G4Material*  mstmDnStreamWallMaterial   = materialFinder.get("mstm.wallDnStr.material");
    const double mstmDnStreamWallHalfLength =  _config.getDouble("mstm.wallDnStr.halfLength");
    const double mstmColl1HalfWidth         =  _config.getDouble("mstm.collimator1.halfWidth");
    const double mstmBeamLeftWallHalfWidth  =  _config.getDouble("mstm.wallBeamLeft.halfWidth");
    const double mstmBeamRightWallHalfWidth =  _config.getDouble("mstm.wallBeamRight.halfWidth");
    // Make the box for the wall
    //const double mstmDnStreamWallHalfLengths[3] = {mstmUpStreamWallHalfWidth,
    const double mstmDnStreamWallHalfLengths[3] = {mstmColl1HalfWidth+mstmBeamLeftWallHalfWidth+mstmBeamRightWallHalfWidth,
                                                   fabs(yExtentLow),
                                                   mstmDnStreamWallHalfLength};
    //G4Box* boxWallDnStream = new G4Box("boxWallDnStream",mstmDnStreamWallHalfLengths[0],mstmDnStreamWallHalfLengths[1],mstmDnStreamWallHalfLengths[2]);

    G4ThreeVector mstmDnStreamWallPositionInMother(0.0,0.0,mstmMotherHalfLength - mstmDnStreamWallHalfLength);
    
    VolumeInfo mstmDnStreamWallInfo = nestBox("boxWallDnStream",
                                                  mstmDnStreamWallHalfLengths,
                                                  mstmDnStreamWallMaterial,
                                                  0x0,
                                                  mstmDnStreamWallPositionInMother,
                                                  mstmMotherInfo,
                                                  0,
                                                  ismstmVisible,
                                                  G4Color::Magenta(),
                                                  ismstmSolid,
                                                  forceAuxEdgeVisible,
                                                  placePV,
                                                  doSurfaceCheck
                                                  );


    //----- Beam Left shielding wall of MSTM area (2 ft thick? concrete wall)-------

    G4Material*  mstmBeamLeftWallMaterial  = materialFinder.get("mstm.wallBeamLeft.material");
    //const double mstmBeamLeftWallHalfWidth =  _config.getDouble("mstm.wallBeamLeft.halfWidth");
    const double mstmCeilingWallHalfHeight =  _config.getDouble("mstm.wallCeiling.halfHeight");
    //const double mstmColl1HalfWidth        =  _config.getDouble("mstm.collimator1.halfWidth");
    
    // Make the box for the wall
    const double mstmBeamLeftWallHalfLengths[3] = {mstmBeamLeftWallHalfWidth,
                                                   fabs(yExtentLow) - mstmCeilingWallHalfHeight,
                                                   mstmMotherHalfLength - 0.5*mstmUpStreamWallUpStrSpace - mstmUpStreamWallHalfLength - mstmDnStreamWallHalfLength};
    //G4Box* boxWallBeamLeft = new G4Box("boxWallBeamLeft",mstmBeamLeftWallHalfLengths[0],mstmBeamLeftWallHalfLengths[1],mstmBeamLeftWallHalfLengths[2]);

    //G4ThreeVector mstmBeamLeftWallPositionInMother = zeroVector + G4ThreeVector(-1.0*(mstmMotherHalfWidth-mstmBeamLeftWallHalfLengths[0]),
    G4ThreeVector mstmBeamLeftWallPositionInMother = zeroVector + G4ThreeVector(-1.0*(mstmColl1HalfWidth+mstmBeamLeftWallHalfLengths[0]),
                                                                                -1.0*mstmCeilingWallHalfHeight,
                                                                                0.5*mstmUpStreamWallUpStrSpace + mstmUpStreamWallHalfLength - mstmDnStreamWallHalfLength);
    
    VolumeInfo mstmBeamLeftWallInfo = nestBox("boxWallBeamLeft",
                                                  mstmBeamLeftWallHalfLengths,
                                                  mstmBeamLeftWallMaterial,
                                                  0x0,
                                                  mstmBeamLeftWallPositionInMother,
                                                  mstmMotherInfo,
                                                  0,
                                                  ismstmVisible,
                                                  G4Color::Magenta(),
                                                  ismstmSolid,
                                                  forceAuxEdgeVisible,
                                                  placePV,
                                                  doSurfaceCheck
                                                  );
    

    //----- Beam Right shielding wall of MSTM area (2 ft thick? concrete wall)-------

    G4Material*  mstmBeamRightWallMaterial  = materialFinder.get("mstm.wallBeamRight.material");
    //const double mstmBeamRightWallHalfWidth =  _config.getDouble("mstm.wallBeamRight.halfWidth");
    //const double mstmCeilingWallHalfHeight =  _config.getDouble("mstm.wallCeiling.halfHeight");

    // Make the box for the wall
    const double mstmBeamRightWallHalfLengths[3] = {mstmBeamRightWallHalfWidth,
                                                   fabs(yExtentLow) - mstmCeilingWallHalfHeight,
                                                   mstmMotherHalfLength - 0.5*mstmUpStreamWallUpStrSpace - mstmUpStreamWallHalfLength - mstmDnStreamWallHalfLength};
    //G4Box* boxWallBeamRight = new G4Box("boxWallBeamRight",mstmBeamRightWallHalfLengths[0],mstmBeamRightWallHalfLengths[1],mstmBeamRightWallHalfLengths[2]);

    //G4ThreeVector mstmBeamRightWallPositionInMother = zeroVector + G4ThreeVector(1.0*(mstmMotherHalfWidth-mstmBeamRightWallHalfLengths[0]),
    G4ThreeVector mstmBeamRightWallPositionInMother = zeroVector + G4ThreeVector(1.0*(mstmColl1HalfWidth+mstmBeamRightWallHalfLengths[0]),                                                   
                                                                                -1.0*mstmCeilingWallHalfHeight,
                                                                                0.5*mstmUpStreamWallUpStrSpace + mstmUpStreamWallHalfLength - mstmDnStreamWallHalfLength);
    
    VolumeInfo mstmBeamRightWallInfo = nestBox("boxWallBeamRight",
                                                  mstmBeamRightWallHalfLengths,
                                                  mstmBeamRightWallMaterial,
                                                  0x0,
                                                  mstmBeamRightWallPositionInMother,
                                                  mstmMotherInfo,
                                                  0,
                                                  ismstmVisible,
                                                  G4Color::Magenta(),
                                                  ismstmSolid,
                                                  forceAuxEdgeVisible,
                                                  placePV,
                                                  doSurfaceCheck
                                                  );

    
    //----- Ceiling shielding wall of MSTM area (2 ft thick? concrete wall)-------

    G4Material*  mstmCeilingWallMaterial  = materialFinder.get("mstm.wallCeiling.material");
    //const double mstmCeilingWallHalfHeight =  _config.getDouble("mstm.wallCeiling.halfHeight");

    // Make the box for the wall
    const double mstmCeilingWallHalfLengths[3] = {mstmDnStreamWallHalfLengths[0],
                                                  mstmCeilingWallHalfHeight,
                                                  mstmMotherHalfLength - 0.5*mstmUpStreamWallUpStrSpace - mstmUpStreamWallHalfLength - mstmDnStreamWallHalfLength};

    G4ThreeVector mstmCeilingWallPositionInMother = zeroVector + G4ThreeVector(0.0,//1.0*(mstmMotherHalfWidth-mstmCeilingWallHalfLengths[0]),
                                                                               mstmMotherHalfHeight - mstmCeilingWallHalfHeight,
                                                                               0.5*mstmUpStreamWallUpStrSpace + mstmUpStreamWallHalfLength - mstmDnStreamWallHalfLength);
    
    VolumeInfo mstmCeilingWallInfo = nestBox("boxWallCeiling",
                                                  mstmCeilingWallHalfLengths,
                                                  mstmCeilingWallMaterial,
                                                  0x0,
                                                  mstmCeilingWallPositionInMother,
                                                  mstmMotherInfo,
                                                  0,
                                                  ismstmVisible,
                                                  G4Color::Magenta(),
                                                  ismstmSolid,
                                                  forceAuxEdgeVisible,
                                                  placePV,
                                                  doSurfaceCheck
                                                  );    
    
    //----- Magnet ----------------------------
    
    //Just use a block of material for now (maybe stainless steel?, specified in fcl configuration)
    
    G4Material*  mstmMagnetMaterial       = materialFinder.get("mstm.magnet.material");
    const double mstmMagnetUpStrSpace     =  _config.getDouble("mstm.magnet.UpStrSpace");
    const double mstmMagnetHalfLength     =  _config.getDouble("mstm.magnet.halfLength");
    const double mstmMagnetHalfWidth      =  _config.getDouble("mstm.magnet.halfWidth");
    const double mstmMagnetHalfHeight     =  _config.getDouble("mstm.magnet.halfHeight");
    const double mstmMagnetHoleHalfHeight =  _config.getDouble("mstm.magnet.holeHalfHeight");
    const double mstmMagnetHoleHalfWidth  =  _config.getDouble("mstm.magnet.holeHalfWidth");
    
    //TODO: Throw if mstmMagnetHalfHeight is larger than distance to the floor.

    // This box has a window.  implemented as a
    // G4SubtractionSolid to alow for another volume placement
    // through it
    
    // Make the magnet
    G4Box* boxMagnet  = new G4Box("boxMagnet",mstmMagnetHalfWidth,mstmMagnetHalfHeight,mstmMagnetHalfLength);
    // Make the rectangular window
    G4Box* windowRect = new G4Box( "window", mstmMagnetHoleHalfWidth, mstmMagnetHoleHalfHeight, mstmMagnetHalfLength);

    VolumeInfo boxWithRectWindow;
    boxWithRectWindow.name = "boxMagnetWithWindow";
          
    // We need to put both the magnet and its window on the z axis
    // For now just put the magnet 10cm downstream of the wall
    G4ThreeVector mstmMagnetPositionInMu2e   = mstmUpStreamWallPositionInMu2e + G4ThreeVector(0.0, 0.0, mstmUpStreamWallHalfLength+mstmMagnetHalfLength + mstmMagnetUpStrSpace);
    G4ThreeVector mstmMagnetPositionInMother = mstmMagnetPositionInMu2e - mstmMotherPositionInMu2e;
                     
    boxWithRectWindow.solid = new G4SubtractionSolid(boxWithRectWindow.name,boxMagnet,windowRect,0,zeroVector);

    finishNesting(boxWithRectWindow,
                  mstmMagnetMaterial,
                  0,
                  mstmMagnetPositionInMother, 
                  mstmMotherInfo.logical,
                  0,
                  ismstmVisible,
                  G4Colour::Magenta(),
                  ismstmSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);
    
    // Make another rectangular volume for the magnetic field
    const double mstmMagnetHoleHalfLengths[3] = {mstmMagnetHoleHalfWidth,
                                                 mstmMagnetHoleHalfHeight,
                                                 mstmMagnetHalfLength};
    //G4Box* magRect = new G4Box( "mstmMagneticField", mstmMagnetHoleHalfWidth, mstmMagnetHoleHalfHeight, mstmMagnetHalfLength);
    VolumeInfo mstmMagneticFieldBoxInfo = nestBox("mstmMagneticField",
                                                  mstmMagnetHoleHalfLengths,
                                                  hallAirMaterial, //findMaterialOrThrow("G4_Galactic"),
                                                  0x0,
                                                  mstmMagnetPositionInMother,
                                                  mstmMotherInfo,
                                                  0,
                                                  ismstmMagnetFieldVisible,            //magnet visible
                                                  G4Color::Blue(),
                                                  false,           //ismstmSolid (this is just a field, not a solid)
                                                  forceAuxEdgeVisible,
                                                  placePV,         //must be true
                                                  doSurfaceCheck
                                                  );    
    
    // Create a magnetic field inside the window (hole) of the magnet box
    // Note the local values for the stepper etc...
    // Geant4 should take ownership of the objects created here

    const double mstmMagnetField = _config.getDouble("mstm.magnet.field");

    G4MagneticField        *localMagField        = new G4UniformMagField(G4ThreeVector(mstmMagnetField*CLHEP::tesla,0.0,0.0));//This makes negatively charged particles go towards the floor
    G4Mag_EqRhs            *MagRHS               = new G4Mag_UsualEqRhs(localMagField);
    G4MagIntegratorStepper *localMagStepper      = new G4ExactHelixStepper(MagRHS); // we use a specialized stepper
    G4ChordFinder          *localMagChordFinder  = new G4ChordFinder(localMagField,1.0e-2*CLHEP::mm,localMagStepper);
    G4FieldManager         *localMagFieldManager = new G4FieldManager(localMagField,localMagChordFinder,false);// pure magnetic filed does not change energy

    mstmMagneticFieldBoxInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
    
    G4UserLimits* mstmMagStepLimit = new G4UserLimits(5.*CLHEP::mm);
    mstmMagneticFieldBoxInfo.logical->SetUserLimits(mstmMagStepLimit);    
    

    
    //----- pipe0 -----------------------------
    // This pipe starts right after the VD at the entrance of the MSTM area
    // and goes inside the upstream shielding wall and inside the magnet.

    G4Material*  mstmPipe0Material              = materialFinder.get("mstm.pipe0.material");
    G4Material*  mstmPipe0Gas                   = materialFinder.get("mstm.pipe0.gas");
    const double mstmPipe0RIn                   =  _config.getDouble("mstm.pipe0.rIn");     
    const double mstmPipe0ROut                  =  _config.getDouble("mstm.pipe0.rOut");      
    //const double mstmPipe0HalfLength            =  _config.getDouble("mstm.pipe0.halfLength");
    const double mstmPipe0HalfLength            =  _config.getDouble("mstm.magnet.halfLength");//must be same size so we can place it inside the magnetic field volume 
    G4Material*  mstmPipe0UpStrWindowMaterial   = materialFinder.get("mstm.pipe0.UpStrWindowMaterial");
    const double mstmPipe0UpStrWindowHalfLength =  _config.getDouble("mstm.pipe0.UpStrWindowHalfLength");
    G4Material*  mstmPipe0DnStrWindowMaterial   = materialFinder.get("mstm.pipe0.DnStrWindowMaterial");
    const double mstmPipe0DnStrWindowHalfLength =  _config.getDouble("mstm.pipe0.DnStrWindowHalfLength");
    
    //Throw if pipe is too big to fit through the holes in the wall and magnet
    if ( mstmPipe0ROut > mstmUpStreamWallHoleROut ){
        throw cet::exception("GEOM")<< " MSTM: Pipe0 radius is too big to fit through upstream shielding wall. \n" ;
    }
    if ( mstmPipe0ROut > mstmMagnetHoleHalfWidth ){
      throw cet::exception("GEOM")<< " MSTM: Pipe0 radius is too big to fit through Magnet width. \n" ;
    }
    if ( mstmPipe0ROut > mstmMagnetHoleHalfHeight ){
      throw cet::exception("GEOM")<< " MSTM: Pipe0 radius is too big to fit through Magnet height. \n" ;
    }

    //TODO: Throw if pipe is not longer than the Wall+Magnet length
    
    const TubsParams mstmPipe0Params(0., mstmPipe0ROut, mstmPipe0HalfLength);
    const TubsParams mstmPipe0GasParams(0.,  mstmPipe0RIn, mstmPipe0HalfLength - 2.0*mstmPipe0UpStrWindowHalfLength - 2.0*mstmPipe0DnStrWindowHalfLength);
    const TubsParams mstmPipe0UpStrWindowParams(0., mstmPipe0RIn, mstmPipe0UpStrWindowHalfLength);
    const TubsParams mstmPipe0DnStrWindowParams(0., mstmPipe0RIn, mstmPipe0DnStrWindowHalfLength);
    
//     G4ThreeVector mstmPipe0PositionInMu2e   = mstmReferencePositionInMu2e + G4ThreeVector(0.0,0.0,mstmUpStreamWallUpStrSpace + mstmPipe0HalfLength);
//     G4ThreeVector mstmPipe0PositionInMother = mstmPipe0PositionInMu2e - mstmMotherPositionInMu2e;
    G4ThreeVector mstmPipe0PositionInMu2e   = mstmReferencePositionInMu2e + G4ThreeVector(0.0,0.0,0.0);
    G4ThreeVector mstmPipe0PositionInMother = mstmPipe0PositionInMu2e - mstmMagnetPositionInMu2e;
    
    VolumeInfo mstmPipe0Info = nestTubs( "mstmPipe0",
                                         mstmPipe0Params,
                                         mstmPipe0Material,
                                         0x0,
                                         G4ThreeVector(0.0,0.0,0.0), //we put the pipe centered on the magnetic field
                                         mstmMagneticFieldBoxInfo, //mstmMotherInfo,
                                         0,
                                         ismstmVisible,
                                         G4Color::Red(),
                                         ismstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe0GasInfo = nestTubs( "mstmPipe0Gas",
                                            mstmPipe0GasParams,
                                            mstmPipe0Gas,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0,2.0*(mstmPipe0UpStrWindowHalfLength-mstmPipe0DnStrWindowHalfLength)),
                                            mstmPipe0Info,
                                            0,
                                            ismstmVisible,
                                            G4Color::Yellow(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    VolumeInfo mstmPipe0UpStrWindowInfo = nestTubs( "mstmPipe0UpStrWindow",
                                            mstmPipe0UpStrWindowParams,
                                            mstmPipe0UpStrWindowMaterial,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0,-1.0*mstmPipe0HalfLength + mstmPipe0UpStrWindowHalfLength),
                                            mstmPipe0Info,
                                            0,
                                            ismstmVisible,
                                            G4Color::Red(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    
    VolumeInfo mstmPipe0DnStrWindowInfo = nestTubs( "mstmPipe0DnStrWindow",
                                            mstmPipe0DnStrWindowParams,
                                            mstmPipe0DnStrWindowMaterial,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0, mstmPipe0HalfLength - mstmPipe0DnStrWindowHalfLength),
                                            mstmPipe0Info,
                                            0,
                                            ismstmVisible,
                                            G4Color::Red(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    
    
    //----- Collimator 1 --------------------------------------------------------
    
    G4Material*  mstmColl1Material   = materialFinder.get("mstm.collimator1.material");
    const double mstmColl1UpStrSpace =  _config.getDouble("mstm.collimator1.UpStrSpace");
    const double mstmColl1HalfLength =  _config.getDouble("mstm.collimator1.halfLength");
    //const double mstmColl1HalfWidth  =  _config.getDouble("mstm.collimator1.halfWidth");
    const double mstmColl1HalfHeight =  _config.getDouble("mstm.collimator1.halfHeight");
    const double mstmColl1HoleROut   =  _config.getDouble("mstm.collimator1.holeRadius");

    // This box has a window.  implemented as a
    // G4SubtractionSolid to alow for another volume placement
    // through it

    // Make the box for the collimator
    G4Box* boxColl1 = new G4Box("boxColl1",mstmColl1HalfWidth,mstmColl1HalfHeight,mstmColl1HalfLength);

    //Make the tube for the hole
    const TubsParams windColl1Params(0.0,                 //inner radius
                                     mstmColl1HoleROut,   //outer raius
                                     mstmColl1HalfLength, //
                                     0.0,                 //start angle
                                     CLHEP::twopi         //end angle
                                    );

    G4Tubs* windowColl1 = new G4Tubs( "window", 
                                   windColl1Params.data()[0], 
                                   windColl1Params.data()[1], 
                                   windColl1Params.data()[2],
                                   windColl1Params.data()[3], 
                                   windColl1Params.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxColl1WithWindow;
    boxColl1WithWindow.name = "boxColl1WithWindow";
          
    // We need to put the window on the z axis
    // Leave a gap between the end of the pipe and the collimator
    G4ThreeVector mstmColl1PositionInMu2e   = mstmMagnetPositionInMu2e + G4ThreeVector(0.0,0.0, mstmMagnetHalfLength + mstmColl1UpStrSpace + mstmColl1HalfLength);
    G4ThreeVector mstmColl1PositionInMother = mstmColl1PositionInMu2e - mstmMotherPositionInMu2e;
        
    
    boxColl1WithWindow.solid = new G4SubtractionSolid(boxColl1WithWindow.name,boxColl1,windowColl1,0,zeroVector);

    finishNesting(boxColl1WithWindow,
                  mstmColl1Material,
                  0,
                  mstmColl1PositionInMother,
                  mstmMotherInfo.logical,
                  0,
                  ismstmVisible,
                  G4Colour::Magenta(),
                  ismstmSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    
    //----- shutter of variable number of segments  -----------------------------
    
    const int mstmShutterNumberSegments = _config.getInt("mstm.shutter.numberSegments");
    const double mstmShutterUpStrSpace  = _config.getDouble("mstm.shutter.UpStrSpace");
    const double mstmShutterHalfHeight  = _config.getDouble("mstm.shutter.halfHeight");
    double mstmShutterSegmentLastHalfLength = mstmShutterUpStrSpace; //this starts as the inital offset between Coll1 and the Shutter
    G4ThreeVector mstmShutterSegmentPositionInMu2e   = mstmColl1PositionInMu2e + G4ThreeVector(0.0,0.0,mstmColl1HalfLength);
    G4ThreeVector mstmShutterSegmentPositionInMother = mstmShutterSegmentPositionInMu2e - mstmMotherPositionInMu2e;
    
    for (int segment = 1; segment <= mstmShutterNumberSegments; ++segment) {
      std::stringstream material_config_name, halfLength_config_name;
      material_config_name << "mstm.shutter.segment" << segment << ".material";
      halfLength_config_name << "mstm.shutter.segment" << segment << ".halfLength";
      
      G4Material* mstmShutterSegmentMaterial = materialFinder.get(material_config_name.str());

      const double mstmShutterSegmentHalfLength = _config.getDouble(halfLength_config_name.str());
      const double mstmShutterSegmentHalfLengths[3] = {mstmShutterHalfHeight,
						       mstmShutterHalfHeight,
						       mstmShutterSegmentHalfLength};
      mstmShutterSegmentPositionInMu2e   += G4ThreeVector(0.0, 0.0, mstmShutterSegmentLastHalfLength + mstmShutterSegmentHalfLength);
      mstmShutterSegmentPositionInMother += G4ThreeVector(0.0, 0.0, mstmShutterSegmentLastHalfLength + mstmShutterSegmentHalfLength);
      
      std::stringstream volume_info_name;
      volume_info_name << "mstmShutterSegment" << segment;
      VolumeInfo mstmShutterSegmentInfo = nestBox(volume_info_name.str(),
						  mstmShutterSegmentHalfLengths,
						  mstmShutterSegmentMaterial,
						  0x0,
						  mstmShutterSegmentPositionInMother,
						  mstmMotherInfo,
						  0,
						  ismstmVisible,
						  G4Color::Yellow(), //Gray
						  ismstmSolid,
						  forceAuxEdgeVisible,
						  placePV,
						  doSurfaceCheck
						  );
      mstmShutterSegmentLastHalfLength = mstmShutterSegmentHalfLength;
    }

    
    //----- Collimator 2 --------------------------------------------------------
    
    G4Material*  mstmColl2Material   = materialFinder.get("mstm.collimator2.material");
    const double mstmColl2UpStrSpace =  _config.getDouble("mstm.collimator2.UpStrSpace");
    const double mstmColl2HalfLength =  _config.getDouble("mstm.collimator2.halfLength");
    const double mstmColl2HalfWidth  =  _config.getDouble("mstm.collimator2.halfWidth");
    const double mstmColl2HalfHeight =  _config.getDouble("mstm.collimator2.halfHeight");
    const double mstmColl2HoleROut   =  _config.getDouble("mstm.collimator2.holeRadius");

    // This box has a window.  implemented as a
    // G4SubtractionSolid to alow for another volume placement
    // through it

    // Make the box for the collimator
    G4Box* boxColl2 = new G4Box("boxColl2",mstmColl2HalfWidth,mstmColl2HalfHeight,mstmColl2HalfLength);

    //Make the tube for the hole
    const TubsParams windColl2Params(0.0,                 //inner radius
                                     mstmColl2HoleROut,   //outer raius
                                     mstmColl2HalfLength, //
                                     0.0,                 //start angle
                                     CLHEP::twopi         //end angle                                     
                                    );

    G4Tubs* windowColl2 = new G4Tubs( "window", 
                                   windColl2Params.data()[0], 
                                   windColl2Params.data()[1], 
                                   windColl2Params.data()[2],
                                   windColl2Params.data()[3], 
                                   windColl2Params.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxColl2WithWindow;
    boxColl2WithWindow.name = "boxColl2WithWindow";
          
    // We need to put the window on the z axis
    // Leave a gap between the this and the previous element
    G4ThreeVector mstmColl2PositionInMu2e   = mstmShutterSegmentPositionInMu2e + G4ThreeVector(0.0,0.0,mstmShutterSegmentLastHalfLength) + G4ThreeVector(0.0,0.0, mstmColl2UpStrSpace + mstmColl2HalfLength);
    //G4ThreeVector mstmColl2PositionInMother   = mstmShutterSegmentPositionInMother + G4ThreeVector(0.0,0.0,mstmShutterSegmentLastHalfLength) + G4ThreeVector(0.0,0.0, mstmColl2UpStrSpace + mstmColl2HalfLength);
    G4ThreeVector mstmColl2PositionInMother = mstmColl2PositionInMu2e - mstmMotherPositionInMu2e;
        
    boxColl2WithWindow.solid = new G4SubtractionSolid(boxColl2WithWindow.name,boxColl2,windowColl2,0,zeroVector);

    finishNesting(boxColl2WithWindow,
                  mstmColl2Material,
                  0,
                  mstmColl2PositionInMother,
                  mstmMotherInfo.logical,
                  0,
                  ismstmVisible,
                  G4Colour::Magenta(),
                  ismstmSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    
    //----- Pipe 1 -----------------------------------------------------------
    // This pipe goes between collimators 2 and 3.

    G4Material*  mstmPipe1Material              = materialFinder.get("mstm.pipe1.material");
    G4Material*  mstmPipe1Gas                   = materialFinder.get("mstm.pipe1.gas");
    const double mstmPipe1RIn                   =  _config.getDouble("mstm.pipe1.rIn");     
    const double mstmPipe1ROut                  =  _config.getDouble("mstm.pipe1.rOut");      
    const double mstmPipe1UpStrSpace            =  _config.getDouble("mstm.pipe1.UpStrSpace");
    const double mstmPipe1HalfLength            =  _config.getDouble("mstm.pipe1.halfLength");
    G4Material*  mstmPipe1UpStrWindowMaterial   = materialFinder.get("mstm.pipe1.UpStrWindowMaterial");
    const double mstmPipe1UpStrWindowHalfLength =  _config.getDouble("mstm.pipe1.UpStrWindowHalfLength");
    G4Material*  mstmPipe1DnStrWindowMaterial   = materialFinder.get("mstm.pipe1.DnStrWindowMaterial");
    const double mstmPipe1DnStrWindowHalfLength =  _config.getDouble("mstm.pipe1.DnStrWindowHalfLength");
    
    //TODO: Throw if pipe is too big to fit through the holes in the wall
    
    const TubsParams mstmPipe1Params(0., mstmPipe1ROut, mstmPipe1HalfLength);
    const TubsParams mstmPipe1GasParams(0.,  mstmPipe1RIn, mstmPipe1HalfLength - 2.0*mstmPipe1UpStrWindowHalfLength - 2.0*mstmPipe1DnStrWindowHalfLength);
    const TubsParams mstmPipe1UpStrWindowParams(0., mstmPipe1RIn, mstmPipe1UpStrWindowHalfLength);
    const TubsParams mstmPipe1DnStrWindowParams(0., mstmPipe1RIn, mstmPipe1DnStrWindowHalfLength);
    
    // Leave a gap between this and the previous element
    G4ThreeVector mstmPipe1PositionInMu2e   = mstmColl2PositionInMu2e + G4ThreeVector(0.0,0.0,mstmColl2HalfLength + mstmPipe1UpStrSpace + mstmPipe1HalfLength);
    G4ThreeVector mstmPipe1PositionInMother = mstmPipe1PositionInMu2e - mstmMotherPositionInMu2e;
    
    VolumeInfo mstmPipe1Info = nestTubs( "mstmPipe1",
                                         mstmPipe1Params,
                                         mstmPipe1Material,
                                         0x0,
                                         mstmPipe1PositionInMother,
                                         mstmMotherInfo,
                                         0,
                                         ismstmVisible,
                                         G4Color::Red(),
                                         ismstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe1GasInfo = nestTubs( "mstmPipe1Gas",
                                            mstmPipe1GasParams,
                                            mstmPipe1Gas,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0,2.0*(mstmPipe1UpStrWindowHalfLength-mstmPipe1DnStrWindowHalfLength)),
                                            mstmPipe1Info,
                                            0,
                                            ismstmVisible,
                                            G4Color::Yellow(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    VolumeInfo mstmPipe1UpStrWindowInfo = nestTubs( "mstmPipe1UpStrWindow",
                                            mstmPipe1UpStrWindowParams,
                                            mstmPipe1UpStrWindowMaterial,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0,-1.0*mstmPipe1HalfLength + mstmPipe1UpStrWindowHalfLength),
                                            mstmPipe1Info,
                                            0,
                                            ismstmVisible,
                                            G4Color::Red(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    
    VolumeInfo mstmPipe1DnStrWindowInfo = nestTubs( "mstmPipe1DnStrWindow",
                                            mstmPipe1DnStrWindowParams,
                                            mstmPipe1DnStrWindowMaterial,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0, mstmPipe1HalfLength - mstmPipe1DnStrWindowHalfLength),
                                            mstmPipe1Info,
                                            0,
                                            ismstmVisible,
                                            G4Color::Red(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    

    //----- Collimator 3 --------------------------------------------------------
    
    G4Material*  mstmColl3Material   = materialFinder.get("mstm.collimator3.material");
    const double mstmColl3UpStrSpace =  _config.getDouble("mstm.collimator3.UpStrSpace");
    const double mstmColl3HalfLength =  _config.getDouble("mstm.collimator3.halfLength");
    const double mstmColl3HalfWidth  =  _config.getDouble("mstm.collimator3.halfWidth");
    const double mstmColl3HalfHeight =  _config.getDouble("mstm.collimator3.halfHeight");
    const double mstmColl3HoleROut   =  _config.getDouble("mstm.collimator3.holeRadius");

    // This box has a window.  implemented as a
    // G4SubtractionSolid to alow for another volume placement through it

    // Make the box for the collimator
    G4Box* boxColl3 = new G4Box("boxColl3",mstmColl3HalfWidth,mstmColl3HalfHeight,mstmColl3HalfLength);

    //Make the tube for the hole
    const TubsParams windColl3Params(0.0,                 //inner radius
                                     mstmColl3HoleROut,   //outer raius
                                     mstmColl3HalfLength, //
                                     0.0,                 //start angle
                                     CLHEP::twopi         //end angle                                     
                                    );

    G4Tubs* windowColl3 = new G4Tubs( "window", 
                                   windColl3Params.data()[0], 
                                   windColl3Params.data()[1], 
                                   windColl3Params.data()[2],
                                   windColl3Params.data()[3], 
                                   windColl3Params.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxColl3WithWindow;
    boxColl3WithWindow.name = "boxColl3WithWindow";
          
    // We need to put the window on the z axis
    // Leave a gap between the this and the next element
    G4ThreeVector mstmColl3PositionInMu2e   = mstmPipe1PositionInMu2e + G4ThreeVector(0.0,0.0,mstmPipe1HalfLength + mstmColl3UpStrSpace + mstmColl3HalfLength);
    G4ThreeVector mstmColl3PositionInMother = mstmColl3PositionInMu2e - mstmMotherPositionInMu2e;
        
    boxColl3WithWindow.solid = new G4SubtractionSolid(boxColl3WithWindow.name,boxColl3,windowColl3,0,zeroVector);

    finishNesting(boxColl3WithWindow,
                  mstmColl3Material,
                  0,
                  mstmColl3PositionInMother, /*  should is be mstmColl3PositionInMu2e ??? */
                  mstmMotherInfo.logical,
                  0,
                  ismstmVisible,
                  G4Colour::Magenta(),
                  ismstmSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);


    //----- Can (Surrounds Ge Crystal) --------------------------------------------------
    // (it has it's plugs "inside" the main part)

    G4Material*  mstmCanMaterial              = materialFinder.get("mstm.can.material");
    G4Material*  mstmCanGas                   = materialFinder.get("mstm.can.gas");
    const double mstmCanRIn                   =  _config.getDouble("mstm.can.rIn");     
    const double mstmCanROut                  =  _config.getDouble("mstm.can.rOut");      
    const double mstmCanUpStrSpace            =  _config.getDouble("mstm.can.UpStrSpace");
    const double mstmCanHalfLength            =  _config.getDouble("mstm.can.halfLength");
    G4Material*  mstmCanUpStrWindowMaterial   = materialFinder.get("mstm.can.UpStrWindowMaterial");
    const double mstmCanUpStrWindowHalfLength =  _config.getDouble("mstm.can.UpStrWindowHalfLength");
    G4Material*  mstmCanDnStrWindowMaterial   = materialFinder.get("mstm.can.DnStrWindowMaterial");
    const double mstmCanDnStrWindowHalfLength =  _config.getDouble("mstm.can.DnStrWindowHalfLength");

    
    const TubsParams mstmCanParams(0., mstmCanROut, mstmCanHalfLength);
    const TubsParams mstmCanGasParams(0.,  mstmCanRIn, mstmCanHalfLength - 2.0*mstmCanUpStrWindowHalfLength - 2.0*mstmCanDnStrWindowHalfLength);
    const TubsParams mstmCanUpStrWindowParams(0., mstmCanRIn, mstmCanUpStrWindowHalfLength);
    const TubsParams mstmCanDnStrWindowParams(0., mstmCanRIn, mstmCanDnStrWindowHalfLength);
    
    // Leave a gap between this and the previous element
    G4ThreeVector mstmCanPositionInMu2e   = mstmColl3PositionInMu2e + G4ThreeVector(0.0,0.0,mstmColl3HalfLength + mstmCanUpStrSpace + mstmCanHalfLength);
    G4ThreeVector mstmCanPositionInMother = mstmCanPositionInMu2e - mstmMotherPositionInMu2e;
    
    VolumeInfo mstmCanInfo = nestTubs( "mstmCan",
                                         mstmCanParams,
                                         mstmCanMaterial,
                                         0x0,
                                         mstmCanPositionInMother,
                                         mstmMotherInfo,
                                         0,
                                         ismstmVisible,
                                         G4Color::Green(),
                                         ismstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

        VolumeInfo mstmCanGasInfo = nestTubs( "mstmCanGas",
                                                mstmCanGasParams,
                                                mstmCanGas,
                                                0x0,
                                                zeroVector + G4ThreeVector(0.0,0.0,2.0*(mstmCanUpStrWindowHalfLength-mstmCanDnStrWindowHalfLength)),
                                                mstmCanInfo,
                                                0,
                                                ismstmVisible,
                                                G4Color::Yellow(),
                                                ismstmSolid,
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck
                                                );

    VolumeInfo mstmCanUpStrWindowInfo = nestTubs( "mstmCanUpStrWindow",
                                            mstmCanUpStrWindowParams,
                                            mstmCanUpStrWindowMaterial,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0,-1.0*mstmCanHalfLength + mstmCanUpStrWindowHalfLength),
                                            mstmCanInfo,
                                            0,
                                            ismstmVisible,
                                            G4Color::Green(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    
    VolumeInfo mstmCanDnStrWindowInfo = nestTubs( "mstmCanDnStrWindow",
                                            mstmCanDnStrWindowParams,
                                            mstmCanDnStrWindowMaterial,
                                            0x0,
                                            zeroVector + G4ThreeVector(0.0,0.0, mstmCanHalfLength - mstmCanDnStrWindowHalfLength),
                                            mstmCanInfo,
                                            0,
                                            ismstmVisible,
                                            G4Color::Green(),
                                            ismstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    
    //----- the sensitive detector ("crystal") -------------------------------------------

    G4Material*  mstmCrystalMaterial   = materialFinder.get("mstm.crystal.material");
    //const double mstmCrystalRIn        =  _config.getDouble("mstm.crystal.rIn");
    const double mstmCrystalROut       =  _config.getDouble("mstm.crystal.rOut");
    const double mstmCrystalHalfLength =  _config.getDouble("mstm.crystal.halfLength");

    //TODO: Throw if crystal doesn't fit inside the can
    if ( mstmCrystalROut > mstmCanRIn ){
      throw cet::exception("GEOM")<< " MSTM: Crystal radius is larger than the inner radius of the can. \n" ;
    }
    if ( mstmCrystalHalfLength > mstmCanHalfLength-mstmCanDnStrWindowHalfLength-mstmCanUpStrWindowHalfLength ){
      throw cet::exception("GEOM")<< " MSTM: Crystal length is larger than the inner length of the can. \n" ;
    }
    
    const TubsParams mstmCrystalParams(0., mstmCrystalROut, mstmCrystalHalfLength);

    G4ThreeVector mstmCrystalPositionInMu2e   = mstmCanPositionInMu2e;
    G4ThreeVector mstmCrystalPositionInMother = mstmCrystalPositionInMu2e - mstmMotherPositionInMu2e;
    
    VolumeInfo mstmCrystal = nestTubs("mstmCrystal",
                                      mstmCrystalParams,
                                      mstmCrystalMaterial,
                                      0x0,
                                      zeroVector,   //mstmCrystalPositionInMother,
                                      mstmCanGasInfo,  //Put the Crystal in the gas=vacuum inside the can
                                      0,
                                      ismstmVisible,
                                      G4Color::Red(),
                                      ismstmSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    if ( verbosityLevel > 0) {
      std::cout << __func__ << " mstmMotherPositionInMu2e = " << mstmMotherPositionInMu2e << endl;       
      std::cout << __func__ << " mstmReferencePositionInMu2e = " << mstmReferencePositionInMu2e << endl; 
    }

  }

} // end of constructMSTM;
