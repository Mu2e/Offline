//
// Free function to create MSTM.
// Muon Stopping Target Monitor
//
// $Id: constructMSTM.cc,v 1.4 2014/07/29 16:23:24 genser Exp $
// $Author: genser $
// $Date: 2014/07/29 16:23:24 $
//
// Original author K.L.Genser 
//

// Mu2e includes.
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "ExternalNeutronShieldingGeom/inc/ExtNeutShieldCendBoxes.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructMSTM.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "GeometryService/inc/VirtualDetector.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

#include "G4SDManager.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

#include <vector>

using namespace std;

namespace mu2e {

  void constructMSTM( const VolumeInfo& parent,
                      SimpleConfig const & _config
                      ){
    MaterialFinder materialFinder(_config);

    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    // Fetch parent (hall) position
    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    int const verbosityLevel = _config.getInt("mstm.verbosityLevel",0);

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    CLHEP::Hep3Vector const & dsP ( ds->position() );

    G4ThreeVector zeroVector(0.,0.,0.);

    const bool mstmVisible = _config.getBool("mstm.visible", true );
    const bool mstmSolid   = _config.getBool("mstm.solid",   false);

    // pipe0

    const double mstmPipe0RIn        = _config.getDouble("mstm.pipe0.rIn");     
    const double mstmPipe0ROut       = _config.getDouble("mstm.pipe0.rOut");      
    const double mstmPipe0HalfLength = _config.getDouble("mstm.pipe0.halfLength");

    G4Material* mstmPipe0Material  = materialFinder.get("mstm.pipe0.material");
    G4Material* mstmPipe0Gas       = materialFinder.get("mstm.pipe0.gas");

    const TubsParams mstmPipe0Params(0., mstmPipe0ROut, mstmPipe0HalfLength);
    const TubsParams mstmPipe0GasParams(0.,  mstmPipe0RIn, mstmPipe0HalfLength);
      
    GeomHandle<DetectorSolenoidShielding> dss;

    // place the MSTM wrt to ENS; account for the vd half lenght
    GeomHandle<ExtNeutShieldCendBoxes> enscendb;

    const std::vector<CLHEP::Hep3Vector>& ENSCBcentersOfBoxes = enscendb->centersOfBoxes();

    size_t nBox = ENSCBcentersOfBoxes.size();
    size_t ib;
    for(ib = 0; ib < nBox; ++ib) {
      if ( enscendb->hasHole(ib) ) break;
    }
    int hID = enscendb->holeIndex(ib);
    if ( verbosityLevel > 0 ) {
      std::cout << __func__ << " ib: " << ib << std::endl;
    }

    GeomHandle<VirtualDetector> vd;

    // locations are wrt HallAir
    // for some reason the location has to be taken from the box not the hole tbd
    //        CLHEP::Hep3Vector holeLocation = enscendb->holeLocation(hID);
    CLHEP::Hep3Vector holeLocation = ENSCBcentersOfBoxes[ib];

    G4ThreeVector mstmPipe0PositionInMu2e(dsP.x(), dsP.y(),
                                          holeLocation.z() + enscendb->holeHalfLength(hID) +
                                          2.*vd->getHalfLength() + mstmPipe0HalfLength);

    VolumeInfo mstmPipe0Info = nestTubs( "mstmPipe0",
                                         mstmPipe0Params,
                                         mstmPipe0Material,
                                         0x0,
                                         mstmPipe0PositionInMu2e - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe0GasInfo = nestTubs( "mstmPipe0Gas",
                                            mstmPipe0GasParams,
                                            mstmPipe0Gas,
                                            0x0,
                                            zeroVector,
                                            mstmPipe0Info,
                                            0,
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );


    // pipe1

    const double mstmPipe1RIn        = _config.getDouble("mstm.pipe1.rIn");     
    const double mstmPipe1ROut       = _config.getDouble("mstm.pipe1.rOut");      
    const double mstmPipe1HalfLength = _config.getDouble("mstm.pipe1.halfLength");

    G4Material* mstmPipe1Material  = materialFinder.get("mstm.pipe1.material");
    G4Material* mstmPipe1Gas       = materialFinder.get("mstm.pipe1.gas");

    const TubsParams mstmPipe1Params(0., mstmPipe1ROut, mstmPipe1HalfLength);
    const TubsParams mstmPipe1GasParams(0.,  mstmPipe1RIn, mstmPipe1HalfLength);

    G4ThreeVector mstmPipe1PositionInMu2e =  mstmPipe0PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe0HalfLength + mstmPipe1HalfLength);

    VolumeInfo mstmPipe1Info = nestTubs( "mstmPipe1",
                                         mstmPipe1Params,
                                         mstmPipe1Material,
                                         0x0,
                                         mstmPipe1PositionInMu2e - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Magenta(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe1GasInfo = nestTubs( "mstmPipe1Gas",
                                            mstmPipe1GasParams,
                                            mstmPipe1Gas,
                                            0x0,
                                            zeroVector,
                                            mstmPipe1Info,
                                            0,
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );


    // creating magnetic field for the mstmPipe1 hierarchy
    // note the local values for the stepper etc...
    // Geant4 should take ownership of the objects created here

    const double mstmPipe1MagField = _config.getDouble("mstm.pipe1.magfield");

    G4MagneticField *localPipe1MagField = 
      new G4UniformMagField(G4ThreeVector(mstmPipe1MagField*CLHEP::tesla,0.0,0.0));
    G4Mag_EqRhs *Pipe1RHS  = new G4Mag_UsualEqRhs(localPipe1MagField);
    G4MagIntegratorStepper *localPipe1Stepper = new G4ExactHelixStepper(Pipe1RHS); // we use a specialized stepper
    G4ChordFinder *localPipe1ChordFinder = new G4ChordFinder(localPipe1MagField,1.0e-2*CLHEP::mm,localPipe1Stepper);
    G4FieldManager *localPipe1FieldManager    
      = new G4FieldManager(localPipe1MagField,localPipe1ChordFinder,false);// pure magnetic filed does not change energy

    mstmPipe1Info.logical->SetFieldManager(localPipe1FieldManager, true); // propagate it down the hierarchy

    // absorber1

    G4Material* mstmAbsorber1Material = materialFinder.get("mstm.absorber1.material");

    const double mstmAbsorber1HalfHeight    = _config.getDouble("mstm.absorber1.halfHeight");
    const double mstmAbsorber1HalfLength    = _config.getDouble("mstm.absorber1.halfLength");

    const double mstmAbsorber1HalfLengths[3] = {mstmAbsorber1HalfHeight,
                                                mstmAbsorber1HalfHeight, 
                                                mstmAbsorber1HalfLength};

    G4ThreeVector mstmAbsorber1PositionInMu2e = mstmPipe1PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe1HalfLength + mstmAbsorber1HalfLength);

    VolumeInfo mstmAbsorber1Info = nestBox("mstmAbsorber1",
                                           mstmAbsorber1HalfLengths,
                                           mstmAbsorber1Material,
                                           0x0,
                                           mstmAbsorber1PositionInMu2e - _hallOriginInMu2e,
                                           parent,
                                           0,
                                           mstmVisible,
                                           G4Color::Gray(),
                                           mstmSolid,
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );



    // pipe2

    const double mstmPipe2RIn        = _config.getDouble("mstm.pipe2.rIn");     
    const double mstmPipe2ROut       = _config.getDouble("mstm.pipe2.rOut");      
    const double mstmPipe2HalfLength = _config.getDouble("mstm.pipe2.halfLength");

    G4Material* mstmPipe2Material  = materialFinder.get("mstm.pipe2.material");
    G4Material* mstmPipe2Gas       = materialFinder.get("mstm.pipe2.gas");

    const TubsParams mstmPipe2Params(0., mstmPipe2ROut, mstmPipe2HalfLength);
    const TubsParams mstmPipe2GasParams(0.,  mstmPipe2RIn, mstmPipe2HalfLength);
      
    G4ThreeVector mstmPipe2PositionInMu2e =  mstmAbsorber1PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmAbsorber1HalfLength + mstmPipe2HalfLength);

    VolumeInfo mstmPipe2Info = nestTubs( "mstmPipe2",
                                         mstmPipe2Params,
                                         mstmPipe2Material,
                                         0x0,
                                         mstmPipe2PositionInMu2e - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe2GasInfo = nestTubs( "mstmPipe2Gas",
                                            mstmPipe2GasParams,
                                            mstmPipe2Gas,
                                            0x0,
                                            zeroVector,
                                            mstmPipe2Info,
                                            0,
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );


    // boxes for pipe3 (they have the length of pipe3 as they are all placed one in the other)
    // and pipe3

    G4Material* mstmPipe3BoxInMaterial  = materialFinder.get("mstm.boxIn.material");

    const double mstmPipe3BoxInHalfHeight     = _config.getDouble("mstm.boxIn.halfHeight");

    G4Material* mstmPipe3BoxOutMaterial = materialFinder.get("mstm.boxOut.material");

    const double mstmPipe3BoxOutHalfHeight    = _config.getDouble("mstm.boxOut.halfHeight");

    const double mstmPipe3RIn        = _config.getDouble("mstm.pipe3.rIn");     
    const double mstmPipe3ROut       = _config.getDouble("mstm.pipe3.rOut");      
    const double mstmPipe3HalfLength = _config.getDouble("mstm.pipe3.halfLength");

    G4Material* mstmPipe3Material  = materialFinder.get("mstm.pipe3.material");
    G4Material* mstmPipe3Gas       = materialFinder.get("mstm.pipe3.gas");

    // constructing boxes

    const double mstmPipe3BoxOutHalfLengths[3] = {mstmPipe3BoxOutHalfHeight, 
                                                  mstmPipe3BoxOutHalfHeight, 
                                                  mstmPipe3HalfLength};

    const double mstmPipe3BoxInHalfLengths[3] = {mstmPipe3BoxInHalfHeight, 
                                                 mstmPipe3BoxInHalfHeight, 
                                                 mstmPipe3HalfLength};

    const TubsParams mstmPipe3Params(0., mstmPipe3ROut, mstmPipe3HalfLength);
    const TubsParams mstmPipe3GasParams(0., mstmPipe3RIn,  mstmPipe3HalfLength);

    G4ThreeVector mstmPipe3PositionInMu2e = mstmPipe2PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe2HalfLength + mstmPipe3HalfLength);

    VolumeInfo mstmPipe3BoxOutInfo = nestBox("mstmPipe3BoxOut",
                                             mstmPipe3BoxOutHalfLengths,
                                             mstmPipe3BoxOutMaterial,
                                             0x0,
                                             mstmPipe3PositionInMu2e - _hallOriginInMu2e,
                                             parent,
                                             0,
                                             mstmVisible,
                                             G4Color::Gray(),
                                             mstmSolid,
                                             forceAuxEdgeVisible,
                                             placePV,
                                             doSurfaceCheck
                                             );

    VolumeInfo mstmPipe3BoxInInfo = nestBox("mstmPipe3BoxIn",
                                            mstmPipe3BoxInHalfLengths,
                                            mstmPipe3BoxInMaterial,
                                            0x0,
                                            zeroVector,
                                            mstmPipe3BoxOutInfo,
                                            0,
                                            mstmVisible,
                                            G4Color::Cyan(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    // constructing pipe3

    VolumeInfo mstmPipe3Info = nestTubs( "mstmPipe3",
                                         mstmPipe3Params,
                                         mstmPipe3Material,
                                         0x0,
                                         zeroVector,
                                         mstmPipe3BoxInInfo,
                                         0,
                                         mstmVisible,
                                         G4Color::Magenta(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    // creating magnetic field for the mstmPipe3 hierarchy
    // note the local values for the stepper etc...
    // Geant4 should take ownership of the objects created here

    const double mstmPipe3MagField = _config.getDouble("mstm.pipe3.magfield");

    G4MagneticField *localPipe3MagField = 
      new G4UniformMagField(G4ThreeVector(mstmPipe3MagField*CLHEP::tesla,0.0,0.0));
    G4Mag_EqRhs *Pipe3RHS  = new G4Mag_UsualEqRhs(localPipe3MagField);
    G4MagIntegratorStepper *localPipe3Stepper = new G4ExactHelixStepper(Pipe3RHS); // we use a specialized stepper
    G4ChordFinder *localPipe3ChordFinder = new G4ChordFinder(localPipe3MagField,1.0e-2*CLHEP::mm,localPipe3Stepper);
    G4FieldManager *localPipe3FieldManager    
      = new G4FieldManager(localPipe3MagField,localPipe3ChordFinder,false);// pure magnetic filed does not change energy

    mstmPipe3Info.logical->SetFieldManager(localPipe3FieldManager, true); // propagate it down the hierarchy

    VolumeInfo mstmPipe3GasInfo = nestTubs("mstmPipe3Gas",
                                           mstmPipe3GasParams,
                                           mstmPipe3Gas,
                                           0x0,
                                           zeroVector,
                                           mstmPipe3Info,
                                           0,
                                           mstmVisible,
                                           G4Color::Yellow(),
                                           mstmSolid,
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );
    
    // window, placed inside the pipe; we need to place it inside the
    // gas at the end of it to avoid volume overlaps

    const double mstmOutWindowRIn        = _config.getDouble("mstm.outWindow.rIn");
    const double mstmOutWindowROut       = _config.getDouble("mstm.outWindow.rOut");
    const double mstmOutWindowHalfLength = _config.getDouble("mstm.outWindow.halfLength");

    G4Material* mstmOutWindowMaterial = materialFinder.get("mstm.outWindow.material");

    const TubsParams mstmOutWindowParams(mstmOutWindowRIn, mstmOutWindowROut, mstmOutWindowHalfLength);
    G4ThreeVector mstmOutWindowPositionInMu2e = mstmPipe3PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe3HalfLength - mstmOutWindowHalfLength);

    VolumeInfo mstmOutWindow = nestTubs("mstmOutWindow",
                                        mstmOutWindowParams,
                                        mstmOutWindowMaterial,
                                        0x0,
                                        mstmOutWindowPositionInMu2e - mstmPipe3PositionInMu2e,
                                        mstmPipe3GasInfo,
                                        0,
                                        mstmVisible,
                                        G4Color::Red(),
                                        mstmSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    // pipe4

    const double mstmPipe4RIn        = _config.getDouble("mstm.pipe4.rIn");     
    const double mstmPipe4ROut       = _config.getDouble("mstm.pipe4.rOut");      
    const double mstmPipe4HalfLength = _config.getDouble("mstm.pipe4.halfLength");

    G4Material* mstmPipe4Material    = materialFinder.get("mstm.pipe4.material");
    G4Material* mstmPipe4Gas         = materialFinder.get("mstm.pipe4.gas");

    const TubsParams mstmPipe4Params(0., mstmPipe4ROut, mstmPipe4HalfLength);
    const TubsParams mstmPipe4GasParams(0.,  mstmPipe4RIn, mstmPipe4HalfLength);

    G4ThreeVector mstmPipe4PositionInMu2e = mstmPipe3PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe3HalfLength + mstmPipe4HalfLength);


    VolumeInfo mstmPipe4Info = nestTubs( "mstmPipe4",
                                         mstmPipe4Params,
                                         mstmPipe4Material,
                                         0x0,
                                         mstmPipe4PositionInMu2e - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe4GasInfo = nestTubs( "mstmPipe4Gas",
                                            mstmPipe4GasParams,
                                            mstmPipe4Gas,
                                            0x0,
                                            zeroVector,
                                            mstmPipe4Info,
                                            0,
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );


    // absorber2

    G4Material* mstmAbsorber2Material = materialFinder.get("mstm.absorber2.material");

    const double mstmAbsorber2HalfHeight    = _config.getDouble("mstm.absorber2.halfHeight");
    const double mstmAbsorber2HalfLength    = _config.getDouble("mstm.absorber2.halfLength");

    const double mstmAbsorber2HalfLengths[3] = {mstmAbsorber2HalfHeight, 
                                                mstmAbsorber2HalfHeight, 
                                                mstmAbsorber2HalfLength};

    G4ThreeVector mstmAbsorber2PositionInMu2e = mstmPipe4PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe4HalfLength + mstmAbsorber2HalfLength);

    VolumeInfo mstmAbsorber2Info = nestBox("mstmAbsorber2",
                                           mstmAbsorber2HalfLengths,
                                           mstmAbsorber2Material,
                                           0x0,
                                           mstmAbsorber2PositionInMu2e - _hallOriginInMu2e,
                                           parent,
                                           0,
                                           mstmVisible,
                                           G4Color::Gray(),
                                           mstmSolid,
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );

    // pipe5

    const double mstmPipe5RIn        = _config.getDouble("mstm.pipe5.rIn");     
    const double mstmPipe5ROut       = _config.getDouble("mstm.pipe5.rOut");      
    const double mstmPipe5HalfLength = _config.getDouble("mstm.pipe5.halfLength");

    G4Material* mstmPipe5Material    = materialFinder.get("mstm.pipe5.material");
    G4Material* mstmPipe5Gas         = materialFinder.get("mstm.pipe5.gas");

    const TubsParams mstmPipe5Params(0., mstmPipe5ROut, mstmPipe5HalfLength);
    const TubsParams mstmPipe5GasParams(0.,  mstmPipe5RIn, mstmPipe5HalfLength);

    G4ThreeVector mstmPipe5PositionInMu2e = mstmAbsorber2PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmAbsorber2HalfLength + mstmPipe5HalfLength);


    VolumeInfo mstmPipe5Info = nestTubs( "mstmPipe5",
                                         mstmPipe5Params,
                                         mstmPipe5Material,
                                         0x0,
                                         mstmPipe5PositionInMu2e - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    VolumeInfo mstmPipe5GasInfo = nestTubs( "mstmPipe5Gas",
                                            mstmPipe5GasParams,
                                            mstmPipe5Gas,
                                            0x0,
                                            zeroVector,
                                            mstmPipe5Info,
                                            0,
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );




    // can
    // (it has it's plugs "inside" the main part)

    G4Material* mstmCanMaterial = materialFinder.get("mstm.can.material");

    const double mstmCanUpSRIn            = _config.getDouble("mstm.canUpS.rIn");
    const double mstmCanUpSROut           = _config.getDouble("mstm.canUpS.rOut");
    const double mstmCanUpSHalfLength     = _config.getDouble("mstm.canUpS.halfLength");

    const double mstmCanRIn               = _config.getDouble("mstm.can.rIn");
    const double mstmCanROut              = _config.getDouble("mstm.can.rOut");
    const double mstmCanHalfLength        = _config.getDouble("mstm.can.halfLength");

    const double mstmCanDownSRIn          = _config.getDouble("mstm.canDownS.rIn");
    const double mstmCanDownSROut         = _config.getDouble("mstm.canDownS.rOut");
    const double mstmCanDownSHalfLength   = _config.getDouble("mstm.canDownS.halfLength");

    // the can upstream part
 
    const TubsParams mstmCanUpSParams(mstmCanUpSRIn, mstmCanUpSROut, mstmCanUpSHalfLength);

    G4ThreeVector mstmCanUpSPositionInMu2e = mstmPipe5PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe5HalfLength + mstmCanUpSHalfLength);

    VolumeInfo mstmCanUpS = nestTubs("mstmCanUpS",
                                     mstmCanUpSParams,
                                     mstmCanMaterial,
                                     0x0,
                                     mstmCanUpSPositionInMu2e - _hallOriginInMu2e,
                                     parent,
                                     0,
                                     mstmVisible,
                                     G4Color::Blue(),
                                     mstmSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    // the main can

    const TubsParams mstmCanParams(mstmCanRIn, mstmCanROut, mstmCanHalfLength);

    G4ThreeVector mstmCanPositionInMu2e = mstmPipe5PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe5HalfLength + mstmCanHalfLength);

    VolumeInfo mstmCan = nestTubs("mstmCan",
                                  mstmCanParams,
                                  mstmCanMaterial,
                                  0x0,
                                  mstmCanPositionInMu2e - _hallOriginInMu2e,
                                  parent,
                                  0,
                                  mstmVisible,
                                  G4Color::Blue(),
                                  mstmSolid,
                                  forceAuxEdgeVisible,
                                  placePV,
                                  doSurfaceCheck
                                  );

    // the can downstream part
 
    const TubsParams mstmCanDownSParams(mstmCanDownSRIn, mstmCanDownSROut, mstmCanDownSHalfLength);

    G4ThreeVector mstmCanDownSPositionInMu2e = mstmPipe5PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe5HalfLength + 2.0*mstmCanHalfLength - mstmCanDownSHalfLength);

    VolumeInfo mstmCanDownS = nestTubs("mstmCanDownS",
                                       mstmCanDownSParams,
                                       mstmCanMaterial,
                                       0x0,
                                       mstmCanDownSPositionInMu2e - _hallOriginInMu2e,
                                       parent,
                                       0,
                                       mstmVisible,
                                       G4Color::Blue(),
                                       mstmSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // the detector ("crystal")

    G4Material* mstmCrystalMaterial = materialFinder.get("mstm.crystal.material");

    const double mstmCrystalRIn           = _config.getDouble("mstm.crystal.rIn");
    const double mstmCrystalROut          = _config.getDouble("mstm.crystal.rOut");
    const double mstmCrystalHalfLength    = _config.getDouble("mstm.crystal.halfLength");

    const TubsParams mstmCrystalParams(mstmCrystalRIn, mstmCrystalROut, mstmCrystalHalfLength);

    G4ThreeVector mstmCrystalPositionInMu2e = mstmCanPositionInMu2e;

    VolumeInfo mstmCrystal = nestTubs("mstmCrystal",
                                      mstmCrystalParams,
                                      mstmCrystalMaterial,
                                      0x0,
                                      mstmCrystalPositionInMu2e - _hallOriginInMu2e,
                                      parent,
                                      0,
                                      mstmVisible,
                                      G4Color::Yellow(),
                                      mstmSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );


    // Make mstmCrystal a sensitive detector.
    G4VSensitiveDetector *sd = 
      G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::MSTMCrystal());
    if(sd) mstmCrystal.logical->SetSensitiveDetector(sd);

    if ( verbosityLevel > 0) {
      std::cout << __func__ << " mstm position:"
                << " dss->getIFBendPlug()->zEnd(): " 
                << dss->getIFBendPlug()->zEnd()
                << " z position: "
                << dss->getIFBendPlug()->zEnd() + mstmPipe0HalfLength <<  std::endl;
      std::cout << __func__ << " ens hole end " 
                << " holeLocation.z() + enscendb->holeHalfLength(hID): " 
                << holeLocation.z() + enscendb->holeHalfLength(hID) 
                << " z position: "
                << holeLocation.z() + enscendb->holeHalfLength(hID) + mstmPipe0HalfLength << std::endl;

      std::cout << __func__ << " dss->getIFBendPlug()->originInMu2e() " 
                << dss->getIFBendPlug()->originInMu2e() <<  std::endl;
      std::cout << __func__ << " dss->getIFBendPlug() position in hall " 
                << dss->getIFBendPlug()->originInMu2e() - _hallOriginInMu2e <<  std::endl;
      std::cout << __func__ << " dsp position in Mu2e " << dsP 
                << " _hallOriginInMu2e " <<  _hallOriginInMu2e <<  std::endl;
      std::cout << __func__ << " dsp position in hall " << dsP - _hallOriginInMu2e <<  std::endl;

      std::cout << __func__ << " mstmPipe0PositionInMu2e " << 
        mstmPipe0PositionInMu2e <<  std::endl;
      std::cout << __func__ << " mstmPipe0PositionInHall " << 
        mstmPipe0PositionInMu2e - _hallOriginInMu2e <<  std::endl;

      std::cout << __func__ << " mstmPipe2PositionInMu2e " << 
        mstmPipe2PositionInMu2e <<  std::endl;
      std::cout << __func__ << " mstmPipe2PositionInHall " << 
        mstmPipe2PositionInMu2e - _hallOriginInMu2e <<  std::endl;

      std::cout << __func__ << " mstmCanPositionInMu2e " << 
        mstmCanPositionInMu2e <<  std::endl;
      std::cout << __func__ << " mstmCanPositionInHall " << 
        mstmCanPositionInMu2e - _hallOriginInMu2e <<  std::endl;
    }

  }

} // end of constructMSTM;
