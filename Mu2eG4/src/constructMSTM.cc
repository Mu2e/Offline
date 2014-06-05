//
// Free function to create MSTM.
// Muon Stopping Target Monitor
//
// $Id: constructMSTM.cc,v 1.1 2014/06/05 21:05:32 genser Exp $
// $Author: genser $
// $Date: 2014/06/05 21:05:32 $
//
// Original author K.L.Genser 
//

// Mu2e includes.
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
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

using namespace std;

namespace mu2e {

  void constructMSTM( const VolumeInfo& parent,
                      SimpleConfig const & _config
                      ){
    MaterialFinder materialFinder(_config);

    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    // bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool doSurfaceCheck      = true;
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

    // pipe2

    //      double mstmPipe2Z          = _config.getDouble("mstm.pipe2.z");
    const double mstmPipe2RIn        = _config.getDouble("mstm.pipe2.rIn");     
    const double mstmPipe2ROut       = _config.getDouble("mstm.pipe2.rOut");      
    const double mstmPipe2HalfLength = _config.getDouble("mstm.pipe2.halfLength");

    G4Material* mstmPipe2Material  = materialFinder.get("mstm.pipe2.material");
    G4Material* mstmPipe2Gas       = materialFinder.get("mstm.pipe2.gas");

    const double mstmPipe2Params[5]     = {0., mstmPipe2ROut, mstmPipe2HalfLength, 0.0, CLHEP::twopi};
    const double mstmPipe2GasParams[5]  = {0.,  mstmPipe2RIn, mstmPipe2HalfLength, 0.0, CLHEP::twopi};
      
    GeomHandle<DetectorSolenoidShielding> dss;

    G4ThreeVector mstmPipe2PositionInMu2e(dsP.x(), dsP.y(),
                                          dss->getIFBendPlug()->zEnd() + mstmPipe2HalfLength);

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


    if ( verbosityLevel > 0) {
      std::cout << " mstm position:"
                << " dss->getIFBendPlug()->zEnd(): " 
                << dss->getIFBendPlug()->zEnd()
                << " z position: " 
                << dss->getIFBendPlug()->zEnd() + mstmPipe2HalfLength
                << "\n dss->getIFBendPlug()->originInMu2e() " 
                << dss->getIFBendPlug()->originInMu2e()
                << "\n dss->getIFBendPlug() position in hall " 
                << dss->getIFBendPlug()->originInMu2e() - _hallOriginInMu2e
                << "\n dsp position in Mu2e " << dsP 
                << " _hallOriginInMu2e " <<  _hallOriginInMu2e
                << "\n dsp position in hall " << dsP - _hallOriginInMu2e
        //                  << " ds3positionInMu2e " << ds3positionInMu2e 
                << " mstmPipe2PositionInMu2e " << mstmPipe2PositionInMu2e
                << " mstmPipe2PositionInHall " << mstmPipe2PositionInMu2e - _hallOriginInMu2e
                <<  std::endl;
    }

    // boxes for pipe3

    G4Material* mstmPipe3BoxInMaterial  = materialFinder.get("mstm.boxIn.material");

    const double mstmPipe3BoxInHalfHeight     = _config.getDouble("mstm.boxIn.halfHeight");
    const double mstmPipe3BoxInHalfLength     = _config.getDouble("mstm.boxIn.halfLength");

    G4Material* mstmPipe3BoxOutMaterial = materialFinder.get("mstm.boxOut.material");

    const double mstmPipe3BoxOutHalfHeight    = _config.getDouble("mstm.boxOut.halfHeight");
    const double mstmPipe3BoxOutHalfLength    = _config.getDouble("mstm.boxOut.halfLength");

    const double mstmPipe3BoxOutHalfLengths[3] = {mstmPipe3BoxOutHalfHeight, 
                                            mstmPipe3BoxOutHalfHeight, 
                                            mstmPipe3BoxOutHalfLength};

    const double mstmPipe3BoxInHalfLengths[3] = {mstmPipe3BoxInHalfHeight, 
                                           mstmPipe3BoxInHalfHeight, 
                                           mstmPipe3BoxInHalfLength};

    // pipe3

    const double mstmPipe3RIn        = _config.getDouble("mstm.pipe3.rIn");     
    const double mstmPipe3ROut       = _config.getDouble("mstm.pipe3.rOut");      
    const double mstmPipe3HalfLength = _config.getDouble("mstm.pipe3.halfLength");

    G4Material* mstmPipe3Material  = materialFinder.get("mstm.pipe3.material");
    G4Material* mstmPipe3Gas       = materialFinder.get("mstm.pipe3.gas");

    const double mstmPipe3Params[5]     = {0., mstmPipe3ROut, mstmPipe3HalfLength, 0.0, CLHEP::twopi};
    const double mstmPipe3GasParams[5]  = {0., mstmPipe3RIn,  mstmPipe3HalfLength, 0.0, CLHEP::twopi};

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

    VolumeInfo mstmPipe3Info = nestTubs( "mstmPipe3",
                                         mstmPipe3Params,
                                         mstmPipe3Material,
                                         0x0,
                                         zeroVector,
                                         mstmPipe3BoxInInfo,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    // creating magnetic field for the mstmPipe3 hierarchy
    // note the local values for the stepper etc...
    // Geant4 should take ownership of the objects created here

    const double mstmPipe3MagField = _config.getDouble("mstm.pipe3.magfield");

    G4MagneticField *localMagField = 
      new G4UniformMagField(G4ThreeVector(mstmPipe3MagField*CLHEP::tesla,0.0,0.0));
    G4Mag_EqRhs *rhs  = new G4Mag_UsualEqRhs(localMagField);
    G4MagIntegratorStepper *stepper = new G4ExactHelixStepper(rhs); // we use a specialized stepper
    G4ChordFinder *chordFinder      = new G4ChordFinder(localMagField,1.0e-2*CLHEP::mm,stepper);
    G4FieldManager *localFieldManager    
      = new G4FieldManager(localMagField,chordFinder,false);// pure magnetic filed does not change energy

    mstmPipe3Info.logical->SetFieldManager(localFieldManager, true); // propagate it down the hierarchy

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
    // window

    const double mstmOutWindowRIn           = _config.getDouble("mstm.outWindow.rIn");
    const double mstmOutWindowROut          = _config.getDouble("mstm.outWindow.rOut");
    const double mstmOutWindowHalfLength    = _config.getDouble("mstm.outWindow.halfLength");

    G4Material* mstmOutWindowMaterial = materialFinder.get("mstm.outWindow.material");

    const double mstmOutWindowParams[5]  = {mstmOutWindowRIn, mstmOutWindowROut, mstmOutWindowHalfLength,
                                      0.0, CLHEP::twopi};

    G4ThreeVector mstmOutWindowPositionInMu2e = mstmPipe3PositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmPipe3HalfLength + mstmOutWindowHalfLength);

    VolumeInfo mstmOutWindow = nestTubs( "mstmOutWindow",
                                         mstmOutWindowParams,
                                         mstmOutWindowMaterial,
                                         0x0,
                                         mstmOutWindowPositionInMu2e - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    // absorber

    G4Material* mstmAbsorberMaterial = materialFinder.get("mstm.absorber.material");

    const double mstmAbsorberHalfHeight    = _config.getDouble("mstm.absorber.halfHeight");
    const double mstmAbsorberHalfLength    = _config.getDouble("mstm.absorber.halfLength");

    const double mstmAbsorberHalfLengths[3] = {mstmAbsorberHalfHeight, 
                                         mstmAbsorberHalfHeight, 
                                         mstmAbsorberHalfLength};

    G4ThreeVector mstmAbsorberPositionInMu2e = mstmOutWindowPositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmOutWindowHalfLength + mstmAbsorberHalfLength);

    VolumeInfo mstmAbsorberInfo = nestBox("mstmAbsorber",
                                          mstmAbsorberHalfLengths,
                                          mstmAbsorberMaterial,
                                          0x0,
                                          mstmAbsorberPositionInMu2e - _hallOriginInMu2e,
                                          parent,
                                          0,
                                          mstmVisible,
                                          G4Color::Gray(),
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
 
    const double mstmCanUpSParams[5] =  {mstmCanUpSRIn, mstmCanUpSROut, mstmCanUpSHalfLength,
                                   0.0, CLHEP::twopi};

    G4ThreeVector mstmCanUpSPositionInMu2e = mstmAbsorberPositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmAbsorberHalfLength + mstmCanUpSHalfLength);

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

    const double mstmCanParams[5] = {mstmCanRIn, mstmCanROut, mstmCanHalfLength,
                               0.0, CLHEP::twopi};

    G4ThreeVector mstmCanPositionInMu2e = mstmAbsorberPositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmAbsorberHalfLength + mstmCanHalfLength);

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
 
    const double mstmCanDownSParams[5] = {mstmCanDownSRIn, mstmCanDownSROut, mstmCanDownSHalfLength,
                                    0.0, CLHEP::twopi};

    G4ThreeVector mstmCanDownSPositionInMu2e = mstmAbsorberPositionInMu2e + 
      G4ThreeVector(0.0, 0.0, mstmAbsorberHalfLength + 2.0*mstmCanHalfLength - mstmCanDownSHalfLength);

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

    const double mstmCrystalParams[5] = {mstmCrystalRIn, mstmCrystalROut, mstmCrystalHalfLength,
                                   0.0, CLHEP::twopi};

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

  }

} // end of constructMSTM;
