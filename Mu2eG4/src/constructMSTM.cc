//
// Free function to create MSTM.
// Muon Stopping Target Monitor
//
// $Id: constructMSTM.cc,v 1.5 2014/09/16 21:58:31 jrquirk Exp $
// $Author: jrquirk $
// $Date: 2014/09/16 21:58:31 $
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
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

#include "G4SDManager.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <sstream>

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

    if ( verbosityLevel > 0 ) {
      std::cout << __func__ << " Constructing MSTM..." << std::endl;
    }


    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    CLHEP::Hep3Vector const & dsP ( ds->position() );

    G4ThreeVector zeroVector(0.,0.,0.);

    const bool mstmVisible = _config.getBool("mstm.visible", true );
    const bool mstmSolid   = _config.getBool("mstm.solid",   false);

    GeomHandle<DetectorSolenoidShielding> dss;

    // place the MSTM wrt to ENS; account for the vd half length
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
    
    //Create a reference position (everything in MSTM will be defined w.r.t. this position)
    G4ThreeVector mstmReferencePositionInMu2e(dsP.x(), 
                                              dsP.y(), 
                                              holeLocation.z() + enscendb->holeHalfLength(hID) + 2.*vd->getHalfLength());
    
    GeomHandle<Mu2eEnvelope> env;
        const CLHEP::Hep3Vector hallFormalCenterInMu2e(
                                                       (env->xmax() + env->xmin())/2.,
                                                       (env->ymax() + env->ymin())/2.,
                                                       (env->zmax() + env->zmin())/2.
                                                       );

    mstmReferencePositionInMu2e  = mstmReferencePositionInMu2e - hallFormalCenterInMu2e;
    
    
    //----- Upstream shielding wall of MSTM area (2 ft thick concrete wall?)-------

    //We want the shielding wall to go down to the floor, so get the necessary info
    GeomHandle<Mu2eBuilding> building;
    //const double yExtentLow = building->hallInsideYmin() - building->hallFloorThickness();
    const double yExtentLow = building->hallInsideYmin();
    
    G4Material*  mstmUpStreamWallMaterial   = materialFinder.get("mstm.wallUpStr.material");
    const double mstmUpStreamWallUpStrSpace =  _config.getDouble("mstm.wallUpStr.UpStrSpace");
    const double mstmUpStreamWallHalfLength =  _config.getDouble("mstm.wallUpStr.halfLength");
    const double mstmUpStreamWallHalfWidth  =  _config.getDouble("mstm.wallUpStr.halfWidth");
    const double mstmUpStreamWallHoleROut   =  _config.getDouble("mstm.wallUpStr.holeRadius");

    // This box has a window.  implemented as a
    // G4SubtractionSolid to alow for another volume placement
    // through it

    // Make the box for the wall
    G4Box* boxWallUpStream = new G4Box("boxWallUpStream",mstmUpStreamWallHalfWidth,fabs(yExtentLow),mstmUpStreamWallHalfLength);

    //Make the tube for the hole
    const TubsParams windparams(0.0,                        //inner radius
                                mstmUpStreamWallHoleROut,   //outer raius
                                mstmUpStreamWallHalfLength  //
                               );

    G4Tubs* windowTub = new G4Tubs( "window", 
                                   windparams.data()[0], 
                                   windparams.data()[1], 
                                   windparams.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                   windparams.data()[3], 
                                   windparams.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxWithWindow;
    boxWithWindow.name = "boxWallUpStreamWithWindow";
          
    // we need to put the window on the z axis
    G4ThreeVector mstmUpStreamWallPositionInMu2e = mstmReferencePositionInMu2e + G4ThreeVector(0.0,0.0,+mstmUpStreamWallUpStrSpace + mstmUpStreamWallHalfLength);
        
    boxWithWindow.solid = new G4SubtractionSolid(boxWithWindow.name,boxWallUpStream,windowTub,0,mstmUpStreamWallPositionInMu2e);

    finishNesting(boxWithWindow,
                  mstmUpStreamWallMaterial,
                  0,
                  mstmUpStreamWallPositionInMu2e, /*  should is be zeroVector or mstmUpStreamWallPositionInMu2e ??? */
                  parent.logical,
                  0,
                  _config.getBool("mstm.visible"),
                  G4Colour::Magenta(),
                  _config.getBool("mstm.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    
    //----- Magnet ----------------------------
    
    //Just use a block of material for now (maybe stainless steel?)
    
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
    G4ThreeVector mstmMagnetPositionInMu2e =  mstmUpStreamWallPositionInMu2e + G4ThreeVector(0.0, 0.0, mstmUpStreamWallHalfLength+mstmMagnetHalfLength + mstmMagnetUpStrSpace);
                     
    boxWithRectWindow.solid = new G4SubtractionSolid(boxWithRectWindow.name,boxMagnet,windowRect,0,mstmMagnetPositionInMu2e);

    finishNesting(boxWithRectWindow,
                  mstmMagnetMaterial,
                  0,
                  mstmMagnetPositionInMu2e, /*  should it be zeroVector or mstmMagnetPositionInMu2e ??? */
                  parent.logical,
                  0,
                  _config.getBool("mstm.visible"),
                  G4Colour::Magenta(),
                  _config.getBool("mstm.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);
    
    // Create a magnetic field inside the window (hole) of the magnet box
    // Note the local values for the stepper etc...
    // Geant4 should take ownership of the objects created here

    const double mstmMagnetField = _config.getDouble("mstm.magnet.field");

    G4MagneticField        *localMagField          = new G4UniformMagField(G4ThreeVector(mstmMagnetField*CLHEP::tesla,0.0,0.0));
    G4Mag_EqRhs            *MagRHS               = new G4Mag_UsualEqRhs(localMagField);
    G4MagIntegratorStepper *localMagStepper      = new G4ExactHelixStepper(MagRHS); // we use a specialized stepper
    G4ChordFinder          *localMagChordFinder  = new G4ChordFinder(localMagField,1.0e-2*CLHEP::mm,localMagStepper);
    G4FieldManager         *localMagFieldManager = new G4FieldManager(localMagField,localMagChordFinder,false);// pure magnetic filed does not change energy

    //WARNING: this puts the field in the whole size of the magnet
    //TODO:    change this so the field is just inside the window.
    boxWithRectWindow.logical->SetFieldManager(localMagFieldManager, true); // propagate it down the hierarchy

    
    //----- pipe0 -----------------------------
    // This pipe starts right after the VD at the entrance of the MSTM area
    // and goes inside the upstream shielding wall and inside the magnet.

    G4Material*  mstmPipe0Material              = materialFinder.get("mstm.pipe0.material");
    G4Material*  mstmPipe0Gas                   = materialFinder.get("mstm.pipe0.gas");
    const double mstmPipe0RIn                   =  _config.getDouble("mstm.pipe0.rIn");     
    const double mstmPipe0ROut                  =  _config.getDouble("mstm.pipe0.rOut");      
    const double mstmPipe0HalfLength            =  _config.getDouble("mstm.pipe0.halfLength");
    G4Material*  mstmPipe0UpStrWindowMaterial   = materialFinder.get("mstm.pipe0.UpStrWindowMaterial");
    const double mstmPipe0UpStrWindowHalfLength =  _config.getDouble("mstm.pipe0.UpStrWindowHalfLength");
    G4Material*  mstmPipe0DnStrWindowMaterial   = materialFinder.get("mstm.pipe0.DnStrWindowMaterial");
    const double mstmPipe0DnStrWindowHalfLength =  _config.getDouble("mstm.pipe0.DnStrWindowHalfLength");
    
    //TODO: Throw if pipe is too big to fit through the holes in the wall and magnet

    //TODO: Throw if pipe is not longer than the Wall+Magnet length
    
    const TubsParams mstmPipe0Params(0., mstmPipe0ROut, mstmPipe0HalfLength);
    const TubsParams mstmPipe0GasParams(0.,  mstmPipe0RIn, mstmPipe0HalfLength - 2.0*mstmPipe0UpStrWindowHalfLength - 2.0*mstmPipe0DnStrWindowHalfLength);
    const TubsParams mstmPipe0UpStrWindowParams(0., mstmPipe0RIn, mstmPipe0UpStrWindowHalfLength);
    const TubsParams mstmPipe0DnStrWindowParams(0., mstmPipe0RIn, mstmPipe0DnStrWindowHalfLength);
    
    G4ThreeVector mstmPipe0PositionInMu2e    = mstmReferencePositionInMu2e + G4ThreeVector(0.0,0.0,mstmUpStreamWallUpStrSpace + mstmPipe0HalfLength);
    
    VolumeInfo mstmPipe0Info = nestTubs( "mstmPipe0",
                                         mstmPipe0Params,
                                         mstmPipe0Material,
                                         0x0,
                                         mstmPipe0PositionInMu2e, // - _hallOriginInMu2e,
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
                                            zeroVector + G4ThreeVector(0.0,0.0,2.0*(mstmPipe0UpStrWindowHalfLength-mstmPipe0DnStrWindowHalfLength)),
                                            mstmPipe0Info,
                                            0,
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Red(),
                                            mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Red(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );
    
    
    //----- Collimator 1 --------------------------------------------------------
    
    G4Material*  mstmColl1Material   = materialFinder.get("mstm.collimator1.material");
    const double mstmColl1UpStrSpace =  _config.getDouble("mstm.collimator1.UpStrSpace");
    const double mstmColl1HalfLength =  _config.getDouble("mstm.collimator1.halfLength");
    const double mstmColl1HalfWidth  =  _config.getDouble("mstm.collimator1.halfWidth");
    const double mstmColl1HalfHeight =  _config.getDouble("mstm.collimator1.halfHeight");
    const double mstmColl1HoleROut   =  _config.getDouble("mstm.collimator1.holeRadius");

    // This box has a window.  implemented as a
    // G4SubtractionSolid to alow for another volume placement
    // through it

    // Make the box for the collimator
    G4Box* boxColl1 = new G4Box("boxColl1",mstmColl1HalfWidth,mstmColl1HalfHeight,mstmColl1HalfLength);

    //Make the tube for the hole
    const TubsParams windColl1Params(0.0,                        //inner radius
                                     mstmColl1HoleROut,   //outer raius
                                     mstmColl1HalfLength  //
                                    );

    G4Tubs* windowColl1 = new G4Tubs( "window", 
                                   windColl1Params.data()[0], 
                                   windColl1Params.data()[1], 
                                   windColl1Params.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                   windColl1Params.data()[3], 
                                   windColl1Params.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxColl1WithWindow;
    boxColl1WithWindow.name = "boxColl1WithWindow";
          
    // We need to put the window on the z axis
    // Leave a gap between the end of the pipe and the collimator
    G4ThreeVector mstmColl1PositionInMu2e = mstmMagnetPositionInMu2e + G4ThreeVector(0.0,0.0, mstmMagnetHalfLength + mstmColl1UpStrSpace + mstmColl1HalfLength);
        
    boxColl1WithWindow.solid = new G4SubtractionSolid(boxColl1WithWindow.name,boxColl1,windowColl1,0,mstmColl1PositionInMu2e);

    finishNesting(boxColl1WithWindow,
                  mstmColl1Material,
                  0,
                  mstmColl1PositionInMu2e, /*  should is be zeroVector or mstmColl1PositionInMu2e ??? */
                  parent.logical,
                  0,
                  _config.getBool("mstm.visible"),
                  G4Colour::Magenta(),
                  _config.getBool("mstm.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    
    //----- shutter of variable number of segments  -----------------------------
    
    const int mstmShutterNumberSegments = _config.getInt("mstm.shutter.numberSegments");
    const double mstmShutterUpStrSpace  = _config.getDouble("mstm.shutter.UpStrSpace");
    const double mstmShutterHalfHeight  = _config.getDouble("mstm.shutter.halfHeight");
    double mstmShutterSegmentLastHalfLength = mstmShutterUpStrSpace; //this starts as the inital offset between Coll1 and the Shutter
    G4ThreeVector mstmShutterSegmentPositionInMu2e = mstmColl1PositionInMu2e + G4ThreeVector(0.0,0.0,mstmColl1HalfLength);
    for (int segment = 1; segment <= mstmShutterNumberSegments; ++segment) {
      std::stringstream material_config_name, halfLength_config_name;
      material_config_name << "mstm.shutter.segment" << segment << ".material";
      halfLength_config_name << "mstm.shutter.segment" << segment << ".halfLength";
      
      G4Material* mstmShutterSegmentMaterial = materialFinder.get(material_config_name.str());

      const double mstmShutterSegmentHalfLength = _config.getDouble(halfLength_config_name.str());
      const double mstmShutterSegmentHalfLengths[3] = {mstmShutterHalfHeight,
						       mstmShutterHalfHeight,
						       mstmShutterSegmentHalfLength};
      mstmShutterSegmentPositionInMu2e += G4ThreeVector(0.0, 0.0, mstmShutterSegmentLastHalfLength + mstmShutterSegmentHalfLength);

      std::stringstream volume_info_name;
      volume_info_name << "mstmShutterSegment" << segment;
      VolumeInfo mstmShutterSegmentInfo = nestBox(volume_info_name.str(),
						  mstmShutterSegmentHalfLengths,
						  mstmShutterSegmentMaterial,
						  0x0,
						  mstmShutterSegmentPositionInMu2e,// - _hallOriginInMu2e,
						  parent,
						  0,
						  mstmVisible,
						  G4Color::Yellow(), //Gray
						  mstmSolid,
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
    const TubsParams windColl2Params(0.0,                        //inner radius
                                     mstmColl2HoleROut,   //outer raius
                                     mstmColl2HalfLength  //
                                    );

    G4Tubs* windowColl2 = new G4Tubs( "window", 
                                   windColl2Params.data()[0], 
                                   windColl2Params.data()[1], 
                                   windColl2Params.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                   windColl2Params.data()[3], 
                                   windColl2Params.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxColl2WithWindow;
    boxColl2WithWindow.name = "boxColl2WithWindow";
          
    // We need to put the window on the z axis
    // Leave a gap between the this and the previous element
    G4ThreeVector mstmColl2PositionInMu2e = mstmShutterSegmentPositionInMu2e + G4ThreeVector(0.0,0.0,mstmShutterSegmentLastHalfLength) + G4ThreeVector(0.0,0.0, mstmColl2UpStrSpace + mstmColl2HalfLength);
        
    boxColl2WithWindow.solid = new G4SubtractionSolid(boxColl2WithWindow.name,boxColl2,windowColl2,0,mstmColl2PositionInMu2e);

    finishNesting(boxColl2WithWindow,
                  mstmColl2Material,
                  0,
                  mstmColl2PositionInMu2e,
                  parent.logical,
                  0,
                  _config.getBool("mstm.visible"),
                  G4Colour::Magenta(),
                  _config.getBool("mstm.solid"),
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
    
    //TODO: Throw if pipe is too big to fit through the holes in the wall and magnet
    
    const TubsParams mstmPipe1Params(0., mstmPipe1ROut, mstmPipe1HalfLength);
    const TubsParams mstmPipe1GasParams(0.,  mstmPipe1RIn, mstmPipe1HalfLength - 2.0*mstmPipe1UpStrWindowHalfLength - 2.0*mstmPipe1DnStrWindowHalfLength);
    const TubsParams mstmPipe1UpStrWindowParams(0., mstmPipe1RIn, mstmPipe1UpStrWindowHalfLength);
    const TubsParams mstmPipe1DnStrWindowParams(0., mstmPipe1RIn, mstmPipe1DnStrWindowHalfLength);
    
    // Leave a gap between this and the previous element
    G4ThreeVector mstmPipe1PositionInMu2e    = mstmColl2PositionInMu2e + G4ThreeVector(0.0,0.0,mstmColl2HalfLength + mstmPipe1UpStrSpace + mstmPipe1HalfLength);
    
    VolumeInfo mstmPipe1Info = nestTubs( "mstmPipe1",
                                         mstmPipe1Params,
                                         mstmPipe1Material,
                                         0x0,
                                         mstmPipe1PositionInMu2e,// - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Red(),
                                         mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Red(),
                                            mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Red(),
                                            mstmSolid,
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
    // G4SubtractionSolid to alow for another volume placement
    // through it

    // Make the box for the collimator
    G4Box* boxColl3 = new G4Box("boxColl3",mstmColl3HalfWidth,mstmColl3HalfHeight,mstmColl3HalfLength);

    //Make the tube for the hole
    const TubsParams windColl3Params(0.0,                        //inner radius
                                     mstmColl3HoleROut,   //outer raius
                                     mstmColl3HalfLength  //
                                    );

    G4Tubs* windowColl3 = new G4Tubs( "window", 
                                   windColl3Params.data()[0], 
                                   windColl3Params.data()[1], 
                                   windColl3Params.data()[2]+2.,// to satisfy a G4SubtractionSolid feature
                                   windColl3Params.data()[3], 
                                   windColl3Params.data()[4]);

    // Combine into the Wall with the Hole
    VolumeInfo boxColl3WithWindow;
    boxColl3WithWindow.name = "boxColl3WithWindow";
          
    // We need to put the window on the z axis
    // Leave a gap between the this and the next element
    G4ThreeVector mstmColl3PositionInMu2e = mstmPipe1PositionInMu2e + G4ThreeVector(0.0,0.0,mstmPipe1HalfLength + mstmColl3UpStrSpace + mstmColl3HalfLength);
        
    boxColl3WithWindow.solid = new G4SubtractionSolid(boxColl3WithWindow.name,boxColl3,windowColl3,0,mstmColl3PositionInMu2e);

    finishNesting(boxColl3WithWindow,
                  mstmColl3Material,
                  0,
                  mstmColl3PositionInMu2e, /*  should is be mstmColl3PositionInMu2e ??? */
                  parent.logical,
                  0,
                  _config.getBool("mstm.visible"),
                  G4Colour::Magenta(),
                  _config.getBool("mstm.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);


    //----- Can (What is this?) --------------------------------------------------
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
    G4ThreeVector mstmCanPositionInMu2e    = mstmColl3PositionInMu2e + G4ThreeVector(0.0,0.0,mstmColl3HalfLength + mstmCanUpStrSpace + mstmCanHalfLength);
    
    VolumeInfo mstmCanInfo = nestTubs( "mstmCan",
                                         mstmCanParams,
                                         mstmCanMaterial,
                                         0x0,
                                         mstmCanPositionInMu2e,// - _hallOriginInMu2e,
                                         parent,
                                         0,
                                         mstmVisible,
                                         G4Color::Green(),
                                         mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Yellow(),
                                            mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Green(),
                                            mstmSolid,
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
                                            mstmVisible,
                                            G4Color::Green(),
                                            mstmSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    
    //----- the sensitive detector ("crystal") -------------------------------------------

    G4Material*  mstmCrystalMaterial   = materialFinder.get("mstm.crystal.material");
    const double mstmCrystalRIn        =  _config.getDouble("mstm.crystal.rIn");
    const double mstmCrystalROut       =  _config.getDouble("mstm.crystal.rOut");
    const double mstmCrystalHalfLength =  _config.getDouble("mstm.crystal.halfLength");

    //TODO: Throw if crystal doesn't fit inside the can

    
    const TubsParams mstmCrystalParams(mstmCrystalRIn, mstmCrystalROut, mstmCrystalHalfLength);

    G4ThreeVector mstmCrystalPositionInMu2e = mstmCanPositionInMu2e;

    VolumeInfo mstmCrystal = nestTubs("mstmCrystal",
                                      mstmCrystalParams,
                                      mstmCrystalMaterial,
                                      0x0,
                                      mstmCrystalPositionInMu2e,// - _hallOriginInMu2e,
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
    G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::MSTMCrystal());

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

      std::cout << __func__ << " mstmPipe1PositionInMu2e " << 
        mstmPipe1PositionInMu2e <<  std::endl;
      std::cout << __func__ << " mstmPipe1PositionInHall " << 
        mstmPipe1PositionInMu2e - _hallOriginInMu2e <<  std::endl;

      std::cout << __func__ << " mstmCanPositionInMu2e " << 
        mstmCanPositionInMu2e <<  std::endl;
      std::cout << __func__ << " mstmCanPositionInHall " << 
        mstmCanPositionInMu2e - _hallOriginInMu2e <<  std::endl;
    }

  }

} // end of constructMSTM;
