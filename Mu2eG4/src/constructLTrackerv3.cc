//
// Free function to construct version 3 of the LTracker
//
// $Id: constructLTrackerv3.cc,v 1.18 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Rob Kutschke
//
// Notes
//
// 1) Version 3 builds the LTracker by making physical mother volumes
//    for each sector.
//

// C++ includes
#include <iostream>
#include <string>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructLTracker.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTrp.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e{

  VolumeInfo constructLTrackerv3( G4LogicalVolume* mother, 
                                  double zOff,
                                  SimpleConfig const& config ){

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    int verbosityLevel = config.get<int>("ltracker.verbosityLevel",0);

    // Master geometry for the LTracker.
    // GeomHandle<LTracker> ltracker;
    LTracker const & ltracker = *(GeomHandle<LTracker>());

    // Ugly hack; 
    // add padding to ensure that this encloses the rotated volumes
    // This should be done inside the LTracker object.
    double rOut  = CLHEP::mm * ltracker.rOut()+10.;
    double zHalf = CLHEP::mm * ltracker.zHalfLength()+30.;
    double z0    = CLHEP::mm * ltracker.z0();

    bool const trackerSolid        = config.get<bool>("ltracker.solid",true);
    bool const trackerVisible      = config.get<bool>("ltracker.visible",true);
    
    bool const strawVisible        = config.get<bool>("ltracker.strawVisible",false); 
    bool const strawSolid          = config.get<bool>("ltracker.strawSolid",true);


    bool const forceAuxEdgeVisible = config.get<bool>("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = config.get<bool>("g4.doSurfaceCheck",false);
    bool const placePV = true;

    // Make the mother volume for the LTracker.

    G4Material* fillMaterial = findMaterialOrThrow(ltracker.fillMaterial());
    G4ThreeVector trackerOffset(0.,0.,z0-zOff);

    verbosityLevel > 0 && 
      cout << "Tracker Offset: z0, zOff, z0-zOff: " 
           << z0 << " "
           << zOff << " "
           << z0-zOff << " "
           << endl;


    TubsParams envelopeParams( 0., rOut, zHalf);

    VolumeInfo trackerInfo = nestTubs( "TrackerMother",
                                       envelopeParams,
                                       fillMaterial,
                                       0,
                                       trackerOffset,
                                       mother,
                                       0,
                                       trackerVisible,
                                       G4Colour::Green(),
                                       trackerSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    Straw const& straw        = ltracker.getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

    G4ThreeVector const zeroVector(0.0,0.0,0.0);

    // Build logical volumes related to the straw: wall, gas, wire

    // Build logical volume for straw wall with Gas & Wire inside it

    TubsParams strawWallTubeParams(0., detail.outerRadius() * CLHEP::mm, detail.halfLength() * CLHEP::mm);

    VolumeInfo strawWall = nestTubs( "StrawWall",
                                     strawWallTubeParams,
                                     findMaterialOrThrow(detail.wallMaterialName()),
                                     0,
                                     zeroVector,
                                     0,
                                     0,
                                     strawVisible,
                                     G4Colour::Yellow(),
                                     strawSolid,
                                     forceAuxEdgeVisible,
                                     false,
                                     doSurfaceCheck
                                     );
                                     
    // Place straw gas inside StrawWall

    TubsParams strawGasTubeParams(0., detail.innerRadius() * CLHEP::mm, detail.halfLength() * CLHEP::mm);

    VolumeInfo strawGas = nestTubs( "StrawGas", 
                                    strawGasTubeParams,
                                    findMaterialOrThrow(detail.gasMaterialName()),
                                    0,
                                    zeroVector,
                                    strawWall.logical,
                                    0,
                                    strawVisible,
                                    G4Colour::Green(),
                                    strawSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    // Place straw wire inside StrawGas

    TubsParams strawWireTubeParams(0., detail.wireRadius() * CLHEP::mm, detail.halfLength() * CLHEP::mm);

    VolumeInfo strawWire = nestTubs( "StrawWire",
                                     strawWireTubeParams,
                                     findMaterialOrThrow(detail.wireMaterialName()),
                                     0,
                                     zeroVector,
                                     strawGas.logical,
                                     0,
                                     strawVisible,
                                     G4Colour::Red(),
                                     strawSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    vector<VolumeInfo> vinfo;

    // device loop
    for ( std::size_t idev = 0; idev<ltracker.getDevices().size(); ++idev){
      Device const& device = ltracker.getDevice(idev);

      for ( std::size_t isec =0; isec<device.getSectors().size(); ++isec){
        Sector const& sector = device.getSector(isec);

        // Name of this sector as string.
        string name = sector.name("LTrackerSector_");

        // Construct the rotation.  
        // This rotation is the inverse of the one in v2.
        // Note the sign and the reversed order : active/passive  confusion.
        // Need to understand if this causes memory leak.

        double const trapezoidRxAngle = M_PI_2;
        CLHEP::HepRotationX RForTrapezoids(trapezoidRxAngle);

        G4RotationMatrix* rotft = reg.add(G4RotationMatrix(RForTrapezoids.inverse()));

        G4RotationMatrix* rotfb = 0;

        CLHEP::HepRotationX RX(-sector.boxRxAngle());
        CLHEP::HepRotationY RY(-sector.boxRyAngle());
        CLHEP::HepRotationZ RZ(-sector.boxRzAngle());
        G4RotationMatrix* rot   = reg.add(G4RotationMatrix( RY*RX*RZ));
        G4RotationMatrix* rottr = reg.add(G4RotationMatrix( RForTrapezoids*RY*RX*RZ));

        // Make a physical volume for this sector.  Same material as the 
        // main LTracker volume ( some sort of vacuum ).

        VolumeInfo tmp = ( device.Id() ==  LTracker::vane) ?
          //VolumeInfo tmp = ( true ) ?
          nestBox( name,
                   sector.boxHalfLengths(),
                   fillMaterial,
                   rot,
                   sector.boxOffset(),
                   trackerInfo.logical,
                   0,
                   trackerVisible,
                   G4Colour::Magenta(),
                   trackerSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck)
          :
          nestTrp( name,
                   sector.boxHalfLengths(),
                   fillMaterial,
                   rottr,
                   sector.boxOffset(),
                   trackerInfo.logical,
                   0,
                   trackerVisible,
                   G4Colour::Yellow(),
                   trackerSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck);

        // the CheckOverlaps code below changes the recorded hits
        // if (tmp.physical->CheckOverlaps(1000,0.0,false)){
        //     std::cout << "Overlap while placing " << name << std::endl;
        //  }

        vinfo.push_back(tmp);
        VolumeInfo const& sectorBoxInfo = vinfo.back();

        CLHEP::Hep3Vector const& delta  = sector.getBaseDelta();

        // for trapezoidal wedges the coordinate system is different
        // due to the Geant4 conventions for Trapezoids
        CLHEP::Hep3Vector const trapezoiddelta = CLHEP::Hep3Vector(delta.x(),delta.z(),delta.y());

        for ( std::size_t ilay =0; ilay<sector.getLayers().size(); ++ilay){
          Layer const& layer = sector.getLayer(ilay);

          CLHEP::Hep3Vector const& origin = sector.getBasePosition().at(ilay);
          CLHEP::Hep3Vector const trapezoidorigin = CLHEP::Hep3Vector(origin.x(),origin.z(),origin.y());

          for ( int istr =0; istr<layer.nStraws(); ++istr){
            Straw const& straw = layer.getStraw(istr);

            // Position all straw volumes within the sector box.
            CLHEP::Hep3Vector position = origin + istr*delta;
            CLHEP::Hep3Vector trapezoidposition = trapezoidorigin + istr*trapezoiddelta;

            // Name of this physical volume
            // straw.name( "LTrackerStrawWall_");
            // Formatted string embedding the id of the straw.
            // std::string name( std::string const& base ) const;

            // G4 takes ownership of whichever object is created.

            // note that strawWall.logical has Gas & Wire inside it

            if (device.Id() ==  LTracker::vane) {
              new G4PVPlacement( rotfb,
                                 position,
                                 strawWall.logical,
                                 straw.name( "LTrackerStrawWall_"),
                                 sectorBoxInfo.logical, 
                                 0, 
                                 straw.Index().asInt(),
                                 doSurfaceCheck);
            }else{
              new G4PVPlacement( rotft,
                                 trapezoidposition,
                                 strawWall.logical,
                                 straw.name( "LTrackerStrawWall_"),
                                 sectorBoxInfo.logical, 
                                 0, 
                                 straw.Index().asInt(),
                                 doSurfaceCheck);
            }
            
          } // loop over straws
        }   // loop over layers
        
      } // loop over sectors
    }   // loop over devices

    strawGas.logical->
      SetSensitiveDetector(G4SDManager::GetSDMpointer()->
                           FindSensitiveDetector(SensitiveDetectorName::StrawGasVolume()));

    return trackerInfo;

  }

} // end namespace mu2e
