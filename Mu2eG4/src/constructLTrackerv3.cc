//
// Free function to construct version 3 of the LTracker
//
// $Id: constructLTrackerv3.cc,v 1.14 2010/12/02 17:49:10 genser Exp $
// $Author: genser $
// $Date: 2010/12/02 17:49:10 $
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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructLTracker.hh"
#include "Mu2eG4/inc/G4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTrp.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e{

  VolumeInfo constructLTrackerv3( G4LogicalVolume* mother, 
                                  double zOff,
                                  SimpleConfig const& config ){

    G4Helper    & _helper = *(edm::Service<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    int verbosityLevel = config.getInt("ltracker.verbosityLevel",0);

    // Master geometry for the LTracker.
    // GeomHandle<LTracker> ltracker;
    LTracker const & ltracker = *(GeomHandle<LTracker>());

    // Ugly hack; 
    // add padding to ensure that this encloses the rotated volumes
    // This should be done inside the LTracker object.
    double rOut  = CLHEP::mm * ltracker.rOut()+10.;
    double zHalf = CLHEP::mm * ltracker.zHalfLength()+30.;
    double z0    = CLHEP::mm * ltracker.z0();

    bool trackerSolid        = config.getBool("ltracker.solid",true);
    bool trackerVisible      = config.getBool("ltracker.visible",true);
    
    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV = true;

    VolumeInfo trackerInfo;

    // Make the mother volume for the LTracker.
    trackerInfo.name = "TrackerMother";
    G4Material* fillMaterial = findMaterialOrThrow(ltracker.fillMaterial());
    G4ThreeVector trackerOffset(0.,0.,z0-zOff);

    verbosityLevel > 0 && 
      cout << "Tracker Offset: z0, zOff, z0-zOff: " 
           << z0 << " "
           << zOff << " "
           << z0-zOff << " "
           << endl;

    trackerInfo.solid  = new G4Tubs( trackerInfo.name,
                                     0., rOut, zHalf, 0., 2.*M_PI );
    
    trackerInfo.logical = new G4LogicalVolume( trackerInfo.solid, fillMaterial, trackerInfo.name); 
    
    trackerInfo.physical =  new G4PVPlacement( 0, 
                                               trackerOffset, 
                                               trackerInfo.logical, 
                                               trackerInfo.name, 
                                               mother, 
                                               0, 
                                               0,
                                               doSurfaceCheck);

    _helper.addVolInfo(trackerInfo);


    // Visualization attributes of the the mother volume.

    if (!trackerVisible) {
      trackerInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Green() ));
      visAtt->SetForceSolid(trackerSolid);
      visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      trackerInfo.logical->SetVisAttributes(visAtt);
    }

    Straw const& straw        = ltracker.getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

    G4ThreeVector const zeroVector(0.0,0.0,0.0);

    // Build logical volumes related to the straw: wall, gas, wire

    // Build logical volume for straw wall.
    string strawWallName("StrawWall");
    VolumeInfo strawWall;
      
    strawWall.solid  = new G4Tubs(strawWallName
                                  ,0.0
                                  ,detail.outerRadius() * CLHEP::mm
                                  ,detail.halfLength() * CLHEP::mm
                                  ,0.
                                  ,CLHEP::twopi*CLHEP::radian
                                  );
    
    strawWall.logical = new G4LogicalVolume( strawWall.solid
                                             ,findMaterialOrThrow(detail.wallMaterialName())
                                             ,strawWallName
                                             );

    // Build logical volume for straw gas.
    string strawGasName("StrawGas");
    VolumeInfo strawGas;
      
    strawGas.solid  = new G4Tubs(strawGasName
                                 ,0.0
                                 ,detail.innerRadius() * CLHEP::mm
                                 ,detail.halfLength() * CLHEP::mm
                                 ,0.
                                 ,CLHEP::twopi*CLHEP::radian
                                 );
    
    strawGas.logical = new G4LogicalVolume( strawGas.solid
                                            ,findMaterialOrThrow(detail.gasMaterialName())
                                            ,strawGasName
                                            );

    strawGas.physical =  new G4PVPlacement( 0, zeroVector, 
                                            strawGas.logical, strawGasName,
                                            strawWall.logical, false, 0, 
                                            doSurfaceCheck);

    // Build logical volume for straw wire.
    string strawWireName("StrawWire");
    VolumeInfo strawWire;
      
    strawWire.solid  = new G4Tubs(strawWireName
                                  ,0.0
                                  ,detail.wireRadius() * CLHEP::mm
                                  ,detail.halfLength() * CLHEP::mm
                                  ,0.
                                  ,CLHEP::twopi*CLHEP::radian
                                  );
    
    strawWire.logical = new G4LogicalVolume( strawWire.solid
                                             ,findMaterialOrThrow(detail.wireMaterialName())
                                             ,strawWireName
                                             );

    strawWire.physical =  new G4PVPlacement( 0, zeroVector, 
                                             strawWire.logical, strawWireName , 
                                             strawGas.logical, false, 0, 
                                             doSurfaceCheck);

    if (config.getBool("ltracker.strawVisible",false)) {
      
      bool strawSolid = config.getBool("ltracker.strawSolid",true);

      G4VisAttributes* strawWallVisAtt = reg.add(G4VisAttributes(true, G4Colour::Yellow() ));
      strawWallVisAtt->SetForceSolid(strawSolid);
      strawWallVisAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      strawWall.logical->SetVisAttributes(strawWallVisAtt);

      G4VisAttributes* strawGasVisAtt = reg.add(G4VisAttributes(true, G4Colour::Green() ));
      strawGasVisAtt->SetForceSolid(strawSolid);
      strawGasVisAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      strawGas.logical->SetVisAttributes(strawGasVisAtt);

      G4VisAttributes* strawWireVisAtt = reg.add(G4VisAttributes(true, G4Colour::Red() ));
      strawWireVisAtt->SetForceSolid(strawSolid);
      strawWireVisAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      strawWire.logical->SetVisAttributes(strawWireVisAtt);

    }
    
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

    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4String strawSDname = "StrawGasVolume";
    StrawSD* strawSD     = new StrawSD( strawSDname, config );
    SDman->AddNewDetector( strawSD );
    strawGas.logical->SetSensitiveDetector( strawSD );

    return trackerInfo;

  }

} // end namespace mu2e
