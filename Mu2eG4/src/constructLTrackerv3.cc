//
// Free function to construct version 3 of the LTracker
//
// $Id: constructLTrackerv3.cc,v 1.4 2010/06/17 20:31:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/17 20:31:13 $
//
// Original author Rob Kutschke
//
// Notes
//
// 1) Version 3 builds the LTracker by making physical mother volumes
//    for each sector.
// 2) There is bug in this version.  The physical volumes that bound
//    each sector of the octagon should be trapezoids, not boxes.
//

// C++ includes
#include <iostream>
#include <string>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructLTracker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "Mu2eG4/inc/StrawPlacer.hh"
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

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    double rOut  = CLHEP::mm * ltracker->rOut();
    double zHalf = CLHEP::mm * ltracker->zHalfLength();
    double z0    = CLHEP::mm * ltracker->z0();

    VolumeInfo trackerInfo;

    // Make the mother volume for the LTracker.
    string trackerName("LTrackerMother");
    G4Material* fillMaterial = findMaterialOrThrow(ltracker->fillMaterial());
    G4ThreeVector trackerOffset(0.,0.,z0-zOff);

    /*
      cout << "Tracker Offset: z0, zOff, z0-zOff: " 
      << z0 << " "
      << zOff << " "
      << z0-zOff << " "
      << endl;
    */

    trackerInfo.solid  = new G4Tubs( trackerName,
                                     0., rOut, zHalf, 0., 2.*M_PI );
    
    trackerInfo.logical = new G4LogicalVolume( trackerInfo.solid, fillMaterial, trackerName); 
    
    trackerInfo.physical =  new G4PVPlacement( 0, 
                                               trackerOffset, 
                                               trackerInfo.logical, 
                                               trackerName, 
                                               mother, 
                                               0, 
                                               0);

    // Visualization attributes of the the mother volume.
    {
      G4VisAttributes* visAtt = new G4VisAttributes(true, G4Colour::Green() );
      visAtt->SetForceSolid(true);
      // visAtt->SetForceSolid(false);
      visAtt->SetForceAuxEdgeVisible (false);
      visAtt->SetVisibility(true);
      trackerInfo.logical->SetVisAttributes(visAtt);
    }

    Straw const& straw        = ltracker->getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

    // Build logical volume for a straw.
    string strawName("Straw");
    VolumeInfo strawInfo;
      
    G4Material* strawMaterial = findMaterialOrThrow( detail.materialName(1) );
    strawInfo.solid  = new G4Tubs(strawName
                                  ,0.
                                  ,detail.outerRadius() * CLHEP::mm
                                  ,detail.halfLength()  * CLHEP::mm
                                  ,0.
                                  ,CLHEP::twopi*CLHEP::radian
                                  );
    
    strawInfo.logical = new G4LogicalVolume( strawInfo.solid
                                             , strawMaterial
                                             , strawName
                                             );

    vector<VolumeInfo> vinfo;

    for ( std::size_t idev = 0; idev<ltracker->getDevices().size(); ++idev){
      Device const& device = ltracker->getDevice(idev);

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

        G4RotationMatrix* rotft = new G4RotationMatrix(RForTrapezoids.inverse());

	G4RotationMatrix* rotfb = 0;

        CLHEP::HepRotationX RX(-sector.boxRxAngle());
        CLHEP::HepRotationY RY(-sector.boxRyAngle());
        CLHEP::HepRotationZ RZ(-sector.boxRzAngle());
        G4RotationMatrix* rot   = new G4RotationMatrix( RY*RX*RZ);
        G4RotationMatrix* rottr = new G4RotationMatrix( RForTrapezoids*RY*RX*RZ);

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
                   G4Colour::Magenta())
          :
          nestTrp( name,
                   sector.boxHalfLengths(),
                   fillMaterial,
                   rottr,
                   sector.boxOffset(),
                   trackerInfo.logical,
                   0,
                   G4Colour::Yellow());


//      the code below changes the recorded hits
// 	if (tmp.physical->CheckOverlaps(1000,0.0,false)){
// 	  std::cout << "Overlap while placing " << name << std::endl;
// 	}

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

          for ( std::size_t istr =0; istr<layer.nStraws(); ++istr){
            Straw const& straw = layer.getStraw(istr);

            // Position within the sector box.
            CLHEP::Hep3Vector position = origin + istr*delta;
            CLHEP::Hep3Vector trapezoidposition = trapezoidorigin + istr*trapezoiddelta;

            // Name of this physical volume.
            string sname = straw.name( "LTrackerStraw_");

            G4VPhysicalVolume* phys = (device.Id() ==  LTracker::vane) ?
              // G4VPhysicalVolume* phys = ( true ) ?
              new G4PVPlacement( rotfb,
                                 position,
                                 strawInfo.logical,
                                 sname, 
                                 sectorBoxInfo.logical, 
                                 0, 
                                 straw.Index().asInt()
                                 ) 
              :
              new G4PVPlacement( rotft,
                                 trapezoidposition,
                                 strawInfo.logical,
                                 sname, 
                                 sectorBoxInfo.logical, 
                                 0, 
                                 straw.Index().asInt()
                                 );
            
          } // loop over straws
        }   // loop over layers
        
      } // loop over sectors
    }   // loop over devices

    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4String strawSDname = "StrawGasVolume";
    StrawSD* strawSD     = new StrawSD( strawSDname, config );
    SDman->AddNewDetector( strawSD );
    strawInfo.logical->SetSensitiveDetector( strawSD );


    {
      //This leaks strawVisAtt.  
      //G4VisAttributes* strawVisAtt = new G4VisAttributes(true, G4Colour::Green() );
      //strawVisAtt->SetForceSolid(true);
      //strawVisAtt->SetForceAuxEdgeVisible (false);
      //strawVisAtt->SetVisibility(true);
      //strawInfo.logical->SetVisAttributes(strawVisAtt);

      // Do not draw straws since it is too time consuming.
      strawInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);

    }
    
    return trackerInfo;

  }

} // end namespace mu2e
