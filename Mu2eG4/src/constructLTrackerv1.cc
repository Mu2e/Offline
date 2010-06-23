//
// Free function to construct version 1 of the LTracker
//
// $Id: constructLTrackerv1.cc,v 1.3 2010/06/23 23:29:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/23 23:29:16 $
//
// Original author Rob Kutschke
//
// Notes
// 1) 
//   Version 1 places each straw inside the tracker mother volume.
//   There is no substructure.
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
#include "Mu2eG4/inc/nestTubs.hh"

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

  VolumeInfo constructLTrackerv1( G4LogicalVolume* mother, 
                                  double zOff,
                                  SimpleConfig const& config ){

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    double rOut  = CLHEP::mm * ltracker->rOut();
    double zHalf = CLHEP::mm * ltracker->zHalfLength();
    double z0    = CLHEP::mm * ltracker->z0();

    VolumeInfo trackerInfo;

    // Make the mother volume for the LTracker.
    string trackerName("TrackerMother");
    G4Material* fillMaterial = findMaterialOrThrow(ltracker->fillMaterial());
    G4ThreeVector trackerOffset(0.,0.,z0-zOff);

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
      visAtt->SetForceAuxEdgeVisible (false);
      visAtt->SetVisibility(true);
      trackerInfo.logical->SetVisAttributes(visAtt);
    }

    // For now cheat and assume that all straws are the same and are just made of gas with
    // no walls or wires.
    Straw const& straw = ltracker->getStraw( StrawId( LTracker::wedge, 0, 0, 0) );
    StrawDetail const& detail = straw.getDetail();

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

    // Define the straws to be sensitive detectors.
    // Does this leak the SDman?
    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4String strawSDname = "StrawGasVolume";
    StrawSD* strawSD     = new StrawSD( strawSDname, config );
    SDman->AddNewDetector( strawSD );
    strawInfo.logical->SetSensitiveDetector( strawSD );


    // Does this leak strawVisAtt ??
    {
      G4VisAttributes* strawVisAtt = new G4VisAttributes(true, G4Colour::Green() );
      strawVisAtt->SetForceSolid(false);
      strawVisAtt->SetForceAuxEdgeVisible (false);
      strawVisAtt->SetVisibility(false);
      strawInfo.logical->SetVisAttributes(strawVisAtt);
    }

    // Cheat again and do not bother to segment the tracker volume.
    // Just place the straws in the final position.
    // Need to properly segment the volume on a per sector basis in
    // order to improve G4 statup speed.
    StrawPlacer placer( "StrawPhys", strawInfo.logical, trackerInfo.logical );
    ltracker->forAllStraws( placer);

    return trackerInfo;

  }

} // end namespace mu2e
