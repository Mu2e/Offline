//
// Free function to construct version 2 of the LTracker
//
// $Id: constructLTrackerv2.cc,v 1.3 2010/06/23 23:29:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/06/23 23:29:16 $
//
// Original author Rob Kutschke
//
// Notes
//
// 1) Version 2 builds the LTracker using assembly volumes.
//      a) Make one assembly volume for the octant sides and one for the vanes.
//      b) Use imprint to make 8 copies of each and put each in the right location.
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
#include "G4AssemblyVolume.hh"
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

  VolumeInfo constructLTrackerv2( G4LogicalVolume* mother, 
                                  double zOff,
                                  SimpleConfig const& config ){

    // Cannot make std::vector of auto_ptr. 
    // So use statically dimensioned array instead.  Remember to check dimensions.
    static int const ndevices = 2;
    static std::auto_ptr<G4AssemblyVolume>    _lTrackerAssemblyVols[ndevices];


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

    if ( ltracker->getDevices().size() != ndevices ){
      throw cms::Exception("GEOM")
        << "Unexpected number of devices in the LTracker: "
        << ltracker->getDevices().size()
        << "\n";
    }

    // Build an assembly volume for each of the octants and vanes.
    for ( std::size_t idev=0; idev<ltracker->getDevices().size(); ++idev ){

      // Assume all sectors are the same as sector 0.
      Sector const& sec = ltracker->getSector(SectorId(idev,0));
        
      _lTrackerAssemblyVols[idev] = auto_ptr<G4AssemblyVolume> (new G4AssemblyVolume());

      for ( std::size_t ilay =0; ilay<sec.getLayers().size(); ++ilay){
        Layer const& lay = sec.getLayer(ilay);
        CLHEP::Hep3Vector const& origin = sec.getBasePosition().at(ilay);
        CLHEP::Hep3Vector const& delta  = sec.getBaseDelta();

        StrawId id(idev,0,ilay,0);
        
        for ( std::size_t istr = 0; istr<lay.nStraws(); ++istr ){
          CLHEP::Hep3Vector position = origin + istr*delta;
          _lTrackerAssemblyVols[idev]->AddPlacedVolume( strawInfo.logical, position, 0); 
        }
      }
      
    }

    // Imprint the assembly volumes.
    for ( std::size_t idev = 0; idev<ltracker->getDevices().size(); ++idev){
      Device const& device = ltracker->getDevice(idev);

      for ( std::size_t isec =0; isec<device.getSectors().size(); ++isec){
        Sector const& sector = device.getSector(isec);

        CLHEP::HepRotationX RX(sector.boxRxAngle());
        CLHEP::HepRotationY RY(sector.boxRyAngle());
        CLHEP::HepRotationZ RZ(sector.boxRzAngle());
     
        // Need to understand if this causes memory leak.
        G4RotationMatrix* rot = new G4RotationMatrix( RZ*RX*RY);

        // MakeImprint requires non-const argument.
        G4ThreeVector offset = sector.boxOffset();


        // Copy numbers start at base+1
        StrawId id(idev,isec,0,0);
        int baseCopyNumber = ltracker->getStraw(id).Index().asInt()-1;

        _lTrackerAssemblyVols[idev]->MakeImprint( trackerInfo.logical, offset, rot, baseCopyNumber); 

      }
    }

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
    
    return trackerInfo;

  }

} // end namespace mu2e
