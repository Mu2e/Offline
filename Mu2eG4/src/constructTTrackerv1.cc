//
// Free function to construct version 1 of the TTracker
//
// $Id: constructTTrackerv1.cc,v 1.5 2010/07/07 16:44:14 genser Exp $
// $Author: genser $
// $Date: 2010/07/07 16:44:14 $
//
// Original author Rob Kutschke
//
// Notes
// 1) This version makes mother volumes per device and places
//    straws within that volume.  There is no per sector or per
//    manifold substructure.
// 2) Need to figure out how to manage the ownership of the
//    sensitive detector objects.
//

// C++ includes
#include <iostream>
#include <string>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructTTracker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"
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

  VolumeInfo constructTTrackerv1( G4LogicalVolume* mother, 
                                  double zOff,
                                  SimpleConfig const& config ){
    
    // Control of graphics for debugging the geometry.
    // Only instantiate sectors to be drawn.
    int devDraw = config.getInt("ttracker.devDraw",-1);
    int secDraw = config.getInt("ttracker.secDraw",-1);

    // Master geometry for the TTracker.
    GeomHandle<TTracker> ttracker;

    // Make the envelope volume that holds the full tracker.
    TubsParams envelopeParams = ttracker->getTrackerEnvelopeParams();

    G4ThreeVector trackerOffset( 0., 0., ttracker->z0()-zOff );

    G4Material* envelopeMaterial = findMaterialOrThrow(ttracker->envelopeMaterial());

    VolumeInfo info = nestTubs( "TrackerMother",
                                envelopeParams,
                                envelopeMaterial,
                                0,
                                trackerOffset,
                                mother,
                                0,
                                G4Color::Blue(),
                                config.getBool("ttracker.envelopeSolid",true)
                                );

    if (!config.getBool("ttracker.envelopeVisible",false)) {
      info.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }


    // Define the straws to be sensitive detectors.  Does this leak the StrawSD?
    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4String strawSDname = "StrawGasVolume";
    StrawSD* strawSD     = new StrawSD( strawSDname, config );
    SDman->AddNewDetector( strawSD );

    TubsParams deviceEnvelopeParams = ttracker->getDeviceEnvelopeParams();

    bool ttrackerDeviceEnvelopeVisible = config.getBool("ttracker.deviceEnvelopeVisible",true);
    bool ttrackerDeviceEnvelopeSolid = config.getBool("ttracker.deviceEnvelopeSolid",true);
    bool ttrackerStrawVisible          = config.getBool("ttracker.strawVisible",false);
    bool ttrackerStrawSolid          = config.getBool("ttracker.strawSolid",true);

    for ( size_t idev=0; idev<ttracker->nDevices(); ++idev ){

      if ( idev != devDraw && devDraw > -1 ) continue;

      const Device& device = ttracker->getDevice(idev);

      VolumeInfo devInfo  = nestTubs("TTrackerDeviceEnvelope",
                                     deviceEnvelopeParams,
                                     envelopeMaterial,
                                     0,
                                     device.origin(),
                                     info.logical,
                                     idev,
                                     G4Color::Magenta(),
                                     ttrackerDeviceEnvelopeSolid
                                     );

      if (!ttrackerDeviceEnvelopeVisible) {
        devInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
      }

      for ( size_t isec = 0; isec<device.nSectors(); ++isec){
        if ( isec != secDraw && secDraw > -1 ) continue;

        const Sector& sector = device.getSector(isec);

        for ( size_t ilay =0; ilay<sector.nLayers(); ++ilay ){
          const Layer& layer = sector.getLayer(ilay);
          
          for ( int istr=0; istr<layer.nStraws(); ++istr ){
            const Straw& straw = layer.getStraw(istr);
            StrawDetail const& detail = straw.getDetail();

            TubsParams strawParams( 0., detail.outerRadius(), detail.halfLength() );
            G4ThreeVector mid    = straw.getMidPoint() - device.origin();
            G4Material* strawGas = findMaterialOrThrow(detail.gasMaterialName());

            const G4ThreeVector& w = straw.direction();

            double theta = w.theta()*CLHEP::radian;
            double phi   = w.phi()*CLHEP::radian;
            double alpha = M_PI/2.-phi;
          
            // This choice of Euler angles ensures that the z axis fixed to the straw
            // object points in the direction of w, as measured in the world frame
            G4RotationMatrix* rot = new G4RotationMatrix( -alpha, -theta, alpha );

            VolumeInfo strawInfo  = nestTubs("TTrackerStraw",
                                             strawParams,
                                             strawGas,
                                             rot,
                                             mid,
                                             devInfo.logical,
                                             straw.Index().asInt(),
                                             G4Color::Green(),
                                             ttrackerStrawSolid
                                             );

            // Make this straw a sensitive detector.
            strawInfo.logical->SetSensitiveDetector( strawSD );
            if (!ttrackerStrawVisible) {
              strawInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
            }

          }   // end loop over straws
        }     // end loop over layers
      }       // end loop over sectors
    }         // end loop over devices

    return info;
  }

} // end namespace mu2e
