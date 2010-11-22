//
// Free function to construct version 3 of the TTracker
//
// $Id: constructTTrackerv3.cc,v 1.10 2010/11/22 05:21:22 genser Exp $
// $Author: genser $
// $Date: 2010/11/22 05:21:22 $
//
// Original author KLG based on RKK using different methodology
//
// Notes

//     This version makes logical mother volumes per device and per
//     sector and places sectors in device and straws in sector
//     It has only one sector/device logical volume placed several times
//     This versoin has a negligeable construction time and a much smaler memory footprint


// C++ includes
#include <iostream>
#include <iomanip>
#include <string>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructTTracker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4String.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e{

  VolumeInfo constructTTrackerv3( G4LogicalVolume* mother, 
                                  double zOff,
                                  SimpleConfig const& config ){
    
    // Control of graphics for debugging the geometry.
    // Only instantiate sectors to be drawn.
    int devDraw = config.getInt("ttracker.devDraw",-1);
    int secDraw = config.getInt("ttracker.secDraw",-1);
    G4bool doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false);


    // Master geometry for the TTracker.
    //    GeomHandle<TTracker> ttracker;
    TTracker const & ttracker = *(GeomHandle<TTracker>());

    // Make the envelope volume that holds the full tracker.
    TubsParams envelopeParams = ttracker.getTrackerEnvelopeParams();

    //     int const oldp = cout.precision();
    //     int const oldw = cout.width();
    //     cout << "Debugging tracker env envelopeParams ir,or,zhl,phi0,phimax:            " <<
    //       "   " << 
    //       fixed << setprecision(8) << setw(14) << envelopeParams.innerRadius << ", " <<
    //       fixed << setprecision(8) << setw(14) << envelopeParams.outerRadius << ", " <<
    //       fixed << setprecision(8) << setw(14) << envelopeParams.zHalfLength << ", " <<
    //       fixed << setprecision(8) << setw(14) << envelopeParams.phi0        << ", " <<
    //       fixed << setprecision(8) << setw(14) << envelopeParams.phiMax      << ", " <<
    //       endl;

    G4ThreeVector trackerOffset( 0., 0., ttracker.z0()-zOff );

    G4Material* envelopeMaterial = findMaterialOrThrow(ttracker.envelopeMaterial());

    VolumeInfo motherInfo = nestTubs( "TrackerMother",
                                      envelopeParams,
                                      envelopeMaterial,
                                      0,
                                      trackerOffset,
                                      mother,
                                      0,
                                      G4Color::Blue(),
                                      config.getBool("ttracker.envelopeSolid",true),
                                      doSurfaceCheck
                                      );

    if (!config.getBool("ttracker.envelopeVisible",false)) {
      motherInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
    } else {
      // leak?
      G4VisAttributes* visAtt = new G4VisAttributes(motherInfo.logical->GetVisAttributes());
      visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
      motherInfo.logical->SetVisAttributes(visAtt);
    }


    // Define the straws to be sensitive detectors.  Does this leak the StrawSD?
    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4String strawSDname = "StrawGasVolume";
    StrawSD* strawSD     = new StrawSD( strawSDname, config );
    SDman->AddNewDetector( strawSD );

    TubsParams deviceEnvelopeParams = ttracker.getDeviceEnvelopeParams();

    bool ttrackerDeviceEnvelopeVisible = config.getBool("ttracker.deviceEnvelopeVisible",false);
    bool ttrackerDeviceEnvelopeSolid = config.getBool("ttracker.deviceEnvelopeSolid",true);
    bool ttrackerSectorEnvelopeVisible = config.getBool("ttracker.sectorEnvelopeVisible",false);
    bool ttrackerSectorEnvelopeSolid = config.getBool("ttracker.sectorEnvelopeSolid",true);
    bool ttrackerStrawVisible          = config.getBool("ttracker.strawVisible",false);
    bool ttrackerStrawSolid          = config.getBool("ttracker.strawSolid",true);

    // construct one logical device (device # 0) with straws inside it

    //nestSomething create a physical volume, we need a logical one first

    size_t idev = 0;
    const Device& device = ttracker.getDevice(idev);

    VolumeInfo devInfo;

    G4String const devName = "TTrackerDeviceEnvelope";

    //     cout << "Debugging device env idev, deviceEnvelopeParams ir,or,zhl,phi0,phimax: " <<
    //       idev << ", " << 
    //       fixed << setprecision(8) << setw(14) << deviceEnvelopeParams.innerRadius << ", " <<
    //       fixed << setprecision(8) << setw(14) << deviceEnvelopeParams.outerRadius << ", " <<
    //       fixed << setprecision(8) << setw(14) << deviceEnvelopeParams.zHalfLength << ", " <<
    //       fixed << setprecision(8) << setw(14) << deviceEnvelopeParams.phi0        << ", " <<
    //       fixed << setprecision(8) << setw(14) << deviceEnvelopeParams.phiMax      << ", " <<
    //       endl;

    devInfo.solid    = new G4Tubs(devName, 
                                  deviceEnvelopeParams.innerRadius,
                                  deviceEnvelopeParams.outerRadius,
                                  deviceEnvelopeParams.zHalfLength, 
                                  deviceEnvelopeParams.phi0, 
                                  deviceEnvelopeParams.phiMax      
                                  );
    
    devInfo.logical  = new G4LogicalVolume( devInfo.solid, envelopeMaterial, devName);
    
    if (!ttrackerDeviceEnvelopeVisible) {
      devInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
    } 
    else {
      G4VisAttributes* visAtt = new G4VisAttributes(true, G4Color::Magenta());
      visAtt->SetForceSolid(ttrackerDeviceEnvelopeSolid);
      visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
      devInfo.logical->SetVisAttributes(visAtt);
    }
    // place straws etc... wrt the envelope

    // create a "sector" volume

    G4String const secName = "TTrackerSectorEnvelope";

    // Construct One sector logical volume (and then place it N times)

    const size_t isec = 0;

    const Sector& sector = device.getSector(isec);

    // constructing sector envelope

    // Make a logical volume for this sector.  

    // reuse device attributes for now

    // we need to "unfold" nestTrp;

    VolumeInfo secInfo;

    secInfo.solid   = new G4Trd ( secName,
                                  sector.boxHalfLengths().at(4),
                                  sector.boxHalfLengths().at(3),
                                  sector.boxHalfLengths().at(2),
                                  sector.boxHalfLengths().at(2),
                                  sector.boxHalfLengths().at(1)
                                  );

    //     cout << "Debugging sector box isec, sector.boxHalfLengths().at(4,3,2,2,1): " <<
    //       isec << ", " << 
    //       fixed << setprecision(8) << setw(14) << sector.boxHalfLengths().at(4) << ", " <<
    //       fixed << setprecision(8) << setw(14) << sector.boxHalfLengths().at(3) << ", " <<
    //       fixed << setprecision(8) << setw(14) << sector.boxHalfLengths().at(2) << ", " <<
    //       fixed << setprecision(8) << setw(14) << sector.boxHalfLengths().at(2) << ", " <<
    //       fixed << setprecision(8) << setw(14) << sector.boxHalfLengths().at(1) << ", " <<
    //       endl << setprecision(oldp) << setw(oldw);

    secInfo.logical  = new G4LogicalVolume( secInfo.solid, envelopeMaterial, secName); 
    
    if (!ttrackerSectorEnvelopeVisible) {
      secInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
    } 
    else {
      // it will be allways Cyan if done this way
      G4Color color = (isec%2 == 1) ? G4Color::Blue() : G4Color::Cyan();
      G4VisAttributes* visAtt = new G4VisAttributes(true, color);
      visAtt->SetForceSolid(ttrackerSectorEnvelopeSolid);
      visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
      secInfo.logical->SetVisAttributes(visAtt);
    }

    static double const tRAngle  = M_PI_2;
    static double const tRAngle2 = M_PI;
    CLHEP::HepRotationX RXForTrapezoids(tRAngle);
    CLHEP::HepRotationX RX2ForTrapezoids(tRAngle2);
    CLHEP::HepRotationY RYForTrapezoids(tRAngle);
    CLHEP::HepRotationZ RZForTrapezoids(tRAngle);
    
    //leak?
    G4RotationMatrix* rotTub = new G4RotationMatrix(RYForTrapezoids);

    for ( int ilay =0; ilay<sector.nLayers(); ++ilay ){

      //      cout << "Debugging constructTTrackerv3 ilay: " << ilay << endl;

      const Layer& layer = sector.getLayer(ilay);
          
      for ( int istr=0; istr<layer.nStraws(); ++istr ){

        // "second" layer will have fewer straws (for now) also see TTrackerMaker
        // no, it complicates StrawSD and TTrackerMaker 
        // if( ilay%2==1 && istr+1 == layer.nStraws() ) break;

        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& detail = straw.getDetail();

        TubsParams strawWallParams( 0.0,                  detail.outerRadius(), detail.halfLength() );
        TubsParams strawGasParams ( 0.0,                  detail.innerRadius(), detail.halfLength() );
        TubsParams strawWireParams( 0.0,                  detail.wireRadius(),  detail.halfLength() );

        // we are placing the straw w.r.t the trapezoid...
        // the trapezoid aka device envelope has a different coordinate system x->z, z->y, y->x

        G4ThreeVector const mid(straw.getMidPoint().y() - sector.boxOffset().y(),
                                straw.getMidPoint().z() - sector.boxOffset().z(),
                                straw.getMidPoint().x() - sector.boxOffset().x());

        G4ThreeVector const zeroVector(0.0,0.0,0.0);

        //         cout << "Debugging istr: " << istr << 
        //           " mid: " << mid << 
        //           ", straw.MidPoint " << straw.getMidPoint() << 
        //           ", sector.boxOffset " <<  sector.boxOffset() << 
        //           ", device.origin " << device.origin() <<
        //           endl;

        //         cout << "Debugging istr: " << istr << " mid: " << 
        //           mid << ", halflenght " << detail.halfLength() << endl;

        // look at StrawSD to see how the straw index is reconstructed

        //         cout << "Debugging straw.Id(), straw.Index() " << 
        //           straw.Id() << ", " << straw.Index() << endl;

        // make the straws more distinguishable when displayed
        G4Color wallColor = (ilay%2 == 1) ? 
          ((istr%2 == 0) ? G4Color::Green() : G4Color::Yellow()) :
          ((istr%2 == 0) ? G4Color::Red() : G4Color::Blue());

        G4Color gasColor = (ilay%2 == 0) ? 
          ((istr%2 == 0) ? G4Color::Green() : G4Color::Yellow()) :
          ((istr%2 == 0) ? G4Color::Red() : G4Color::Blue());

        G4Color wireColor = G4Color::Cyan();

        //         cout << "Debugging Straw istr, RYForTrapezoids, midpoint: " <<
        //           istr << ", " << RYForTrapezoids << ", " <<
        //           fixed << setprecision(8) << setw(14) << mid << ", " <<
        //           endl << setprecision(oldp) << setw(oldw);

        bool strawAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);

        VolumeInfo strawWallInfo  = nestTubs("TTrackerStrawWall",
                                             strawWallParams,
                                             findMaterialOrThrow(detail.wallMaterialName() ),
                                             rotTub,
                                             mid,
                                             secInfo.logical,
                                             straw.Index().asInt(),
                                             wallColor,
                                             ttrackerStrawSolid,
                                             doSurfaceCheck
                                             );

        if (!ttrackerStrawVisible) {
          strawWallInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          // leak?
          G4VisAttributes* visAtt = new G4VisAttributes(*strawWallInfo.logical->GetVisAttributes());
          visAtt->SetForceAuxEdgeVisible(strawAuxEdgeVisible);
          strawWallInfo.logical->SetVisAttributes(visAtt);
        }

        VolumeInfo strawGasInfo  = nestTubs("TTrackerStrawGas",
                                            strawGasParams,
                                            findMaterialOrThrow(detail.gasMaterialName()),
                                            0,
                                            zeroVector,
                                            strawWallInfo.logical,
                                            straw.Index().asInt(),
                                            gasColor,
                                            ttrackerStrawSolid,
                                            doSurfaceCheck
                                            );

        if (!ttrackerStrawVisible) {
          strawGasInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          // leak?
          G4VisAttributes* visAtt = new G4VisAttributes(*strawGasInfo.logical->GetVisAttributes());
          visAtt->SetForceAuxEdgeVisible(strawAuxEdgeVisible);
          strawGasInfo.logical->SetVisAttributes(visAtt);
        }

        VolumeInfo strawWireInfo  = nestTubs("TTrackerStrawWire",
                                             strawWireParams,
                                             findMaterialOrThrow(detail.wireMaterialName()),
                                             0,
                                             zeroVector,
                                             strawGasInfo.logical,
                                             straw.Index().asInt(),
                                             wireColor,
                                             ttrackerStrawSolid,
                                             doSurfaceCheck
                                             );

        if (!ttrackerStrawVisible) {
          strawWireInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          // leak?
          G4VisAttributes* visAtt = new G4VisAttributes(*strawWireInfo.logical->GetVisAttributes());
          visAtt->SetForceAuxEdgeVisible(strawAuxEdgeVisible);
          strawWireInfo.logical->SetVisAttributes(visAtt);
        }

        // Make gas of this straw a sensitive detector.
        strawGasInfo.logical->SetSensitiveDetector( strawSD );

      }   // end loop over straws
    }     // end loop over layers

    // We have constructed one sector, Now place the sectors (in the logical device)

    for ( int isec = 0; isec<device.nSectors(); ++isec){

      if ( secDraw > -1 && isec > secDraw ) continue;

      //       cout << "Debugging sector: " << isec << " " << secName << " secDraw: " << secDraw << endl;

      const Sector& sector = device.getSector(isec);

      // place the trapezoid in its position ready for the RZ rotation

      CLHEP::HepRotationZ RZ(sector.boxRzAngle() - device.rotation()); // well we know it is only arround z...

      //       cout << "Debugging sector.boxRzAngle(), device.rotation(): " << sector.boxRzAngle() << " " << 
      //         device.rotation() << endl;      
      
      // a leak?
      // we add an 180deg rotation for even sectors
      G4RotationMatrix* secRotation = ((isec%2)==1) ? 
        new G4RotationMatrix(RXForTrapezoids*RZForTrapezoids*RZ.inverse()):
        new G4RotationMatrix(RXForTrapezoids*RZForTrapezoids*RX2ForTrapezoids*RZ.inverse());

      // origin a.k.a offset wrt current mother volume
      CLHEP::Hep3Vector origin = sector.boxOffset() - device.origin();

      //       cout << "Debugging sector.origin:      "   << isec << " " << secName << origin << endl;
      //       cout << "Debugging sector.boxOffset(): " << isec << " " << secName << sector.boxOffset() << endl;
      
      // we may need to keep those pointers somewhre... (this is only the last one...)

      secInfo.physical =  new G4PVPlacement(secRotation,
                                            origin,
                                            secInfo.logical,
                                            secName,
                                            devInfo.logical,
                                            false,
                                            isec,
                                            doSurfaceCheck);

    }       // end loop over sectors

    // we have constructed one logical device above, we need to place it a few times now...

    for ( int idev=0; idev<ttracker.nDevices(); ++idev ){

      // changes here affect StrawSD

      if ( devDraw > -1 && idev > devDraw ) continue;

      //      cout << "Debugging dev: " << idev << " " << devName << " devDraw: " << devDraw << endl;

      const Device& device = ttracker.getDevice(idev);

      CLHEP::HepRotationZ RZ(-device.rotation()); //It is arround z

      // a leak?
      G4RotationMatrix* devRotation  = new G4RotationMatrix(RZ);

      // cout << "Debugging -device.rotation(): " << -device.rotation() << " " << endl;
      // cout << "Debugging device.origin(): " << device.origin() << " " << endl;

      // could we descend the final hierarchy and set the "true" copy numbers?

      // we may need to keep those pointers somewhre... (this is only the last one...)
      devInfo.physical =  new G4PVPlacement(devRotation,
                                            device.origin(),
                                            devInfo.logical,
                                            devName,
                                            motherInfo.logical,
                                            false,
                                            idev,
                                            doSurfaceCheck);


    } // end loop over devices

    return motherInfo;

  } // end of constructTTrackerv3

} // end namespace mu2e
