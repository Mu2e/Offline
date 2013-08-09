//
// Free function to construct version 3 of the TTracker
//
// $Id: constructTTrackerv3.cc,v 1.31 2013/08/09 01:37:41 genser Exp $
// $Author: genser $
// $Date: 2013/08/09 01:37:41 $
//
// Original author KLG based on RKK's version using different methodology
//
// Notes
//
// 1)  The v3 in this function name says that this is the third way we
//     have implemented a single TTracker design in G4.  It does not refer
//     to alternate designs of the TTracker.
//
//     This version makes logical mother volumes per device and per
//     sector and places sectors in device and straws in sector
//     It has only one sector/device logical volume placed several times
//     This version has a negligeable construction time and a much smaler memory footprint
//
// 2) This function can build the TTracker designs described in:
//      Mu2eG4/test/ttracker_meco.txt - The MECO design, uniform plane spacing
//      Mu2eG4/test/ttracker_v0.txt   - The first Aseet version, pairs of planes form stations
//                                      but one layer of straws per panel (called a sector in this code)
//      Mu2eG4/test/ttracker_v1.txt   - v0 but with with two layers of straws per panel
//      Mu2eG4/test/ttracker_v2.txt   - Adjust spacings to match Mu2e-doc-888-v2.
//
// 3) This function does not know how to build the TTracker described in:
//       Mu2eG4/test/ttracker_v3.txt - Detail support model and detailed layering of straws
//    This geometry can be detected by the method by
//
//      if ( ttracker.getSupportModel() == SupportModel::detailedv0 ) ....
//
//    If this geometry is detected, this function call through to constructTTrackerv3Detailed.cc

// C++ includes
#include <iomanip>
#include <iostream>
#include <string>

// Framework includes
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/constructTTracker.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "TTrackerGeom/inc/TTracker.hh"

// G4 includes
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"


using namespace std;

namespace mu2e{

  VolumeInfo constructTTrackerv3( VolumeInfo const& mother,
                                  double zOff,
                                  SimpleConfig const& config ){

    // Master geometry for the TTracker.
    TTracker const & ttracker = *(GeomHandle<TTracker>());

    // The more detailed version has its own function.
    if ( ttracker.getSupportModel() == SupportModel::detailedv0 ) {
      return constructTTrackerv3Detailed(mother, zOff, config);
    }

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    int verbosityLevel = config.getInt("ttracker.verbosityLevel",0);

    // Control of graphics for debugging the geometry.
    // Only instantiate sectors to be drawn.
    int deviceDraw = config.getInt("ttracker.devDraw",-1);
    int sectorDraw = config.getInt("ttracker.secDraw",-1);
    bool doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false);
    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);

    G4ThreeVector const zeroVector(0.0,0.0,0.0);


    // Make the envelope volume that holds the tracker devices - not the end rings and staves.

    // The devices are now called planes in the CDR

    TubsParams envelopeParams = ttracker.getInnerTrackerEnvelopeParams();

    static int const newPrecision = 8;
    static int const newWidth = 14;

    if (verbosityLevel > 0) {
      int oldPrecision = cout.precision(newPrecision);
      int oldWidth = cout.width(newWidth);
      std::ios::fmtflags oldFlags = cout.flags();
      cout.setf(std::ios::fixed,std::ios::floatfield);
      cout << "Debugging tracker env envelopeParams ir,or,zhl,phi0,phimax:            " <<
	"   " <<
	envelopeParams.innerRadius() << ", " <<
	envelopeParams.outerRadius() << ", " <<
	envelopeParams.zHalfLength() << ", " <<
	envelopeParams.phi0()        << ", " <<
	envelopeParams.phiMax()      << ", " <<
	endl;
      cout.setf(oldFlags);
      cout.precision(oldPrecision);
      cout.width(oldWidth);
    }

    G4ThreeVector trackerOffset( 0., 0., ttracker.z0()-zOff );

    G4Material* envelopeMaterial = findMaterialOrThrow(ttracker.envelopeMaterial());

    VolumeInfo motherInfo = nestTubs( "TrackerMother",
                                      envelopeParams,
                                      envelopeMaterial,
                                      0,
                                      trackerOffset,
                                      mother,
                                      0,
                                      config.getBool("ttracker.envelopeVisible",false),
                                      G4Colour::Blue(),
                                      config.getBool("ttracker.envelopeSolid",true),
                                      forceAuxEdgeVisible,
                                      true,
                                      doSurfaceCheck
                                      );

    if ( verbosityLevel > 0) {
      double zhl         = static_cast<G4Tubs*>(motherInfo.solid)->GetZHalfLength();
      double motherOffsetInMu2eZ = motherInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
      int oldPrecision = cout.precision(3);
      std::ios::fmtflags oldFlags = cout.flags();
      cout.setf(std::ios::fixed,std::ios::floatfield);
      cout << __func__ << " motherOffsetZ           in Mu2e    : " <<
        motherOffsetInMu2eZ << endl;
      cout << __func__ << " mother         Z extent in Mu2e    : " <<
        motherOffsetInMu2eZ - zhl << ", " << motherOffsetInMu2eZ + zhl << endl;
      cout.setf(oldFlags);
      cout.precision(oldPrecision);
    }

    TubsParams deviceEnvelopeParams = ttracker.getDeviceEnvelopeParams();

    bool ttrackerDeviceEnvelopeVisible = config.getBool("ttracker.deviceEnvelopeVisible",false);
    bool ttrackerDeviceEnvelopeSolid   = config.getBool("ttracker.deviceEnvelopeSolid",true);
    bool ttrackerSupportVisible        = config.getBool("ttracker.supportVisible",false);
    bool ttrackerSupportSolid          = config.getBool("ttracker.supportSolid",true);
    bool ttrackerSectorEnvelopeVisible = config.getBool("ttracker.sectorEnvelopeVisible",false);
    bool ttrackerSectorEnvelopeSolid   = config.getBool("ttracker.sectorEnvelopeSolid",true);
    bool ttrackerStrawVisible          = config.getBool("ttracker.strawVisible",false);
    bool ttrackerStrawSolid            = config.getBool("ttracker.strawSolid",true);

    bool ttrackerActiveWr_Wl_SD        = config.getBool("ttracker.ActiveWr_Wl_SD",false);


    // will construct one panel=sector in its nominal position
    // in the new language the device is called a plane ( with two faces ) 
    // then stations have n=2 planes

    // some specific g4 rotations related to the volume type and direction of their axis
    static double const tRAngle  = M_PI_2;
    static double const tRAngle2 = M_PI;
    CLHEP::HepRotationX RXForTrapezoids(tRAngle);
    CLHEP::HepRotationX RX2ForTrapezoids(tRAngle2);
    CLHEP::HepRotationY RYForTrapezoids(tRAngle);
    CLHEP::HepRotationZ RZForTrapezoids(tRAngle);

    VolumeInfo sectorInfo;

    size_t idev = 0;

    const Device& device = ttracker.getDevice(idev);

    // place straws etc... wrt the envelope

    // create a "sector" volume

    // Construct One sector logical volume (and then place it N times)

    const size_t isec = 0;

    const Sector& sector = device.getSector(isec);

    // constructing sector envelope

    // Make a logical volume for this sector,

    // G4IntersectionSolid of G4Box and G4Trd to avoid overlaps of two envelopes

    // reuse device attributes for now

    // get the length of the innermost straw
    Layer const& layer0        = sector.getLayer(0);
    Straw const& straw0        = layer0.getStraw(0);
    StrawDetail const& detail0 = straw0.getDetail();

    verbosityLevel > 0 &&
      cout << "Debugging sector box isec detail0.halfLength(): " << detail0.halfLength() << endl;


    sectorInfo.name = "TTrackerSectorEnvelope";

    G4Box* secBox = new G4Box(sectorInfo.name+"Box",
                              detail0.halfLength(),
                              sector.boxHalfLengths()[2],
                              sector.boxHalfLengths()[1]
                              );

    G4Trd* secTrd = new G4Trd(sectorInfo.name+"Trd",
                              sector.boxHalfLengths()[4],
                              sector.boxHalfLengths()[3],
                              sector.boxHalfLengths()[2],
                              sector.boxHalfLengths()[2],
                              sector.boxHalfLengths()[1]
                              );

    // one could also intersect it with a ring to decrease its radial spread

    sectorInfo.solid =
      new G4IntersectionSolid(sectorInfo.name, secBox, secTrd);

    // false for placing physical volume, just create a logical one
    finishNesting(sectorInfo,
                  envelopeMaterial,
                  0,
                  zeroVector, // this is the "canonical" position, but it does not matter as there is no placement
                  0,
                  0,
                  ttrackerSectorEnvelopeVisible,
                  G4Colour::Cyan(),
                  ttrackerSectorEnvelopeSolid,
                  forceAuxEdgeVisible,
                  false,
                  doSurfaceCheck
                  );

    if (verbosityLevel > 0 ){
      int oldPrecision = cout.precision(newPrecision);
      int oldWidth = cout.width(newWidth);
      std::ios::fmtflags oldFlags = cout.flags();
      cout.setf(std::ios::fixed,std::ios::floatfield);
      cout << "Debugging sector box isec, sector.boxHalfLengths().at(4,3,2,2,1): " <<
        isec << ", " <<
        sector.boxHalfLengths().at(4) << ", " <<
        sector.boxHalfLengths().at(3) << ", " <<
        sector.boxHalfLengths().at(2) << ", " <<
        sector.boxHalfLengths().at(2) << ", " <<
        sector.boxHalfLengths().at(1) << ", " <<
        endl;
      cout.setf(oldFlags);
      cout.precision(oldPrecision);
      cout.width(oldWidth);
    }

    G4RotationMatrix* rotTub = reg.add(G4RotationMatrix(RYForTrapezoids));

    for ( int ilay =0; ilay<sector.nLayers(); ++ilay ){

      verbosityLevel > 1 &&   cout << "Debugging constructTTrackerv3 ilay: " << ilay << endl;

      const Layer& layer = sector.getLayer(ilay);

      for ( int istr=0; istr<layer.nStraws(); ++istr ){

        // "second" layer will have fewer straws (for now) also see TTrackerMaker
        // no, it complicates StrawSD and TTrackerMaker
        // if( ilay%2==1 && istr+1 == layer.nStraws() ) break;

        const Straw& straw = layer.getStraw(istr);

        StrawDetail const& detail = straw.getDetail();

        TubsParams strawWallParams( 0.0, detail.outerRadius(), detail.halfLength() );
        TubsParams strawGasParams ( 0.0, detail.innerRadius(), detail.halfLength() );
        TubsParams strawWireParams( 0.0, detail.wireRadius(),  detail.halfLength() );

        // we are placing the straw w.r.t the trapezoid...
        // the trapezoid aka device envelope has a different coordinate system x->z, z->y, y->x

        G4ThreeVector const mid(straw.getMidPoint().y() - sector.boxOffset().y(),
                                straw.getMidPoint().z() - sector.boxOffset().z(),
                                straw.getMidPoint().x() - sector.boxOffset().x());

        G4ThreeVector const zeroVector(0.0,0.0,0.0);

        if ( verbosityLevel > 2 ) {

          cout << "Debugging istr: " << istr <<
            " mid: " << mid <<
            ", straw.MidPoint " << straw.getMidPoint() <<
            ", sector.boxOffset " <<  sector.boxOffset() <<
            ", device.origin " << device.origin() <<
            endl;

          cout << "Debugging istr: " << istr << " mid: " <<
            mid << ", halflenght " << detail.halfLength() << endl;

          // look at StrawSD to see how the straw index is reconstructed

          cout << "Debugging straw.id(), straw.index() " <<
            straw.id() << ", " << straw.index() << endl;

        }

        // make the straws more distinguishable when displayed
        G4Colour wallColor = (ilay%2 == 1) ?
          ((istr%2 == 0) ? G4Colour::Green() : G4Colour::Yellow()) :
          ((istr%2 == 0) ? G4Colour::Red() : G4Colour::Blue());

        G4Colour gasColor = (ilay%2 == 0) ?
          ((istr%2 == 0) ? G4Colour::Green() : G4Colour::Yellow()) :
          ((istr%2 == 0) ? G4Colour::Red() : G4Colour::Blue());

        G4Colour wireColor = G4Colour::Cyan();

        if (verbosityLevel > 2) {
          int oldPrecision = cout.precision(newPrecision);
          int oldWidth = cout.width(newWidth);
          std::ios::fmtflags oldFlags = cout.flags();
          cout.setf(std::ios::fixed,std::ios::floatfield);
          cout << "Debugging Straw istr, RYForTrapezoids, midpoint: " <<
            istr << ", " << RYForTrapezoids << ", " <<
            mid << ", " <<
            endl;
          cout.setf(oldFlags);
          cout.precision(oldPrecision);
          cout.width(oldWidth);
        }

        VolumeInfo strawWallInfo  = nestTubs(straw.name("TTrackerStrawWall_"),
                                             strawWallParams,
                                             findMaterialOrThrow(detail.wallMaterialName() ),
                                             rotTub,
                                             mid,
                                             sectorInfo.logical,
                                             straw.index().asInt(),
                                             ttrackerStrawVisible,
                                             wallColor,
                                             ttrackerStrawSolid,
                                             forceAuxEdgeVisible,
                                             true,
                                             doSurfaceCheck
                                             );

        // may be not all straws have to be in the volInfo registry, another param?
        // we use the Straw name facility to make them unique

        VolumeInfo strawGasInfo  = nestTubs(straw.name("TTrackerStrawGas_"),
                                            strawGasParams,
                                            findMaterialOrThrow(detail.gasMaterialName()),
                                            0,
                                            zeroVector,
                                            strawWallInfo.logical,
                                            straw.index().asInt(),
                                            ttrackerStrawVisible,
                                            gasColor,
                                            ttrackerStrawSolid,
                                            forceAuxEdgeVisible,
                                            true,
                                            doSurfaceCheck
                                            );

        VolumeInfo strawWireInfo  = nestTubs(straw.name("TTrackerStrawWire_"),
                                             strawWireParams,
                                             findMaterialOrThrow(detail.wireMaterialName()),
                                             0,
                                             zeroVector,
                                             strawGasInfo.logical,
                                             straw.index().asInt(),
                                             ttrackerStrawVisible,
                                             wireColor,
                                             ttrackerStrawSolid,
                                             forceAuxEdgeVisible,
                                             true,
                                             doSurfaceCheck
                                             );

        // Make gas of this straw a sensitive detector.

        G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->
          FindSensitiveDetector(SensitiveDetectorName::TrackerGas());
        if(sd) strawGasInfo.logical->SetSensitiveDetector(sd);

        if (ttrackerActiveWr_Wl_SD) {
          G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->
            FindSensitiveDetector(SensitiveDetectorName::TrackerSWires());
          if(sd) strawWireInfo.logical->SetSensitiveDetector(sd);

          sd = nullptr;
          sd = G4SDManager::GetSDMpointer()->
            FindSensitiveDetector(SensitiveDetectorName::TrackerWalls());
          if (sd) strawWallInfo.logical->SetSensitiveDetector(sd);
        }

      }   // end loop over straws
    }     // end loop over layers

    // We have constructed one sector, 

    // Now construct the devices and place the sectors in them

    // nestSomething creates a physical volume, we need a logical one first
    // false for placing physical volume, just create a logical one
    VolumeInfo deviceInfo  = nestTubs( "TTrackerDeviceEnvelope",
                                       deviceEnvelopeParams,
                                       envelopeMaterial,
                                       0,
                                       zeroVector, // nominal position
                                       0,
                                       0,
                                       ttrackerDeviceEnvelopeVisible,
                                       G4Colour::Magenta(),
                                       ttrackerDeviceEnvelopeSolid,
                                       forceAuxEdgeVisible,
                                       false, /* we are not placing the volume yet */
                                       doSurfaceCheck
                                       );

  // placing the support material

  TubsParams ttrackerDeviceSupportParams = ttracker.getSupportParams().getTubsParams();

  G4Colour  lightBlue (0.0, 0.0, 0.75);
  VolumeInfo supportInfo = nestTubs( "TTrackerDeviceSupport",
                                     ttrackerDeviceSupportParams,
                                     findMaterialOrThrow(ttracker.getSupportParams().materialName()),
                                     0,
                                     zeroVector,
                                     deviceInfo.logical,
                                     0, // all device supports have the same number here
                                     ttrackerSupportVisible,
                                     lightBlue,
                                     ttrackerSupportSolid,
                                     forceAuxEdgeVisible,
                                     true,
                                     doSurfaceCheck
                                     );

  if ( verbosityLevel > 0) {
    cout << "TTrackerDeviceSupport params: "
         << ttrackerDeviceSupportParams.innerRadius() << " "
         << ttrackerDeviceSupportParams.outerRadius() << " "
         << ttrackerDeviceSupportParams.zHalfLength() << " "
         << endl;
  }

  // Make TTrackerDeviceSupport a sensitive detector for radiation damage studies

  G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->
  FindSensitiveDetector(SensitiveDetectorName::TTrackerDeviceSupport());
  if(sd) supportInfo.logical->SetSensitiveDetector(sd);


  if (verbosityLevel > 0 ) {
    int oldPrecision = cout.precision(newPrecision);
    int oldWidth = cout.width(newWidth);
    std::ios::fmtflags oldFlags = cout.flags();
    cout.setf(std::ios::fixed,std::ios::floatfield);
    cout << "Debugging device env idev, deviceEnvelopeParams ir,or,zhl,phi0,phimax: " <<
      idev << ", " <<
      deviceEnvelopeParams.innerRadius() << ", " <<
      deviceEnvelopeParams.outerRadius() << ", " <<
      deviceEnvelopeParams.zHalfLength() << ", " <<
      deviceEnvelopeParams.phi0()        << ", " <<
      deviceEnvelopeParams.phiMax()      << ", " <<
      endl;
    cout.setf(oldFlags);
    cout.precision(oldPrecision);
    cout.width(oldWidth);
  }

  // we now place sectors in the devices and devices in their mother


  // idev can't be size_t here as deviceDraw can be -1
  for ( int idev=0; idev<ttracker.nDevices(); ++idev ){

    // changes here may affect StrawSD

    if ( deviceDraw > -1 && idev > deviceDraw )  continue;

    verbosityLevel > 1 &&
      cout << "Debugging dev: " << idev << " " << deviceInfo.name << " deviceDraw: " << deviceDraw << endl;

    const Device& device = ttracker.getDevice(idev);

    verbosityLevel > 1 &&
      cout << "Debugging -device.rotation(): " << -device.rotation() << " " << endl;
    verbosityLevel > 1 &&
      cout << "Debugging device.origin(): " << device.origin() << " " << endl;

    // isec can't be size_t here as sectorDraw can be -1
    for ( int isec = 0; isec<device.nSectors(); ++isec){

      if ( sectorDraw > -1 && isec > sectorDraw ) continue;

      verbosityLevel > 1 &&
        cout << "Debugging sector: " << isec << " " << sectorInfo.name << " sectorDraw: " << sectorDraw << endl;

      const Sector& sector = device.getSector(isec);

      // place the trapezoid in its position ready for the RZ rotation

      CLHEP::HepRotationZ sectorRZ(sector.boxRzAngle() - device.rotation()); // we know it is only arround z...
      // this is a relative rotation and this is what we need to calculate relative positions

      // it is probably the safest to recalculate offsets from the
      // nominal horizontal position and rotations and ignore absolute
      // positions provided by the geometry service

      verbosityLevel > 1 &&
        cout << "Debugging sector.boxRzAngle(), device.rotation(), diff:   " 
             << sector.boxRzAngle()/M_PI*180. << ", "
             << device.rotation()/M_PI*180.   << ", " 
             << ((sector.boxRzAngle() - device.rotation())/M_PI)*180. << endl;

//       verbosityLevel > 1 &&
//         cout << "Debugging sector rotation: " << isec << " " << sectorInfo.name << " " << 
//         sector.rotation() << endl;

      // we add a 180deg rotation for even sectors
      G4RotationMatrix* sectorRotation = ((isec%2)==1) ?
        reg.add(G4RotationMatrix(RXForTrapezoids*RZForTrapezoids*sectorRZ.inverse())):
        reg.add(G4RotationMatrix(RXForTrapezoids*RZForTrapezoids*RX2ForTrapezoids*sectorRZ.inverse()));

      // origin a.k.a offset wrt current mother volume
      CLHEP::Hep3Vector sectorOrigin = sector.boxOffset() - device.origin();

      CLHEP::Hep3Vector nominalRelPos(CLHEP::Hep3Vector(sectorOrigin.x(),sectorOrigin.y(),0.).mag(), 
                                      0., sectorOrigin.z());

      CLHEP::Hep3Vector sectorRelOrigin = sectorRZ*nominalRelPos; 
      // we still need to do a complemetary rotation

      verbosityLevel > 1 &&
        cout << "Debugging device.origin:      " << isec << " " << deviceInfo.name << device.origin() << endl;
      verbosityLevel > 1 &&
        cout << "Debugging sector.origin:      " << isec << " " << sectorInfo.name << sectorOrigin << endl;
      verbosityLevel > 1 &&
        cout << "Debugging nominalRelPos:      " << isec << " " << sectorInfo.name << nominalRelPos << endl;
      verbosityLevel > 1 &&
        cout << "Debugging sectorRelOrigin:    " << isec << " " << sectorInfo.name << sectorRelOrigin << endl;
      verbosityLevel > 1 &&
        cout << "Debugging sector.boxOffset(): " << isec << " " << sectorInfo.name << sector.boxOffset() << endl;

      // we may need to keep those pointers somewhre... (this is only the last one...)

       sectorInfo.physical =  new G4PVPlacement(sectorRotation,
                                               sectorRelOrigin,
                                               sectorInfo.logical,
                                               sectorInfo.name,
                                               deviceInfo.logical,
                                               false,
                                               isec,
                                               doSurfaceCheck);

    } // end loop over sectors

    CLHEP::HepRotationZ deviceRZ(-device.rotation()); //It is arround z

    G4RotationMatrix* deviceRotation  = reg.add(G4RotationMatrix(deviceRZ));

    verbosityLevel > 1 &&
      cout << "Debugging placing device: " << idev << " " << device.origin() << " " 
           << deviceInfo.name << endl;

    // could we descend the final hierarchy and set the "true" copy numbers?
    // we may need to keep those pointers somewhre... (this is only the last one...)

    deviceInfo.physical =  new G4PVPlacement(deviceRotation,
                                             device.origin(),
                                             deviceInfo.logical,
                                             deviceInfo.name,
                                             motherInfo.logical,
                                             false,
                                             idev,
                                             doSurfaceCheck);
    
    verbosityLevel > 1 &&
      cout << "Debugging placed device: " << idev << " " << idev%2 << " " 
           << deviceInfo.name << endl;

  } // end loop over devices

  return motherInfo;

} // end of constructTTrackerv3

} // end namespace mu2e
