//
// Free function to construct the TTracker version that:
//  - Uses version v3 of the G4 model
//  - Uses the support model SupportModel::detailedv0.
//  - See comments at the top of constructTTrackerv3.cc for additional details
//    of the meaning of version numbers.
//
// $Id: constructTTrackerv3Detailed.cc,v 1.7 2014/01/06 20:46:39 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/01/06 20:46:39 $
//
// Contact person Rob Kutschke,
//   - Based on constructTTrackerv3 by KLG
//
// Notes
//     This version makes logical mother volumes per plane and per
//     panel and places panels in plane and straws in panel
//     It has only one panel/plane logical volume placed several times
//     This versoin has a negligeable construction time and a much smaller memory footprint
//

// C++ includes
#include <iomanip>
#include <iostream>
#include <string>

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/constructTTracker.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

#include "cetlib/pow.h"

using namespace std;

namespace mu2e{

namespace {

  // A helper structure to for nesting a many tubs.
  struct VolHelper{

    VolHelper( PlacedTubs const& apart, G4Colour const& acolour, std::string const& amotherName, std::deque<VolHelper> const& vols, bool avisible=false ):
      part(&apart), colour(acolour), motherName(amotherName), info(), mother(0), visible(avisible) {
      for ( std::deque<VolHelper>::const_iterator i=vols.begin(), e=vols.end();
            i != e; ++i ){
        if ( i->part->name() == motherName ){
          mother = &(i->info);
        }
      }
      if ( mother == 0 ){
        throw cet::exception("GEOM") << "For volume: "
                                     << part->name()
                                     << "  Could not find mother: "
                                     << motherName
                                     << "\n";
      }
    }

    VolHelper( PlacedTubs const& apart, G4Colour const& acolour, std::string const& amotherName, mu2e::VolumeInfo const& amother, bool avisible=false ):
      part(&apart), colour(acolour), motherName(amotherName), info(), mother(&amother), visible(avisible) {
    }

    PlacedTubs const* part;
    G4Colour          colour;
    std::string       motherName;
    mu2e::VolumeInfo        info;
    mu2e::VolumeInfo const* mother;
    bool                    visible;
  };

  void addPanels( int ipln,
                   VolumeInfo const& plnInfo,
                   VolumeInfo& pnl0Info,
                   double panelCenterPhi,
                   AntiLeakRegistry& reg,
                   SimpleConfig const& config ){

    TTracker const& ttracker = *(GeomHandle<TTracker>());
    Plane const& pln        = ttracker.getPlane(ipln);

    int pnlDraw         = config.getInt ("ttracker.pnlDraw",-1);
    bool doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false)
      || config.getBool("ttracker.doSurfaceCheck",false);
    int verbosityLevel  = config.getInt ("ttracker.verbosityLevel",0);


    // Needed to get the center position in local z.
    PlacedTubs const& chanUp(ttracker.getSupportStructure().innerChannelUpstream());

    // Place the panel envelope into the plane envelope, one placement per panel.
    for ( int ipnl=0; ipnl<pln.nPanels(); ++ipnl){

      // For debugging, only place the one requested panel.
      if ( pnlDraw > -1  && ipnl != pnlDraw ) continue;
      //if ( pnlDraw > -1  && ipnl%2 == 0 ) continue;

      // Choose a representative straw from this this (plane,panel).
      Straw const& straw = ttracker.getStraw( StrawId(ipln,ipnl,0) );

      // Azimuth of the midpoint of the wire.
      CLHEP::Hep3Vector const& mid = straw.getMidPoint();
      double phimid = mid.phi();
      if ( phimid < 0 ) phimid += 2.*M_PI;

      // Is this panel on the front or back face of the plane?
      double sign   = ((mid.z() - pln.origin().z())>0. ? 1.0 : -1.0 );
      CLHEP::Hep3Vector panelPosition(0.,0., sign*std::abs(chanUp.position().z()) );

      // The rotation that will make the straw mid point have the correct azimuth.
      // Panels on the downstream side are flipped front/back by rotation about Y.
      // so the sense of the z rotation changes sign.
      double phi0 = panelCenterPhi-phimid;
      if ( sign > 0 ){
        phi0 = panelCenterPhi + M_PI +phimid;
      }

      CLHEP::HepRotationZ rotZ(phi0);
      CLHEP::HepRotationY rotY(M_PI);
      G4RotationMatrix* rotation  = (sign < 0 )?
        reg.add(G4RotationMatrix(rotZ)) :
        reg.add(G4RotationMatrix(rotZ*rotY));

      bool many(false);
      pnl0Info.physical =  new G4PVPlacement(rotation,
                                             panelPosition,
                                             pnl0Info.logical,
                                             pnl0Info.name,
                                             plnInfo.logical,
                                             many,
                                             ipnl,
                                             false);
      if ( doSurfaceCheck) {
        checkForOverlaps(pnl0Info.physical, config, verbosityLevel>0);
      }

    }

  } // end addPanels

} // end anonymous namespace


  VolumeInfo constructTTrackerv3Detailed( VolumeInfo const& parent,
                                          SimpleConfig const& config,
                                          SensitiveDetectorHelper const& sdHelper
                                          ){

    cout << "Mark Detailed. " << endl;

    G4Helper& _helper      = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    int verbosityLevel = config.getInt("ttracker.verbosityLevel",0);

    // Control of graphics for debugging the geometry.
    // Only instantiate panels to be drawn.
    int plnDraw = config.getInt("ttracker.plnDraw",-1);
    //int pnlDraw = config.getInt("ttracker.pnlDraw",-1);
    bool const doSurfaceCheck = config.getBool("g4.doSurfaceCheck",false)
      || config.getBool("ttracker.doSurfaceCheck",false);
    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);

    G4ThreeVector const zeroVector(0.0,0.0,0.0);
    bool place(true);
    bool doNotPlace(false);
    G4RotationMatrix* noRotation(0);

    // Master geometry for the TTracker.
    TTracker const & ttracker = *(GeomHandle<TTracker>());

    // Parameters of the new style mother volume ( replaces the envelope volume ).
    PlacedTubs const& mother = ttracker.mother();

    // Make the envelope volume that holds the tracker planes - not the end rings and staves.
    //TubsParams envelopeParams = ttracker.getInnerTrackerEnvelopeParams();

    static int const newPrecision = 8;
    static int const newWidth = 14;

    if (verbosityLevel > 0) {
      int oldPrecision = cout.precision(newPrecision);
      int oldWidth = cout.width(newWidth);
      std::ios::fmtflags oldFlags = cout.flags();
      cout.setf(std::ios::fixed,std::ios::floatfield);
      cout << "Debugging tracker env envelopeParams ir,or,zhl,phi0,phimax:            " <<
	"   " <<
	mother.innerRadius() << ", " <<
	mother.outerRadius() << ", " <<
	mother.zHalfLength() << ", " <<
	mother.phi0()        << ", " <<
	mother.phiMax()      << ", " <<
	endl;
      cout.setf(oldFlags);
      cout.precision(oldPrecision);
      cout.width(oldWidth);
    }

    //G4ThreeVector trackerOffset( 0., 0., ttracker.z0() );
    // Offset of the center of the tracker within its mother volume.

    CLHEP::Hep3Vector motherOffset(0., 0., ttracker.z0()-mother.position().z() );
    cout << "Centers: " << mother.position() << " "
         << parent.centerInWorld << " " << ttracker.z0() << endl;
    cout << "         " << motherOffset << endl;

    // All mother/envelope volumes are made of this material.
    G4Material* envelopeMaterial = findMaterialOrThrow(ttracker.envelopeMaterial());

    VolumeInfo motherInfo = nestTubs( "TrackerMother",
                                      mother.tubsParams(),
                                      envelopeMaterial,
                                      noRotation,
                                      mother.position() - parent.centerInWorld,
                                      parent,
                                      0,
                                      config.getBool("ttracker.envelopeVisible",false),
                                      G4Colour::Blue(),
                                      config.getBool("ttracker.envelopeSolid",true),
                                      forceAuxEdgeVisible,
                                      place,
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

    // Now place the endRings and Staves.
    // FixME: add the cut-outs for services within the staves.
    if ( ttracker.getSupportModel() == SupportModel::detailedv0 ) {

      SupportStructure const& sup = ttracker.getSupportStructure();

      for ( auto const& ring : sup.stiffRings() ){

        cout << "Ring Position: "
             << ring.position() << " "
             << motherInfo.centerInWorld << " "
             << ring.position()-motherInfo.centerInWorld
             << endl;

        VolumeInfo info = nestTubs( ring.name(),
                                    ring.tubsParams(),
                                    findMaterialOrThrow(ring.materialName()),
                                    noRotation,
                                    ring.position()-motherInfo.centerInWorld,
                                    motherInfo,
                                    0,
                                    config.getBool("ttracker.envelopeVisible",false),
                                    G4Colour::Red(),
                                    config.getBool("ttracker.envelopeSolid",true),
                                    forceAuxEdgeVisible,
                                    place,
                                    doSurfaceCheck
                                    );

      }

      for ( auto const& sbeam : sup.beamBody() ) {

        if ( verbosityLevel > 0 ) {
          cout << "Support Beam Position: "
               << sbeam.name()               << " "
               << sbeam.position()           << " "
               << motherInfo.centerInWorld  << " "
               << sbeam.position()-motherInfo.centerInWorld << " "
               << sbeam.tubsParams()
               << endl;
        }

        nestTubs( sbeam.name(),
                  sbeam.tubsParams(),
                  findMaterialOrThrow(sbeam.materialName()),
                  0x0,
                  sbeam.position()-motherInfo.centerInWorld,
                  motherInfo,
                  0,
                  config.getBool("ttracker.envelopeVisible",false),
                  G4Colour::Cyan(),
                  config.getBool("ttracker.envelopeSolid",true),
                  forceAuxEdgeVisible,
                  place,
                  true
                  // _doSurfaceCheck
                  );

      }

    }

    TubsParams planeEnvelopeParams = ttracker.getPlaneEnvelopeParams();

    bool planeEnvelopeVisible = config.getBool("ttracker.planeEnvelopeVisible",false);
    bool planeEnvelopeSolid   = config.getBool("ttracker.planeEnvelopeSolid",true);
    bool supportVisible        = config.getBool("ttracker.supportVisible",false);
    bool supportSolid          = config.getBool("ttracker.supportSolid",true);
    bool panelEnvelopeVisible = config.getBool("ttracker.panelEnvelopeVisible",false);
    bool panelEnvelopeSolid   = config.getBool("ttracker.panelEnvelopeSolid",true);
    bool strawVisible          = config.getBool("ttracker.strawVisible",false);
    bool strawSolid            = config.getBool("ttracker.strawSolid",true);
    //bool strawLayeringVisible  = config.getBool("ttracker.strawLayeringVisible",false);
    // bool strawLayeringSolid    = config.getBool("ttracker.strawLayeringSolid",false);
    //bool drawAxes              = config.getBool("ttracker.drawAxes",false);

    //bool ttrackerActiveWr_Wl_SD = config.getBool("ttracker.ActiveWr_Wl_SD",false);

    // The planes differ only in the rotations of their panels (panels).
    // Construct one plane but do not place it until its insides have been populated.

    string trackerEnvelopeName("TTrackerPlaneEnvelope");
    VolumeInfo plnInfo = nestTubs( trackerEnvelopeName,
                                   planeEnvelopeParams,
                                   envelopeMaterial,
                                   noRotation,
                                   zeroVector,
                                   0,
                                   0,
                                   planeEnvelopeVisible,
                                   G4Colour::Magenta(),
                                   planeEnvelopeSolid,
                                   forceAuxEdgeVisible,
                                   doNotPlace,
                                   doSurfaceCheck
                                   );
    if ( verbosityLevel > 0 ){
      cout << "Plane Envelope parameters: " << planeEnvelopeParams << endl;
    }

    // Temporary while debugging the new mother volume.
    if ( config.getBool("rkk.disableTT",false ) ){
      cout << "Will not create TTracker because of a user request " << endl;
      return motherInfo;
    }


    // Many parts of the support structure are G4Tubs objects.  Place all of them,
    SupportStructure const& sup = ttracker.getSupportStructure();
    G4Colour  lightBlue (0.0, 0.0, 0.75);
    std::deque<VolHelper> vols;
    vols.push_back( VolHelper(sup.centerPlate(),         lightBlue,           trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.outerRingUpstream(),   G4Colour::Green(),   trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.outerRingDownstream(), G4Colour::Green(),   trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.coverUpstream(),       G4Colour::Cyan(),    trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.coverDownstream(),     G4Colour::Cyan(),    trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.gasUpstream(),         G4Colour::Magenta(), trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.gasDownstream(),       G4Colour::Magenta(), trackerEnvelopeName,        plnInfo, supportVisible ));
    vols.push_back( VolHelper(sup.g10Upstream(),         G4Colour::Blue(),    sup.gasUpstream().name(),   vols, supportVisible ));
    vols.push_back( VolHelper(sup.g10Downstream(),       G4Colour::Blue(),    sup.gasDownstream().name(), vols, supportVisible ));
    vols.push_back( VolHelper(sup.cuUpstream(),          G4Colour::Yellow(),  sup.gasUpstream().name(),   vols, supportVisible ));
    vols.push_back( VolHelper(sup.cuDownstream(),        G4Colour::Yellow(),  sup.gasDownstream().name(), vols, supportVisible ));

    for ( std::deque<VolHelper>::iterator i=vols.begin(), e=vols.end();
          i != e; ++i ){
      PlacedTubs const& part = *i->part;

      bool placePhysicalVolume = true;
      i->info = nestTubs( part.name(),
                          part.tubsParams(),
                          findMaterialOrThrow(part.materialName()),
                          noRotation,
                          part.position(),
                          i->mother->logical,
                          0,
                          i->visible,
                          i->colour,
                          supportSolid,
                          forceAuxEdgeVisible,
                          placePhysicalVolume,
                          doSurfaceCheck
                          );

    }

    // Pick one of the tubs that represents mocked-up electronics and make it a senstive detector.
    G4VSensitiveDetector *sd = (sdHelper.enabled(StepInstanceName::ttrackerDS)) ?
      G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::TTrackerPlaneSupport()) : nullptr;

    for ( std::deque<VolHelper>::iterator i=vols.begin(), e=vols.end();
          i != e; ++i ){
      PlacedTubs const& part = *i->part;
      if ( part.name().find("TTrackerSupportElecCu") != string::npos ){
        if(sd) i->info.logical->SetSensitiveDetector(sd);
      }
    }

    // The last part of the support structure is the inner ring; Construct it as a polycone
    // that includes the two channels to receive the straws.
    // FixME: use the other c'tor of Polycone that takes points in the rz plane.
    //        This will let us get rid of eps,below.
    {

      PlacedTubs const& ring     = sup.innerRing();
      PlacedTubs const& chanUp   = sup.innerChannelUpstream();
      PlacedTubs const& chanDown = sup.innerChannelDownstream();

      // Parameters of the polycone.
      const int n(10);
      double z[n];
      double rIn[n];
      double rOut[n];

      // FixME: see above about alternate c'tor.
      // I am not sure if the zplanes of a polycone are permitted to have zero spacing.
      // So include a minimum z space at the sharp edges.
      double eps=0.1;

      z[0]    = -ring.tubsParams().zHalfLength();
      z[1]    = chanUp.position().z()   - chanUp.tubsParams().zHalfLength() - eps;
      z[2]    = chanUp.position().z()   - chanUp.tubsParams().zHalfLength();
      z[3]    = chanUp.position().z()   + chanUp.tubsParams().zHalfLength();
      z[4]    = chanUp.position().z()   + chanUp.tubsParams().zHalfLength() + eps;
      z[5]    = chanDown.position().z() - chanDown.tubsParams().zHalfLength() - eps;
      z[6]    = chanDown.position().z() - chanDown.tubsParams().zHalfLength();
      z[7]    = chanDown.position().z() + chanDown.tubsParams().zHalfLength();
      z[8]    = chanDown.position().z() + chanDown.tubsParams().zHalfLength() + eps;
      z[9]    = ring.tubsParams().zHalfLength();

      rIn[0]  = ring.tubsParams().innerRadius();
      rIn[1]  = ring.tubsParams().innerRadius();
      rIn[2]  = chanUp.tubsParams().outerRadius();
      rIn[3]  = chanUp.tubsParams().outerRadius();
      rIn[4]  = ring.tubsParams().innerRadius();
      rIn[5]  = ring.tubsParams().innerRadius();
      rIn[6]  = chanUp.tubsParams().outerRadius();
      rIn[7]  = chanUp.tubsParams().outerRadius();
      rIn[8]  = ring.tubsParams().innerRadius();
      rIn[9]  = ring.tubsParams().innerRadius();

      for ( int i=0; i<n; ++i){
        rOut[i] = ring.tubsParams().outerRadius();
      }

      VolumeInfo innerRingInfo( ring.name(), ring.position(), ring.position() + plnInfo.centerInWorld );
      innerRingInfo.solid = new G4Polycone( ring.name(),
                                            0.,
                                            2.*M_PI,
                                            n,
                                            z,
                                            rIn,
                                            rOut);

      bool visible = true;
      bool placePV = true;
      finishNesting( innerRingInfo,
                     findMaterialOrThrow(ring.materialName()),
                     0,
                     ring.position(),
                     plnInfo.logical,
                     0,
                     visible,
                     G4Colour::Red(),
                     supportSolid,
                     forceAuxEdgeVisible,
                     placePV,
                     doSurfaceCheck
                     );


    }

    // Construct one panel envelope.
    // In this model, panel envelopes are arcs of a disk.

    // Panels are identical other than placement - so get required properties from plane 0, panel 0.
    Panel const& pnl00(ttracker.getPanel(StrawId(0,0,0)));

    // Azimuth of the centerline of a the panel enveleope: panelCenterPhi
    // Upon construction and before placement, the panel envelope occupies [0,phiMax].
    double r1              = planeEnvelopeParams.innerRadius();
    double r2              = sup.innerChannelUpstream().tubsParams().outerRadius();
    double halfChord       = sqrt ( r2*r2 - r1*r1 );
    double panelCenterPhi = atan2(halfChord,r1);
    double phiMax          = 2.*panelCenterPhi;

    // Get information about the channel position and depth.
    PlacedTubs const& chanUp(sup.innerChannelUpstream());

    // Internally all panel envelopes are the same.
    // Create one logical panel envelope but, for now, do not place it.
    // Fill it with straws and then place it multiple times.
    TubsParams s0(r1, r2, chanUp.tubsParams().zHalfLength(), 0., phiMax);
    VolumeInfo pnl0Info = nestTubs( "PanelEnvelope0",
                                    s0,
                                    findMaterialOrThrow(chanUp.materialName()),
                                    0,
                                    zeroVector,
                                    plnInfo.logical,
                                    0,                     // copyNo
                                    panelEnvelopeVisible,
                                    G4Colour::Magenta(),
                                    panelEnvelopeSolid,
                                    true,                  // edge visible
                                    doNotPlace,
                                    doSurfaceCheck
                                    );


    // The rotation matrix that will place the straw inside the panel envelope.
    // For straws on the upstream side of the support, the sign of the X rotation was chosen to put the
    // z axis of the straw volume to be positive when going clockwise; this is the same convention used
    // within the TTracker class.  For straws on the downstream side of the support, the z axis of the straw
    // volume is positive when going counter-clockwise.  This won't be important since we never work
    // in the straw-local coordiates within G4.
    CLHEP::HepRotationZ pnlRz(-panelCenterPhi);
    CLHEP::HepRotationX pnlRx(M_PI/2.);
    G4RotationMatrix* pnl00Rotation = reg.add(G4RotationMatrix(pnlRx*pnlRz));

    // The z of the center of the placed panel envelope, in Mu2e coordinates.
    // This carries a sign, depending on upstream/downstream.
    double zPanel(0.);
    for ( int i=0; i<pnl00.nLayers(); ++i){
      zPanel += pnl00.getStraw(StrawId(0,0,i)).getMidPoint().z();
    }
    zPanel /= pnl00.nLayers();

    // A unit vector in the direction from the origin to the wire center within the panel envelope.
    CLHEP::Hep3Vector unit( cos(panelCenterPhi), sin(panelCenterPhi), 0.);

    // Place the straws into the panel envelope.
    for (const auto straw_p : pnl00.getStrawPointers() ) {
        Straw const&       straw(*straw_p);

        // if ( lay.id().getLayer() != 0 ) continue;
        if ( straw.id().getLayer() != 0 ) continue; // fixme: why a single layer?

        StrawDetail const& detail(straw.getDetail());

        //if ( straw.id().getStraw()%8 != 0 ) continue;

        // Mid point of the straw in Mu2e coordinates.
        CLHEP::Hep3Vector const& pos(straw.getMidPoint());

        // Mid point of the straw, within the panel envelope.
        double r = (CLHEP::Hep3Vector( pos.x(), pos.y(), 0.)).mag();
        CLHEP::Hep3Vector mid = r*unit;
        mid.setZ(pos.z() - zPanel);

        int copyNo=straw.index().asInt();
        bool edgeVisible(true);
        bool placeIt(true);

        // The enclosing volume for the straw is made of gas.  The walls and the wire will be placed inside.
        VolumeInfo strawVol =  nestTubs( straw.name("TTrackerStrawGas_"),
                                         detail.getOuterTubsParams(),
                                         findMaterialOrThrow(detail.gasMaterialName()),
                                         pnl00Rotation,
                                         mid,
                                         pnl0Info.logical,
                                         copyNo,
                                         strawVisible,
                                         G4Colour::Red(),
                                         strawSolid,
                                         edgeVisible,
                                         placeIt,
                                         doSurfaceCheck
                                         );
        /*
        // Make the gas volume of this straw a sensitive detector.
        G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->
          FindSensitiveDetector(SensitiveDetectorName::TrackerGas());
        if(sd) strawVol.logical->SetSensitiveDetector(sd);


        // Wall has 4 layers; the three metal layers sit inside the plastic layer.
        // The plastic layer sits inside the gas.
        VolumeInfo wallVol =  nestTubs( straw.name("TTrackerStrawWall_"),
                                        detail.wallMother(),
                                        findMaterialOrThrow(detail.wallMaterialName()),
                                        0,
                                        zeroVector,
                                        strawVol.logical,
                                        copyNo,
                                        strawLayeringVisible,
                                        G4Colour::Green(),
                                        strawLayeringSolid,
                                        edgeVisible,
                                        placeIt,
                                        doSurfaceCheck
                                        );

        VolumeInfo outerMetalVol =  nestTubs( straw.name("TTrackerStrawWallOuterMetal_"),
                                              detail.wallOuterMetal(),
                                              findMaterialOrThrow(detail.wallOuterMetalMaterialName()),
                                              0,
                                              zeroVector,
                                              wallVol.logical,
                                              copyNo,
                                              strawLayeringVisible,
                                              G4Colour::Blue(),
                                              strawLayeringSolid,
                                              edgeVisible,
                                              placeIt,
                                              doSurfaceCheck
                                              );

        VolumeInfo innerMetal1Vol =  nestTubs( straw.name("TTrackerStrawWallInnerMetal1_"),
                                               detail.wallInnerMetal1(),
                                               findMaterialOrThrow(detail.wallInnerMetal1MaterialName()),
                                               0,
                                               zeroVector,
                                               wallVol.logical,
                                               copyNo,
                                               strawLayeringVisible,
                                               G4Colour::Blue(),
                                               strawLayeringSolid,
                                               edgeVisible,
                                               placeIt,
                                               doSurfaceCheck
                                               );

        VolumeInfo innerMetal2Vol =  nestTubs( straw.name("TTrackerStrawWallInnerMetal2_"),
                                               detail.wallInnerMetal2(),
                                               findMaterialOrThrow(detail.wallInnerMetal2MaterialName()),
                                               0,
                                               zeroVector,
                                               wallVol.logical,
                                               copyNo,
                                               strawLayeringVisible,
                                               G4Colour::Blue(),
                                               strawLayeringSolid,
                                               edgeVisible,
                                               placeIt,
                                               doSurfaceCheck
                                               );

        VolumeInfo wireVol =  nestTubs( straw.name("TTrackerWireCore_"),
                                        detail.wireMother(),
                                        findMaterialOrThrow(detail.wireCoreMaterialName()),
                                        0,
                                        zeroVector,
                                        strawVol.logical,
                                        copyNo,
                                        strawLayeringVisible,
                                        G4Colour::Green(),
                                        strawLayeringSolid,
                                        edgeVisible,
                                        placeIt,
                                        doSurfaceCheck
                                        );

        VolumeInfo platingVol =  nestTubs( straw.name("TTrackerWirePlate_"),
                                           detail.wirePlate(),
                                           findMaterialOrThrow(detail.wirePlateMaterialName()),
                                           0,
                                           zeroVector,
                                           wireVol.logical,
                                           copyNo,
                                           strawLayeringVisible,
                                           G4Colour::Red(),
                                           strawLayeringSolid,
                                           edgeVisible,
                                           placeIt,
                                           doSurfaceCheck
                                           );
        if (ttrackerActiveWr_Wl_SD) {
	        G4VSensitiveDetector *sd = G4SDManager::GetSDMpointer()->
                                FindSensitiveDetector(SensitiveDetectorName::TrackerSWires());
                if (sd) {
                        wireVol.logical->SetSensitiveDetector(sd);
                        platingVol.logical->SetSensitiveDetector(sd);
                }

                sd = nullptr;
                sd = G4SDManager::GetSDMpointer()->
                                FindSensitiveDetector(SensitiveDetectorName::TrackerWalls());
                if (sd) {
                        wallVol.logical->SetSensitiveDetector(sd);
                        outerMetalVol.logical->SetSensitiveDetector(sd);
                        innerMetal1Vol.logical->SetSensitiveDetector(sd);
                        innerMetal2Vol.logical->SetSensitiveDetector(sd);
                }
        }
        */

        if ( verbosityLevel > 1 ){
          cout << "Detail for: " << straw.id() << " " << detail.Id()            << endl;
          cout << "           outerTubsParams: " << detail.getOuterTubsParams() << detail.gasMaterialName()             << endl;
          cout << "           wallMother:      " << detail.wallMother()         << detail.wallMotherMaterialName()      << endl;
          cout << "           wallOuterMetal:  " << detail.wallOuterMetal()     << detail.wallOuterMetalMaterialName()  << endl;
          cout << "           wallCore         " << detail.wallCore()           << detail.wallCoreMaterialName()        << endl;
          cout << "           wallInnerMetal1: " << detail.wallInnerMetal1()    << detail.wallInnerMetal1MaterialName() << endl;
          cout << "           wallInnerMetal2: " << detail.wallInnerMetal2()    << detail.wallInnerMetal2MaterialName() << endl;
          cout << "           wireMother:      " << detail.wireMother()         << detail.wireMotherMaterialName()      << endl;
          cout << "           wirePlate:       " << detail.wirePlate()          << detail.wirePlateMaterialName()       << endl;
          cout << "           wireCore:        " << detail.wireCore()           << detail.wireCoreMaterialName()        << endl;
        }

    } // end loop over straws within a panel


    // Place the plane envelopes into the tracker mother.
    for ( int ipln=0; ipln<ttracker.nPlanes(); ++ipln ){

      // changes here affect StrawSD

      //if ( plnDraw > -1 && ipln != plnDraw ) continue;
      if ( plnDraw > -1 && ipln > plnDraw ) continue;

      verbosityLevel > 1 &&
        cout << "Debugging pln: " << ipln << " " << plnInfo.name << " plnDraw: " << plnDraw << endl;

      const Plane& plane = ttracker.getPlane(ipln);

      // plnPosition - is in the coordinate system of the Tracker mother volume.
      // plane.origin() is in the detector coordinate system.
      CLHEP::Hep3Vector plnPosition = plane.origin()+motherOffset;

      if ( verbosityLevel > 1 ){
        cout << "Debugging -plane.rotation(): " << -plane.rotation() << " " << endl;
        cout << "Debugging plane.origin():    " << plane.origin() << " " << endl;
        cout << "Debugging position in mother: " << plnPosition << endl;
      }

      // could we descend the final hierarchy and set the "true" copy numbers?

      // we may need to keep those pointers somewhre... (this is only the last one...)
      bool pMany(false);
      bool pSurfChk(false);
      plnInfo.physical =  new G4PVPlacement(noRotation,
                                            plnPosition,
                                            plnInfo.logical,
                                            plnInfo.name,
                                            motherInfo.logical,
                                            pMany,
                                            ipln,
                                            pSurfChk);
      if ( doSurfaceCheck) {
         checkForOverlaps(plnInfo.physical, config, verbosityLevel>0);
      }

      addPanels ( ipln, plnInfo, pnl0Info, panelCenterPhi, reg, config );

    } // end loop over planes


    // The section below this is for graphical debugging, not part of the TTracker construction.

    // Draw three boxes to mark the coordinate axes.
    // Dimensions and positions chosen to look right when drawing the first plane or the first station.
    cout << "Doing draw axes: " << config.getBool("ttracker.drawAxes",false) << endl;
    if ( config.getBool("ttracker.drawAxes",false) ) {

      double blockWidth(20.);
      double blockLength(100.);
      double extra(10.);
      double z0   = blockWidth + extra;
      double xOff = blockLength + blockWidth + extra;
      double yOff = blockLength + blockWidth + extra;
      double zOff = z0 + blockLength + blockWidth + extra;

      int copyNo(0);
      int isVisible(true);
      int isSolid(true);
      int edgeVisible(true);
      int placeIt(true);

      // A box marking the x-axis.
      double halfDimX[3];
      halfDimX[0] = blockLength;
      halfDimX[1] = blockWidth;
      halfDimX[2] = blockWidth;
      CLHEP::Hep3Vector xAxisPos(xOff,0.,z0);
      CLHEP::Hep3Vector trackerOffset(0.,0.,11800.);
      xAxisPos += trackerOffset;
      nestBox( "xAxis",
               halfDimX,
               envelopeMaterial,
               noRotation,
               xAxisPos,
               parent,
               copyNo,
               isVisible,
               G4Colour::Blue(),
               isSolid,
               edgeVisible,
               placeIt,
               doSurfaceCheck
               );

      // A box marking the y-axis.
      double halfDimY[3];
      halfDimY[0] = blockWidth;
      halfDimY[1] = blockLength;
      halfDimY[2] = blockWidth;
      CLHEP::Hep3Vector yAxisPos(0.,yOff,z0);
      yAxisPos += trackerOffset;
      nestBox( "yAxis",
               halfDimY,
               envelopeMaterial,
               noRotation,
               yAxisPos,
               parent,
               copyNo,
               isVisible,
               G4Colour::Red(),
               isSolid,
               edgeVisible,
               placeIt,
               doSurfaceCheck
               );

      // A box marking the z-axis.
      double halfDimZ[3];
      halfDimZ[0] = blockWidth;
      halfDimZ[1] = blockWidth;
      halfDimZ[2] = blockLength;
      CLHEP::Hep3Vector zAxisPos(0.,0.,zOff);
      zAxisPos += trackerOffset;
      nestBox( "zAxis",
               halfDimZ,
               envelopeMaterial,
               noRotation,
               zAxisPos,
               parent,
               copyNo,
               isVisible,
               G4Colour::Green(),
               isSolid,
               edgeVisible,
               placeIt,
               doSurfaceCheck
               );

      cout << "Pos: " << xAxisPos << " " << yAxisPos << " "<< zAxisPos << endl;
    }

    return motherInfo;

  } // end of constructTTrackerv3Detailed

} // end namespace mu2e
