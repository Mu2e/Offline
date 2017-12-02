//
// Class to construct a modern (Nov 2017) version of the TTracker with
// more details than in the recent past.
//
// Based on Original by Rob Kutschke and Krzysztof Genser

// David Norvil Brown (Louisville), November 2017
//

#include "ConfigTools/inc/SimpleConfig.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "Mu2eG4/inc/ConstructTTrackerDetail5.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"

#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <cmath>

using namespace std;

mu2e::ConstructTTrackerDetail5::ConstructTTrackerDetail5( VolumeInfo   const& ds3Vac,
                                                  SimpleConfig const& config  ):
  // References to arguments
  _ds3Vac(ds3Vac),
  _config(config),

  // Assorted tools
  _helper(*art::ServiceHandle<G4Helper>()),
  _reg(_helper.antiLeakRegistry()),
  _ttracker(*GeomHandle<TTracker>()),

  // Switches that affect the entire tracker
  _verbosityLevel(_config.getInt("ttracker.verbosityLevel",0)),
  _doSurfaceCheck(_config.getBool("g4.doSurfaceCheck",false)||
                  _config.getBool("ttracker.doSurfaceCheck",false)),
  _forceAuxEdgeVisible(_config.getBool("g4.forceAuxEdgeVisible",false)),

  // Coordinate system transformation
  _offset(computeOffset()),

  // The value to be returned to the code that instantiates this object.
  _motherInfo(){

  // Do the work.
  constructMother(); // Unchanged from "TDR" version
  constructMainSupports(); // Unchanged from "TDR" version
  constructPlanes(); // Much different from "TDR" version <==

  // Debugging only.
  if ( _config.getBool("ttracker.drawAxes",false) ) {
    constructAxes();
  }

}

// Some helper functions and classes.
namespace mu2e {

  namespace {

    // Useful named constants.
    G4ThreeVector const zeroVector(0.,0.,0.);
    bool const place(true);
    bool const doNotPlace(false);
    bool const noSurfaceCheck(false);
    bool const doSurfaceCheck(true);
    G4RotationMatrix * noRotation(nullptr);
    string trackerEnvelopeName("TTrackerPlaneEnvelope");


  } // end anonymous namespace

} // end namespace mu2e

// *****************************************************
// The tracker is not centered within its mother volume.
// Compute the offset that moves positions in the tracker coordinate system to positions in the
// coordinate system of the mother.
CLHEP::Hep3Vector
mu2e::ConstructTTrackerDetail5::computeOffset(){
  return CLHEP::Hep3Vector(0., 0., _ttracker.z0()-_ttracker.mother().position().z() );
}
// *****************************************************


// *****************************
// Construct the TTracker Mother
// *****************************
void
mu2e::ConstructTTrackerDetail5::constructMother(){

  // Parameters of the new style mother volume ( replaces the envelope volume ).
  PlacedTubs const& mother = _ttracker.mother();

  static int const newPrecision = 8;
  static int const newWidth = 14;

  if (_verbosityLevel > 0) {

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

  // All mother/envelope volumes are made of this material.
  G4Material* envelopeMaterial = findMaterialOrThrow(_ttracker.envelopeMaterial());

  _motherInfo = nestTubs( "TrackerMother",
                          mother.tubsParams(),
                          envelopeMaterial,
                          noRotation,
                          mother.position() - _ds3Vac.centerInWorld,
                          _ds3Vac,
                          0,
                          _config.getBool("_ttracker.envelopeVisible",false),
                          G4Colour::Blue(),
                          _config.getBool("_ttracker.envelopeSolid",true),
                          _forceAuxEdgeVisible,
                          place,
                          _doSurfaceCheck
                          );

  if ( _verbosityLevel > 0) {
    double zhl                 = static_cast<G4Tubs*>(_motherInfo.solid)->GetZHalfLength();
    double motherOffsetInMu2eZ = _motherInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
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

} // end constructMother
// *********************


// ************************
// Construct Main Supports
// Unfortunately, the term "support" has been multiply susbscribed
// in the TTracker world.  Here, it means the upstream and downstream
// stiffening rings and inner brass rings.  The other meaning, farther
// below, is of the rings holding the straws.
// ************************
void
mu2e::ConstructTTrackerDetail5::constructMainSupports(){

  SupportStructure const& sup = _ttracker.getSupportStructure();

  bool ttrackerSupportSurfaceCheck = _config.getBool("ttrackerSupport.doSurfaceCheck",false);

  // These rings are the same as TDR version
  for ( auto const& ring : sup.stiffRings() ){

    if ( _verbosityLevel > 0 ) {
      cout << __func__
           << " Support ring position: "
           << ring.name()     << " "
           << ring.position() << " "
           << _motherInfo.centerInWorld << " "
           << ring.position()-_motherInfo.centerInWorld
           << endl;
    }

    VolumeInfo info = nestTubs( ring.name(),
                                ring.tubsParams(),
                                findMaterialOrThrow(ring.materialName()),
                                noRotation,
                                ring.position()-_motherInfo.centerInWorld,
                                _motherInfo,
                                0,
                                _config.getBool("ttracker.envelopeVisible",false),
                                G4Colour::Red(),
                                _config.getBool("ttracker.envelopeSolid",true),
                                _forceAuxEdgeVisible,
                                place,
                                _doSurfaceCheck || ttrackerSupportSurfaceCheck
                                );

  } // end of loop over stiff support rings


  // contructing the support beams

  for ( auto const& sbeam : sup.beamBody() ) {

    if ( _verbosityLevel > 0 ) {
      cout  << __func__
            << " Support Beam Position: "
            << sbeam.name()               << " "
            << sbeam.position()           << " "
            << _motherInfo.centerInWorld  << " "
            << sbeam.position()-_motherInfo.centerInWorld << " "
            << sbeam.tubsParams()
            << ", phi range : "
            << sbeam.phi0()/M_PI*180.
            << ", "
            << sbeam.phiMax()/M_PI*180.
            << endl;
    }

    nestTubs( sbeam.name(),
              sbeam.tubsParams(),
              findMaterialOrThrow(sbeam.materialName()),
              0x0,
              sbeam.position()-_motherInfo.centerInWorld,
              _motherInfo,
              0,
              _config.getBool("ttracker.envelopeVisible",false),
              G4Colour::Cyan(),
              _config.getBool("ttracker.envelopeSolid",true),
              _forceAuxEdgeVisible,
              place,
              _doSurfaceCheck || ttrackerSupportSurfaceCheck
              );

  } // end of loop over support beam

  // here construct the supportServices

  VolumeInfo serviceVI;
  VolumeInfo& ttSSE1 = _helper.locateVolInfo("TTrackerSupportServiceEnvelope_11");
  VolumeInfo& ttSSE2 = _helper.locateVolInfo("TTrackerSupportServiceEnvelope_21");

  // each subsection has en envelope; creating them first

  for ( auto const& sbeam : sup.beamServices() ) {

    // placing the services in the right envelope
    VolumeInfo& ttSSE = ( sbeam.name().find("eSectionEnvelope_1") != string::npos ) ? ttSSE1 : ttSSE2;

    if ( sbeam.name().find("TTrackerSupportServiceSectionEnvelope_") != string::npos ) {

      if ( _verbosityLevel > 0 ) {
        cout << __func__
             << " Support Beam Service eSection Envelope Position: "
             << sbeam.name()               << " "
             << sbeam.position()           << " "
             << _motherInfo.centerInWorld  << " "
             << sbeam.position()-ttSSE.centerInWorld << " "
             << sbeam.tubsParams()
             << endl;
      }

      serviceVI =
        nestTubs( sbeam.name(),
                  sbeam.tubsParams(),
                  findMaterialOrThrow(sbeam.materialName()),
                  0x0,
                  sbeam.position()-ttSSE.centerInWorld,
                  ttSSE,
                  0,
                  _config.getBool("ttracker.envelopeVisible",false),
                  G4Colour::White(),
                  _config.getBool("ttracker.envelopeSolid",true),
                  _forceAuxEdgeVisible,
                  place,
                  _doSurfaceCheck || ttrackerSupportSurfaceCheck
                  );

      if ( _verbosityLevel > 0 ) {
        cout << " Material "
             << serviceVI.logical->GetMaterial()->GetName()
             << " Mass in kg: "
             << serviceVI.logical->GetMass()/CLHEP::kg
             << endl;
      }

    }

  } // end of loop over support services

  // now placing services in their envelopes

  for ( auto const& sbeam : sup.beamServices() ) {

    // placing the services in the right envelope

    std::string stf = "TTrackerSupportService_";

    size_t stfp = sbeam.name().find(stf);

    if ( stfp != string::npos ) {

      std::string sse = "TTrackerSupportServiceSectionEnvelope_" + sbeam.name().substr(stfp+stf.size(),2);

      VolumeInfo& ttSSE =  _helper.locateVolInfo(sse);

      if ( _verbosityLevel > 0 ) {
        cout << __func__
             << " Support Beam Service Position: "
             << sbeam.name()               << " "
             << sse                        << " "
             << sbeam.position()           << " "
             << ttSSE.centerInWorld  << " "
             << sbeam.position()-ttSSE.centerInWorld << " "
             << sbeam.tubsParams()
             << endl;
      }

      serviceVI =
        nestTubs( sbeam.name(),
                  sbeam.tubsParams(),
                  findMaterialOrThrow(sbeam.materialName()),
                  0x0,
                  sbeam.position()-ttSSE.centerInWorld,
                  ttSSE,
                  0,
                  _config.getBool("ttracker.envelopeVisible",false),
                  ( sbeam.name().find("_c") != string::npos ) ? G4Colour::Yellow() : G4Colour::Green(),
                  _config.getBool("ttracker.envelopeSolid",true),
                  _forceAuxEdgeVisible,
                  place,
                  _doSurfaceCheck || ttrackerSupportSurfaceCheck
                  );

      if ( _verbosityLevel > 0 ) {
        cout << " Material "
             << serviceVI.logical->GetMaterial()->GetName()
             << " Mass in kg: "
             << serviceVI.logical->GetMass()/CLHEP::kg
             << endl;
      }

    }

  } // end of loop over support beam services to put in envelopes

  // print the final mass per Service Section envelope

  if ( _verbosityLevel > 0 ) {
    for ( auto const& sbeam : sup.beamServices() ) {
      // prinitnt the final mass per Service Section envelope
      if ( sbeam.name().find("TTrackerSupportServiceSectionEnvelope_") != string::npos ) {

        VolumeInfo& ttSSE =  _helper.locateVolInfo(sbeam.name());

        cout << __func__
             << " Support Beam Service Section Envelope: "
             << sbeam.name()               << " "
             << " Material "
             << ttSSE.logical->GetMaterial()->GetName()
             << " Final Mass in kg: "
             << ttSSE.logical->GetMass(true)/CLHEP::kg
             << endl;
      }

    }

  } // end of printout of support service envelope masses


} // end constructMainSupports
// ***************************


// **********************************************************
// Construct all stations and place them in the mother volume.
// We do not yet represent stations as G4 objects.  Each station 
// is representated in G4 as two independent planes.
// **********************************************************
void
mu2e::ConstructTTrackerDetail5::constructPlanes(){

  // If the following is set, turns on detailed drawing of the "named" planes
  int plnDraw(_config.getInt("ttracker.plnDraw",-1));

  // Build logical volume heirarchy for a single panel.  <==
  // Historically this only considered a panel as the location of straws 
  // - not the manifold, electronics, base ring, etc, collectively known
  // as the "supports" (not to be confused with the "main supports").
  // We are changing that with v5.  Now we make a panel in software
  // that is like what we refer to as a panel in hardware.
  VolumeInfo baseStrawPanel = prepareStrawPanel();

  // Build the electronic board aka key <==
  VolumeInfo baseEBKey       = prepareEBKey(true);   // key itself
  VolumeInfo baseEBKeyShield = prepareEBKey(false);  // key shield

  // Build logical volume heirarchy for all elements of the support 
  // structure that live inside each plane envelope.  <==
  // Here is where the support material for each plane has historically been 
  // built. 

  TubsParams planeEnvelopeParams = _ttracker.getPlaneEnvelopeParams();
  bool planeEnvelopeVisible      = _config.getBool("ttracker.planeEnvelopeVisible",false);
  bool planeEnvelopeSolid        = _config.getBool("ttracker.planeEnvelopeSolid",true);

  G4Material* envelopeMaterial = findMaterialOrThrow(_ttracker.envelopeMaterial());

  double dPhiPanel = panelHalfAzimuth();

  for ( int ipln=0; ipln<_ttracker.nPlanes(); ++ipln ){

    if ( plnDraw > -1 && ipln > plnDraw ) continue;

    if (_verbosityLevel > 0 ) { 
      cout << "Debugging pln: " << ipln << " plnDraw: " << plnDraw << endl;
      cout << __func__ << " working on plane:   " << ipln << endl;
    }

    const Plane& plane = _ttracker.getPlane(ipln);

    if (!plane.exists()) continue;
    if (_verbosityLevel > 0 ) { 
      cout << __func__ << " existing   plane:   " << ipln << endl;
    }

    // plane.origin() is in the detector coordinate system.
    // plnPosition - is in the coordinate system of the Tracker mother volume.
    CLHEP::Hep3Vector plnPosition = plane.origin() + _offset;

    if ( _verbosityLevel > 1 ){ 
      cout << "Debugging -plane.rotation(): " << -plane.rotation() << " " << endl;
      cout << "Debugging plane.origin():    " << plane.origin() << " " << endl;
      cout << "Debugging position in mother: " << plnPosition << endl;
    }


    // We need a new logical volume for each plane envelope - because the panels
    // may be placed differently.  We need a distinct name for each logical volume.
    ostringstream os;
    os << "_"  << std::setfill('0') << std::setw(2) << ipln;

    VolumeInfo plnInfo = nestTubs( trackerEnvelopeName + os.str(),
                                   planeEnvelopeParams,
                                   envelopeMaterial,
                                   noRotation,
                                   plnPosition,
                                   _motherInfo.logical,
                                   ipln,
                                   planeEnvelopeVisible,
                                   G4Colour::Magenta(),
                                   planeEnvelopeSolid,
                                   _forceAuxEdgeVisible,
                                   place,
                                   noSurfaceCheck
                                   );

    if ( _doSurfaceCheck) {
      checkForOverlaps(plnInfo.physical, _config, _verbosityLevel>0);
    }


    addPanelsAndEBKeys( baseStrawPanel, ipln, plnInfo, 
			baseEBKey, baseEBKeyShield,
                        _motherInfo, dPhiPanel );

  } // end loop over planes

} // end constructPlanes


// ************************
// ************************


// Some helper functions and classes.
namespace mu2e {

  namespace {

    // A helper structure for the bookkeeping of nesting many tubs.
    struct VHelper{

      VHelper( PlacedTubs const& apart, G4Colour const& acolour, std::string amotherName, bool avisible=false ):
        part(apart), colour(acolour), info(), motherName(amotherName), visible(avisible) {
      }

      PlacedTubs const& part;
      G4Colour          colour;
      mu2e::VolumeInfo  info;
      std::string       motherName;
      bool              visible;
    };

    G4LogicalVolume* findMotherLogical( std::vector<VHelper> const& vols, std::string const& name ){
      for ( auto const& vol : vols ){
        if ( vol.part.name() == name ){
          return vol.info.logical;
        }
      }
      throw cet::exception("GEOM") << __func__
                                   << " could not find requested mother logical volume: "
                                   << name
                                   << "\n";
    }

  } // end anonymous namespace

} // end namespace mu2e


// ************************
// ************************



// Now the code that makes each panel


// ************************
// ************************

mu2e::VolumeInfo
mu2e::ConstructTTrackerDetail5::preparePanel(const int& iPlane, 
					     const int& iPanel, 
					     VolumeInfo& thePlane, 
					     VolumeInfo& strawPanel,	
					     CLHEP::Hep3Vector& pnlPosition, 
					     G4RotationMatrix* rot ){

  // Detail 5 version will take panels very seriously.
  // In the TDR version of this code, a "panel" was just a set of straws.
  // Then all the support structure for six panels was mushed together as
  // one single plane ring assembly.
  // In this version, a panel is the straws and manifold and support structure
  // and electronics associated with what we actually construct and call
  // a "panel."

  TubsParams panEnvParams = _ttracker.getPanelEnvelopeParams();
  //  SupportStructure const& sup     = _ttracker.getSupportStructure();

  // Panels are identical other than placement - so get required properties from plane 0, panel 0.
  Panel const& panel(_ttracker.getPanel(PanelId(iPlane,iPanel)));

  bool panelEnvelopeVisible = _config.getBool("ttracker.panelEnvelopeVisible",false);
  bool panelEnvelopeSolid   = _config.getBool("ttracker.panelEnvelopeSolid",true);

  // Azimuth of the centerline of the panel envelope: panelCenterPhi
  // Upon construction and before placement, the panel envelope occupies [0,phiMax].
  double panelCenterPhi = panelHalfAzimuth();

  if (_verbosityLevel>0) {
    cout << __func__
         << " preparing panel: "
         << " panelCenterPhi "
         << panelCenterPhi/M_PI*180.
         << " panel center from the maker "
         << panel.boxRzAngle()/M_PI*180.
         << endl;
  }

  // Get information about the channel position and depth.
  //  PlacedTubs const& chanUp(sup.innerChannelUpstream());

  // Internally all panel envelopes are the same.
  // Create one logical panel envelope but, for now, do not place it.
  // Fill it with straws and then place it multiple times.
  //  TubsParams panEnvParams(panelEnvelopeParams.innerRadius(),
  //                      sup.outerRingUpstream().tubsParams().outerRadius(),
  //                      planeEnvelopeParams.halfLength(),
  //                      0.,
  //                      phiMax);


  if (_verbosityLevel>0) {
    cout << __func__
         << " preparing panel: "
         << " panel envelope panEnvParams: "
         << panEnvParams.innerRadius()
         << ", "
         << panEnvParams.outerRadius()
         << ", "
         << panEnvParams.zHalfLength()
         << ", "
         << panEnvParams.phi0()/M_PI*180.
         << ", "
         << panEnvParams.phiMax()/M_PI*180.
         << endl;
  }

  G4Material* envelopeMaterial = findMaterialOrThrow(_ttracker.envelopeMaterial());

  // *****************
  // Here we make the envelope for everything in one panel, including 
  // the "panel supports."  (Not "Main supports").
  // *****************
  //  ostringstream pnlName;
  //  pnlName << "Panel_" << std::setfill('0') << std::setw(2) << iPlane << "_" 
  //	  << std::setw(1) << iPanel;

  // *****************************************
  // Here we will get Alignment information
  // AlignmentMap* alignMap = art::ServiceHandle<AlignmentService>()->alignmentMap();
  // AlignmentSequence tmpSeq = AlignMap->find(panel.name());
  // unsigned long time = 0; // Need a way to get time?
  // CLHEP::Hep3Vector alignTrans(0,0,0);
  // rotation here
  // if ( tmpSeq.isValid() ) {
  //    alignTrans = tmpSeq.getAlignment(time).diplacement();
  //    rotation here
  // }
  // Now modify the position vector and rotation with the Alignment info
  // *********************************************************

  int baseCopyNo = iPlane * _ttracker.getPlane(0).nPanels();

  VolumeInfo pnl0Info = nestTubs( panel.name("Panel"),
                                  panEnvParams,
                                  envelopeMaterial,
                                  rot,
                                  pnlPosition,
                                  thePlane.logical, // logical volume of plane
                                  0,               // copyNo
                                  panelEnvelopeVisible,
                                  G4Colour::Magenta(),
                                  panelEnvelopeSolid,
                                  true,                  // edge visible
                                  place,
                                  _doSurfaceCheck
                                  );

  // ***************************************
  // Now put in the straw panel already prepared
  // ***************************************

  CLHEP::Hep3Vector spOffset(0., 0., 0. ); 
  G4VPhysicalVolume* physVol = new G4PVPlacement( 0, spOffset,
						  strawPanel.logical,
						  strawPanel.name,
						  pnl0Info.logical,
						  false,
						  baseCopyNo + iPanel,
						  false );
  if ( _doSurfaceCheck ) {
    checkForOverlaps( physVol, _config, _verbosityLevel>0);
  }


  return pnl0Info;

} // end preparePanel
// *******************


mu2e::VolumeInfo 
mu2e::ConstructTTrackerDetail5::prepareStrawPanel() {

 // ************************
  // Now we put in the Straws!
  // ************************

  // What we are doing is building one set of straws that we'll use in
  // all of the panels in all of the planes in all the world.


  //  TubsParams planeEnvelopeParams = _ttracker.getPlaneEnvelopeParams();
  SupportStructure const& sup     = _ttracker.getSupportStructure();

  // Straw Panels are identical other than placement - so get required 
  // properties from plane 0, panel 0.
  Plane const& plane(_ttracker.getPlane(PlaneId(0)));
  Panel const& panel(_ttracker.getPanel(PanelId(0,0)));

  bool panelEnvelopeVisible = _config.getBool("ttracker.panelEnvelopeVisible",
					      false);
  bool panelEnvelopeSolid   = _config.getBool("ttracker.panelEnvelopeSolid",
					      true );
  bool strawVisible          = _config.getBool("ttracker.strawVisible",false);
  bool strawSolid            = _config.getBool("ttracker.strawSolid",true);
  bool strawLayeringVisible  = _config.getBool("ttracker.strawLayeringVisible",
					       false);
  bool strawLayeringSolid    = _config.getBool("ttracker.strawLayeringSolid",
					       false);
  bool partialStraws         = _config.getBool("ttracker.partialStraws",false);

  // Azimuth of the centerline of a the panel enveleope: panelCenterPhi
  // Upon construction and before placement, the panel envelope 
  // occupies [0,phiMax].
  double panelCenterPhi = panelHalfAzimuth();

  if (_verbosityLevel>1) {
    cout << __func__
         << " preparing straw panel: "
         << " panelCenterPhi "
         << panelCenterPhi/M_PI*180.
         << " panel center from the maker "
         << panel.boxRzAngle()/M_PI*180.
         << endl;
  }


  // Internally all panel envelopes are the same.
  // Create one logical panel envelope but, for now, do not place it.
  // Fill it with straws and then place it multiple times.
  //  TubsParams panEnvParams(planeEnvelopeParams.innerRadius(),
  //                          sup.innerChannelUpstream().tubsParams().outerRadius(),
  //                          chanUp.tubsParams().zHalfLength(),
  //                          0.,
  //                          phiMax);

  TubsParams panEnvParams = _ttracker.getPanelEnvelopeParams();
  double zCorrection = panEnvParams.zHalfLength() - _ttracker.panelOffset();

  if (_verbosityLevel>0) {
    cout << __func__
         << " preparing straw panel: "
         << " panel envelope panEnvParams: "
         << panEnvParams.innerRadius()
         << ", "
         << panEnvParams.outerRadius()
         << ", "
         << panEnvParams.zHalfLength()
         << ", "
         << panEnvParams.phi0()/M_PI*180.
         << ", "
         << panEnvParams.phiMax()/M_PI*180.
         << endl;
  }


  G4Material* envelopeMaterial = findMaterialOrThrow(
						_ttracker.envelopeMaterial());

  VolumeInfo spnl0Info = nestTubs( "StrawPanelEnvelope",
                                  panEnvParams,
                                  envelopeMaterial,
                                  noRotation,
                                  zeroVector,
                                  nullptr,               // logical volume - not needed since no placement.
                                  0,                     // copyNo
                                  panelEnvelopeVisible,
                                  G4Colour::Magenta(),
                                  panelEnvelopeSolid,
                                  true,                  // edge visible
                                  doNotPlace,
                                  _doSurfaceCheck
                                  );

  // The rotation matrix that will place the straw inside the panel envelope.
  // For straws on the upstream side of the support, the sign of the 
  // X rotation was chosen to put the z axis of the straw volume to be 
  // positive when going clockwise; this is the same convention used
  // within the TTracker class.  For straws on the downstream side of 
  // the support, the z axis of the straw volume is positive when going 
  // counter-clockwise.  This won't be important since we never work
  // in the straw-local coordinates within G4.

  CLHEP::HepRotationZ pnlRz(-panelCenterPhi);
  CLHEP::HepRotationX pnlRx(M_PI/2.);
  G4RotationMatrix* panelRotation = _reg.add(G4RotationMatrix(pnlRx*pnlRz));

  // The z of the center of the placed panel envelope, in Mu2e coordinates.
  // This carries a sign, depending on upstream/downstream.
  double zPanel(0.);
  for ( int i=0; i<panel.nLayers(); ++i){
    zPanel += panel.getStraw(StrawId(0,0,i,0)).getMidPoint().z();
  }
  zPanel /= panel.nLayers();
  // Is panel 0 on the upstream(+1) or downstream(-z) side of the plane.
  double side = (zPanel-plane.origin().z()) > 0. ? -1. : 1.;

  // A unit vector in the direction from the origin to the wire center within the panel envelope.
  CLHEP::Hep3Vector unit( cos(panelCenterPhi), sin(panelCenterPhi), 0.);

  // Sensitive detector object for the straw gas.
  G4VSensitiveDetector *strawSD =
    G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::TrackerGas());
  if ( strawSD == nullptr ){
    cout << __func__
         << " Warning: there is no sensitive detector object for the straw gas.  Continuing ..."
         << endl;
  }

  // Sensitive detector object for the straw sense wires.
  G4VSensitiveDetector *senseWireSD = G4SDManager::GetSDMpointer()->
    FindSensitiveDetector(SensitiveDetectorName::TrackerSWires());

  // Sensitive detector object for the straw walls.
  G4VSensitiveDetector *strawWallSD = G4SDManager::GetSDMpointer()->
    FindSensitiveDetector(SensitiveDetectorName::TrackerWalls());

  // Place the straws into the panel envelope.
  for ( std::vector<Layer>::const_iterator i=panel.getLayers().begin(); i != panel.getLayers().end(); ++i ){

    Layer const& lay(*i);

    for ( std::vector<Straw const*>::const_iterator j=lay.getStraws().begin();
          j != lay.getStraws().end(); ++j ){

      Straw const&       straw(**j);
      StrawDetail const& detail(straw.getDetail());

      if (_verbosityLevel>2) {
        cout << __func__ << " constructing straw "
             << straw.id().getStraw()
             << " id: "
             << straw.id()
             << " with detailIndex of: "
             << straw.detailIndex()
             << endl;
      }

      if ( _verbosityLevel > 2 ){
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

      // To make graphical debugging less busy, create only a subset of the straws.
      if ( partialStraws ){
        if ( straw.id().getStraw()%8 != 0 && straw.id().getStraw() !=47 ) continue;
      }

      // Mid point of the straw in Mu2e coordinates.
      CLHEP::Hep3Vector const& pos(straw.getMidPoint());
      // Mid point of the straw, within the panel envelope.
      double r = (CLHEP::Hep3Vector( pos.x(), pos.y(), 0.)).mag();
      CLHEP::Hep3Vector mid = r*unit;
      mid.setZ(side*(pos.z() - zPanel - zCorrection));

      int copyNo=straw.index().asInt();
      bool edgeVisible(true);

      // The enclosing volume for the straw is made of gas.  The walls and the wire will be placed inside.
      VolumeInfo strawVol =  nestTubs( straw.name("TTrackerStrawGas_"),
                                       detail.getOuterTubsParams(),
                                       findMaterialOrThrow(detail.gasMaterialName()),
                                       panelRotation,
                                       mid,
                                       spnl0Info.logical,
                                       copyNo,
                                       strawVisible,
                                       G4Colour::Red(),
                                       strawSolid,
                                       edgeVisible,
                                       place,
                                       _doSurfaceCheck
                                       );

      if (strawSD) {
        strawVol.logical->SetSensitiveDetector(strawSD);
      }

      // Wall has 4 layers; the three metal layers sit inside the plastic layer.
      // The plastic layer sits inside the gas.
      VolumeInfo wallVol =  nestTubs( straw.name("TTrackerStrawWall_"),
                                      detail.wallMother(),
                                      findMaterialOrThrow(detail.wallMaterialName()),
                                      noRotation,
                                      zeroVector,
                                      strawVol.logical,
                                      copyNo,
                                      strawLayeringVisible,
                                      G4Colour::Green(),
                                      strawLayeringSolid,
                                      edgeVisible,
                                      place,
                                      _doSurfaceCheck
                                      );


      VolumeInfo outerMetalVol =  nestTubs( straw.name("TTrackerStrawWallOuterMetal_"),
                                            detail.wallOuterMetal(),
                                            findMaterialOrThrow(detail.wallOuterMetalMaterialName()),
                                            noRotation,
                                            zeroVector,
                                            wallVol.logical,
                                            copyNo,
                                            strawLayeringVisible,
                                            G4Colour::Blue(),
                                            strawLayeringSolid,
                                            edgeVisible,
                                            place,
                                            _doSurfaceCheck
                                            );


      VolumeInfo innerMetal1Vol =  nestTubs( straw.name("TTrackerStrawWallInnerMetal1_"),
                                             detail.wallInnerMetal1(),
                                             findMaterialOrThrow(detail.wallInnerMetal1MaterialName()),
                                             noRotation,
                                             zeroVector,
                                             wallVol.logical,
                                             copyNo,
                                             strawLayeringVisible,
                                             G4Colour::Blue(),
                                             strawLayeringSolid,
                                             edgeVisible,
                                             place,
                                             _doSurfaceCheck
                                             );

      VolumeInfo innerMetal2Vol =  nestTubs( straw.name("TTrackerStrawWallInnerMetal2_"),
                                             detail.wallInnerMetal2(),
                                             findMaterialOrThrow(detail.wallInnerMetal2MaterialName()),
                                             noRotation,
                                             zeroVector,
                                             wallVol.logical,
                                             copyNo,
                                             strawLayeringVisible,
                                             G4Colour::Blue(),
                                             strawLayeringSolid,
                                             edgeVisible,
                                             place,
                                             _doSurfaceCheck
                                             );

      VolumeInfo wireVol =  nestTubs( straw.name("TTrackerWireCore_"),
                                      detail.wireMother(),
                                      findMaterialOrThrow(detail.wireCoreMaterialName()),
                                      noRotation,
                                      zeroVector,
                                      strawVol.logical,
                                      copyNo,
                                      strawLayeringVisible,
                                      G4Colour::Green(),
                                      strawLayeringSolid,
                                      edgeVisible,
                                      place,
                                      _doSurfaceCheck
                                      );

      VolumeInfo platingVol =  nestTubs( straw.name("TTrackerWirePlate_"),
                                         detail.wirePlate(),
                                         findMaterialOrThrow(detail.wirePlateMaterialName()),
                                         noRotation,
                                         zeroVector,
                                         wireVol.logical,
                                         copyNo,
                                         strawLayeringVisible,
                                         G4Colour::Red(),
                                         strawLayeringSolid,
                                         edgeVisible,
                                         place,
                                         _doSurfaceCheck
                                         );

      if (senseWireSD) {
        wireVol.logical->SetSensitiveDetector(senseWireSD);
        platingVol.logical->SetSensitiveDetector(senseWireSD);
      }

      if (strawWallSD) {
        wallVol.logical->SetSensitiveDetector(strawWallSD);
        outerMetalVol.logical->SetSensitiveDetector(strawWallSD);
        innerMetal1Vol.logical->SetSensitiveDetector(strawWallSD);
        innerMetal2Vol.logical->SetSensitiveDetector(strawWallSD);
      }

    } // end loop over straws within a layer
  } // end loop over layers


  // We have now placed all the straws in the panel.


						 
  // ***************************************
   // Now we put in the so-called support structure.
  // A lot of this code is taken from the original preparePlaneSupports
  // function.  Many of the pieces here have an "upstream"
  // and "downstream" versions because of the way the code
  // for v3 of the TTracker was developed.  We only need one
  // or the other for v5, so arbitrarily choose downstream.
  // In the future, when we retire the ability to choose
  // v3 of the ttracker, we can simplify the support structure class.
  // ***************************************

  bool supportVisible = _config.getBool("ttracker.supportVisible",false);
  bool supportSolid   = _config.getBool("ttracker.supportSolid",true);

  // Many parts of the support structure are G4Tubs objects.
  G4Colour  lightBlue (0.0, 0.0, 0.75);
  std::vector<VHelper> vols;
  vols.reserve(6);

  vols.push_back( VHelper(sup.centerPlate(),         lightBlue,           "",  
			  supportVisible ));
  vols.push_back( VHelper(sup.outerRingDownstream(), G4Colour::Green(),   "",  
			  supportVisible ));
  vols.push_back( VHelper(sup.coverDownstream(),     G4Colour::Cyan(),    "", 
			  supportVisible ));
  vols.push_back( VHelper(sup.gasDownstream(),       G4Colour::Magenta(), "", 
			  supportVisible ));
  vols.push_back( VHelper(sup.g10Downstream(),         G4Colour::Blue(),   
			  sup.gasDownstream().name(),   supportVisible ));
  vols.push_back( VHelper(sup.cuDownstream(),          G4Colour::Yellow(),  
			  sup.gasDownstream().name(),   supportVisible ));


  for ( auto& vol : vols ){

    PlacedTubs const& part = vol.part;

    if ( _verbosityLevel > 0 ) {
      cout << "Building panel support element " << part.name() 
	   << ", with params: " << part.tubsParams() 
	   << ", at position:  " << part.position() << endl;
    }


    if ( vol.motherName.empty() ){

      // These volumes are top level volumes of the support structure 
      // and will be placed inside each panel envelope.
      ostringstream aName;
      aName << part.name();
      vol.info = nestTubs( aName.str(),
			   part.tubsParams(),
			   findMaterialOrThrow(part.materialName()),
			   noRotation,
			   part.position(),
			   spnl0Info.logical,
			   0,
			   vol.visible,
			   vol.colour,
			   supportSolid,
			   _forceAuxEdgeVisible,
			   place,
			   _doSurfaceCheck
			   );

    } else{

      // These volumes are lower level volumes and are placed, at this time, 
      // inside their parents.
      ostringstream aName;
      aName << part.name();

      G4LogicalVolume* motherLogical = findMotherLogical( vols, 
							  vol.motherName);

      vol.info = nestTubs( aName.str(),
                           part.tubsParams(),
                           findMaterialOrThrow(part.materialName()),
                           noRotation,
                           part.position(),
                           motherLogical,
                           0,
                           vol.visible,
                           vol.colour,
                           supportSolid,
                           _forceAuxEdgeVisible,
                           place,
                           _doSurfaceCheck
                           );
    }
  } // end loop over support vols

  {  // Now add the inner ring.

    PlacedTubs const& ring     = sup.innerRing();
    PlacedTubs const& chanDown = sup.innerChannelDownstream();

    // Parameters of the polycone.
    const int n(6);
    double z[n];
    double rIn[n];
    double rOut[n];

    // I am not sure if I am allowed to have two coincident zplanes so slightly
    // slope the edges of the channel.
    double eps=0.1;

    z[0]    = -ring.tubsParams().zHalfLength();
    z[1]    = chanDown.position().z()  - chanDown.tubsParams().zHalfLength() - eps;
    z[2]    = chanDown.position().z() - chanDown.tubsParams().zHalfLength();
    z[3]    = chanDown.position().z() + chanDown.tubsParams().zHalfLength();
    z[4]    = chanDown.position().z() + chanDown.tubsParams().zHalfLength() + eps;
    z[5]    = ring.tubsParams().zHalfLength();

    rIn[0]  = ring.tubsParams().innerRadius();
    rIn[1]  = ring.tubsParams().innerRadius();
    rIn[2]  = chanDown.tubsParams().outerRadius();
    rIn[3]  = chanDown.tubsParams().outerRadius();
    rIn[4]  = ring.tubsParams().innerRadius();
    rIn[5]  = ring.tubsParams().innerRadius();

    for ( int i=0; i<n; ++i){
      rOut[i] = ring.tubsParams().outerRadius();
    }

    ostringstream aName;
    aName << ring.name() << "_a";

    // Inner ring is divided into three parts - two channels where the
    // straws hit it, and a solid ring where they don't.  First
    // part is from phi = 0 to first rib...
    VolumeInfo innerRingInfoa( aName.str(), ring.position(), zeroVector );
    double maxPhi1 = (sup.panelPhiRange() - sup.panelPhiRibs())/2.;
    innerRingInfoa.solid = new G4Polycone( ring.name(),
                                          0.,
                                          maxPhi1,
                                          n,
                                          z,
                                          rIn,
                                          rOut);

    finishNesting( innerRingInfoa,
                   findMaterialOrThrow(ring.materialName()),
                   noRotation,
                   zeroVector,
                   spnl0Info.logical,
                   0,
                   supportVisible,
                   G4Colour::Red(),
                   supportSolid,
                   _forceAuxEdgeVisible,
                   place,
                   false
                   );

    // Second part is from second rib to phi = phiMax...
    ostringstream bName;
    bName << ring.name() << "_b";

    VolumeInfo innerRingInfob( bName.str(), ring.position(), zeroVector );
    innerRingInfob.solid = new G4Polycone( ring.name(),
                                          maxPhi1+sup.panelPhiRibs(),
                                          maxPhi1,
                                          n,
                                          z,
                                          rIn,
                                          rOut);

    finishNesting( innerRingInfob,
                   findMaterialOrThrow(ring.materialName()),
                   noRotation,
                   zeroVector,
                   spnl0Info.logical,
                   0,
                   supportVisible,
                   G4Colour::Red(),
                   supportSolid,
                   _forceAuxEdgeVisible,
                   place,
                   false
                   );

    // Now the third, solid part, from first rib to second rib.
    ostringstream cName;
    cName << ring.name() << "_c";
    TubsParams middlePartParams ( ring.innerRadius(),ring.outerRadius(),
				  ring.zHalfLength(), maxPhi1,
				  sup.panelPhiRibs());
    VolumeInfo ribRingPart = nestTubs( cName.str(),
				       middlePartParams,
				       findMaterialOrThrow(sup.outerRingDownstream().materialName()),
				       noRotation,
				       ring.position(),
				       spnl0Info.logical,
				       0,
				       supportVisible,
				       G4Colour::Red(),
				       supportSolid,
				       _forceAuxEdgeVisible,
				       place,
				       _doSurfaceCheck
				       );
  }
  // Did the inner ring, now we've got it all for the panel!

  return spnl0Info;
} // end of prepareStrawPanel()
  // ***************************************


mu2e::VolumeInfo
mu2e::ConstructTTrackerDetail5::prepareEBKey(bool keyItself){

  // the keys and their shields are almost the same; the boolean decides which one to make

  // place they keys in the addPannels...; use different mother volume, i.e. tracker mother
  // write the corresponding TTrackeMaker code first

  // keys are identical other than placement - so get required properties from plane 0, panel 0.
  // they are placed wrt to the panels; they could be constructed independently though
  //

  bool keyVisible  = _config.getBool("ttrackerSupport.electronics.key.visible",true);
  bool keySolid    = _config.getBool("ttrackerSupport.electronics.key.solid",false);

  // Internally all keys are the same.
  // Create one logical volume for now, do not place it.
  Panel const& panel(_ttracker.getPanel(PanelId(0,0)));

  TubsParams  keyParams  = keyItself ? panel.getEBKeyParams() : panel.getEBKeyShieldParams();

  G4Material* keyMaterial = findMaterialOrThrow( keyItself ?
                                                 panel.getEBKeyMaterial() :
                                                 panel.getEBKeyShieldMaterial() );

  VolumeInfo key0Info = keyItself ? nestTubs("PanelEBKey",
                                             keyParams,
                                             keyMaterial,
                                             noRotation,
                                             zeroVector,
                                             nullptr,       // logical volume - not needed since no placement.
                                             0,             // copyNo
                                             keyVisible,
                                             G4Colour::Yellow(),
                                             keySolid,
                                             true,          // edge visible
                                             doNotPlace,
                                             _doSurfaceCheck
                                             ) :
                                    nestTubs("PanelEBKeyShield",
                                             keyParams,
                                             keyMaterial,
                                             noRotation,
                                             zeroVector,
                                             nullptr,       // logical volume - not needed since no placement.
                                             0,             // copyNo
                                             keyVisible,
                                             G4Colour::Green(),
                                             keySolid,
                                             true,          // edge visible
                                             doNotPlace,
                                             _doSurfaceCheck
                                             );

  if (_verbosityLevel>1) {
    cout << __func__
         << " preparing EBKey: "
         << " EBKey Params: "
         << keyParams.innerRadius()
         << ", "
         << keyParams.outerRadius()
         << ", "
         << keyParams.zHalfLength()
         << ", "
         << keyParams.phi0()/M_PI*180.
         << ", "
         << keyParams.phiMax()/M_PI*180.
         << endl;

  }

  // Make it sensitive here
  if (keyItself) {
    G4VSensitiveDetector* EBKeySD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::panelEBKey());
    if (EBKeySD) {
      key0Info.logical->SetSensitiveDetector(EBKeySD);
    }
  }

  return key0Info;

} // end prepareEBKey


// For debugging graphically, add three boxes that are aligned with the coordinate axes,
// on the tracker axis, just downstream of the tracker.
void
mu2e::ConstructTTrackerDetail5::constructAxes(){

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

  G4Material* envelopeMaterial = findMaterialOrThrow(_ttracker.envelopeMaterial());

  // A box marking the x-axis.
  double halfDimX[3];
  halfDimX[0] = blockLength;
  halfDimX[1] = blockWidth;
  halfDimX[2] = blockWidth;
  CLHEP::Hep3Vector xAxisPos(xOff,0.,z0);
  CLHEP::Hep3Vector axisOffset(0.,0.,12000.);
  xAxisPos += axisOffset;
  nestBox( "xAxis",
           halfDimX,
           envelopeMaterial,
           noRotation,
           xAxisPos,
           _ds3Vac,
           copyNo,
           isVisible,
           G4Colour::Blue(),
           isSolid,
           edgeVisible,
           place,
           _doSurfaceCheck
           );

  // A box marking the y-axis.
  double halfDimY[3];
  halfDimY[0] = blockWidth;
  halfDimY[1] = blockLength;
  halfDimY[2] = blockWidth;
  CLHEP::Hep3Vector yAxisPos(0.,yOff,z0);
  yAxisPos += axisOffset;
  nestBox( "yAxis",
           halfDimY,
           envelopeMaterial,
           noRotation,
           yAxisPos,
           _ds3Vac,
           copyNo,
           isVisible,
           G4Colour::Red(),
           isSolid,
           edgeVisible,
           place,
           _doSurfaceCheck
           );

  // A box marking the z-axis.
  double halfDimZ[3];
  halfDimZ[0] = blockWidth;
  halfDimZ[1] = blockWidth;
  halfDimZ[2] = blockLength;
  CLHEP::Hep3Vector zAxisPos(0.,0.,zOff);
  zAxisPos += axisOffset;
  nestBox( "zAxis",
           halfDimZ,
           envelopeMaterial,
           noRotation,
           zAxisPos,
           _ds3Vac,
           copyNo,
           isVisible,
           G4Colour::Green(),
           isSolid,
           edgeVisible,
           place,
           _doSurfaceCheck
           );

} // end constructAxes


// ******************************
// Here panels get added to plane
// ******************************
void
mu2e::ConstructTTrackerDetail5::addPanelsAndEBKeys(VolumeInfo& baseStrawPanel,
						   int ipln,
						   VolumeInfo& plane,
						   VolumeInfo& baseEBKey,
						   VolumeInfo& baseEBKeyShield,
						   VolumeInfo& trackerMother,
						   double panelCenterPhi ){

  // panel EBKeys are placed inside the tracker mother volume as they
  // are outside the planes
  // we somehow need to avoid the overlaps

  Plane const& pln = _ttracker.getPlane(ipln);
  // to get the key info from the base panel
  Panel const& panel(_ttracker.getPanel(PanelId(0,0)));

  // to prevent the overlaps
  SupportStructure const& sup = _ttracker.getSupportStructure();

  int pnlDraw = _config.getInt ("ttracker.pnlDraw",-1);

  // Assume all planes have the same number of panels.
  int baseCopyNo = ipln * _ttracker.getPlane(0).nPanels();

  // The panel will be placed in the middle of the channel.
  //  PlacedTubs const& chanUp(_ttracker.getSupportStructure().innerChannelUpstream());

  // Place the panel envelope into the plane envelope, one placement per panel.
  for ( int ipnl=0; ipnl<pln.nPanels(); ++ipnl){

    // For debugging, only place the one requested panel.
    if ( pnlDraw > -1  && ipnl > pnlDraw ) continue;
    //if ( pnlDraw > -1  && ipnl%2 == 0 ) continue;

    // Choose a representative straw from this  (plane,panel).
    Straw const& straw = _ttracker.getStraw( StrawId(ipln,ipnl,0,0) );

    // Azimuth of the midpoint of the wire.
    CLHEP::Hep3Vector const& mid = straw.getMidPoint();
    double phimid = mid.phi();
    if ( phimid < 0 ) phimid += 2.*M_PI;

    //    double theZOffset = _ttracker.panelOffset();
    double theZOffset = _ttracker.getPanelEnvelopeParams().zHalfLength();

    // Is this panel on the front(+) or back(-) face of the plane?
    double sign   = ((mid.z() - pln.origin().z())>0. ? 1.0 : -1.0 );
    CLHEP::Hep3Vector panelPosition(0.,0., sign*std::abs(theZOffset) );

    // The rotation that will make the straw mid point have the correct azimuth.
    // Panels on the downstream side are flipped front/back by rotation about Y.
    // so the sense of the z rotation changes sign.

    // the specific rotation below is a consequence of the orientation
    // and origin of the local coordinate system of the G4Tub

    double phi0 = panelCenterPhi-phimid;
    if ( sign > 0 ){
      phi0 = panelCenterPhi + M_PI +phimid;
    }

    CLHEP::HepRotationZ rotZ(phi0);
    CLHEP::HepRotationY rotY(M_PI);
    G4RotationMatrix* rotation  = (sign < 0 )?
      _reg.add(G4RotationMatrix(rotZ)) :
      _reg.add(G4RotationMatrix(rotZ*rotY));

    if (_verbosityLevel>1) {
      cout << __func__
           << " placing panel: "
           << ipnl
           << " inside plane: "
           << ipln
           << " sign "
           << sign
           << " phimid "
           << phimid/M_PI*180.
           << " phi0 "
           << phi0/M_PI*180.
           << " rel position "
           << panelPosition
           << " "
           << endl;
    }

    VolumeInfo placedPan = preparePanel( ipln, ipnl, plane, 
					 baseStrawPanel,
					 panelPosition, rotation );
    // ostringstream disPanelName;
    // disPanelName << "Panel_" << ipln << "_" << ipnl;
    bool many(false);
    // G4VPhysicalVolume* physVol = new G4PVPlacement(rotation,
    //                                                panelPosition,
    //                                                basePanel.logical,
    //                                                disPanelName.str(),
    //                                                plane.logical,
    //                                                many,
    //                                                0,
    //                                                false);
    if ( _doSurfaceCheck) {
      checkForOverlaps( placedPan.physical, _config, _verbosityLevel>0);
    }

    // Now placing the EBkey

    // see comment above about panel phi0
    double keyAngShift = panel.getEBKeyParams().phiMax()*0.5;

    // additional rotation to place some of the keys at the top of the tracker

    double keyPhi0 = keyAngShift - phimid - panel.getEBKeyPhiExtraRotation(); // as it needs to placed per geant4
    if ( sign > 0 ) {
      keyPhi0 = keyAngShift + M_PI + phimid + panel.getEBKeyPhiExtraRotation();
    }

    keyPhi0 = getWithinZeroTwoPi(keyPhi0);

    bool willOverlap = false;

    double keyNominalPhi0 =  getWithinZeroTwoPi(phimid + panel.getEBKeyPhiExtraRotation());

    for ( auto const& sbeam : sup.beamBody() ) {

      // check if there will be overlap with the support

      double phiEnd    =  getWithinZeroTwoPi(sbeam.phi0() + sbeam.phiMax());
      double diffEdge1 = diffWithinZeroTwoPi(keyNominalPhi0,sbeam.phi0());
      double diffEdge2 = diffWithinZeroTwoPi(keyNominalPhi0,phiEnd);

      if ( ( diffEdge1 < panel.getEBKeyParams().phiMax() ) ||
           ( diffEdge2 < panel.getEBKeyParams().phiMax() ) ) {
        willOverlap = true;
      }

      if ( willOverlap && (_verbosityLevel > 0 ) ) {
        cout  << __func__
              << " Support Beam/key overalp, will skip the key: "
              << " plane: "
              << ipln
              << ", key "
              << ipnl
              << ", "
              << sbeam.name()               << " "
              << ", phi0, "
              << sbeam.phi0()/M_PI*180.
              << ", span "
              << sbeam.phiMax()/M_PI*180.
              << ", keyNominalPhi0 "
              << keyNominalPhi0/M_PI*180.
              << ", keyPhi0 "
              << keyPhi0/M_PI*180.
              << ", span "
              << panel.getEBKeyParams().phiMax()
              << ", keyAngShift "
              << keyAngShift/M_PI*180.
               << endl;
      }

      if (willOverlap) break;

    }

    // override willOverlap
    //willOverlap = false;

    if (!willOverlap) {

      CLHEP::HepRotationZ keyRotZ(keyPhi0);
      CLHEP::HepRotationY keyRotY(M_PI);
      G4RotationMatrix* keyRotation  = (sign < 0 )?
        _reg.add(G4RotationMatrix(keyRotZ)) :
        _reg.add(G4RotationMatrix(keyRotZ*keyRotY));

      // it is placed inside the trackerMother not inside the plane

      // panelPosition is in the plane coordinate system
      // pln.origin()  is in the detector/tracker coordinate system.
      // plnPosition - is in the coordinate system of the Tracker mother volume.

      // we need to shift the positions by half of the width of the
      // keys/shield depending on the key location wrt to the plane

      CLHEP::Hep3Vector keyPosition =  panelPosition + pln.origin()+ _offset +
        CLHEP::Hep3Vector(0., 0., panel.getEBKeyShieldParams().zHalfLength());

      if (_verbosityLevel>1) {
        cout << __func__
             << " placing key: "
             << ipnl
             << " inside tracker mother volume "
             << " sign "
             << sign
             << " phi0 " << phi0 << ", "
             << phi0/M_PI*180.
             << " keyPhi0 " << keyPhi0 << ", "
             << keyPhi0/M_PI*180.
             << " keyAngShift " << keyAngShift << ", "
             << keyAngShift/M_PI*180.
             << " positions: "
             <<  _offset << ", "
             << pln.origin() << ", "
             << panelPosition  << ", "
             << keyPosition
             << " tracker Mother position in the world "
             << trackerMother.centerInWorld
             << " tracker Mother position in Mu2e "
             << trackerMother.centerInMu2e()
             << endl;
      }

      // the copy number allows to distinguish the actual key
      // EBKey itself
      G4VPhysicalVolume* keyPhysVol = new G4PVPlacement(keyRotation,
                                                        keyPosition,
                                                        baseEBKey.logical,
                                                        baseEBKey.name,
                                                        trackerMother.logical,
                                                        many,
                                                        baseCopyNo + ipnl,
                                                        false);
      if ( _doSurfaceCheck) {
        checkForOverlaps( keyPhysVol, _config, _verbosityLevel>0);
      }

      CLHEP::Hep3Vector keyShieldPosition =  panelPosition + pln.origin()+ _offset -
        CLHEP::Hep3Vector(0., 0., panel.getEBKeyShieldParams().zHalfLength());

      // EBKeyShield
      G4VPhysicalVolume* keyShieldPhysVol = new G4PVPlacement(keyRotation,
                                                              keyShieldPosition,
                                                              baseEBKeyShield.logical,
                                                              baseEBKeyShield.name,
                                                              trackerMother.logical,
                                                              many,
                                                              baseCopyNo + ipnl,
                                                              false);
      if ( _doSurfaceCheck) {
        checkForOverlaps( keyShieldPhysVol, _config, _verbosityLevel>0);
      }

    }

  }

} // end of addPanelsAndEBKeys 
// ***************************



double
mu2e::ConstructTTrackerDetail5::panelHalfAzimuth(){

  SupportStructure const& sup     = _ttracker.getSupportStructure();
  return sup.panelPhiRange()/2.0;

}


double mu2e::ConstructTTrackerDetail5::getWithinZeroTwoPi (double phi0) {
  phi0 = fmod(phi0,2.*M_PI);
  if ( phi0<0. ) {
    phi0 += 2.*M_PI;
  }
  return phi0;
}

double mu2e::ConstructTTrackerDetail5::diffWithinZeroTwoPi (double phi1, double phi2) {
  phi1 = getWithinZeroTwoPi(phi1);
  phi2 = getWithinZeroTwoPi(phi2);
  double diff = fabs(phi1 - phi2);
  return fmin(diff,2.*M_PI - diff);
}
