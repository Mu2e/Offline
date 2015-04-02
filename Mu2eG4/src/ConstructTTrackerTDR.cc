// Todo:
// x  0) Refactor support structure so that one logical is placed multiple times.
// x  1) Check sensitive detector for straw gas volume - move code outside of the loop
// x  2) Check sensitive detector for supports
// x  3) Uncomment all of straw structure.
// z  4) Use proper polycone c'tor
//    5) Make figures for docs
//    6) CE fit test
// x  7) Add verbosity switch to all printout
// x  8) Check dimensions
// x  9) Make ttracker_v3.txt work on its own.
//
// Class to construct the TDR version of the TTracker
//
// $Id: ConstructTTrackerTDR.cc,v 1.3 2014/04/11 04:42:06 genser Exp $
// $Author: genser $
// $Date: 2014/04/11 04:42:06 $
//
// Original author Rob Kutschke
//

#include "ConfigTools/inc/SimpleConfig.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"
#include "Mu2eG4/inc/ConstructTTrackerTDR.hh"
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

using namespace std;

mu2e::ConstructTTrackerTDR::ConstructTTrackerTDR( VolumeInfo   const& ds3Vac,
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
  _doSurfaceCheck(_config.getBool("g4.doSurfaceCheck",false)),
  _forceAuxEdgeVisible(_config.getBool("g4.forceAuxEdgeVisible",false)),

  // Coordinate system transformation
  _offset(computeOffset()),

  // The value to be returned to the code that instantiates this object.
  _motherInfo(){

  // Do the work.
  constructMother();
  constructMainSupports();
  constructStations();

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
    string trackerEnvelopeName("TTrackerDeviceEnvelope");


  } // end anonymous namespace

} // end namespace mu2e

// The tracker is not centered within its mother volume.
// Compute the offset that moves positions in the tracker coordinate system to positions in the
// coordinate system of the mother.
CLHEP::Hep3Vector
mu2e::ConstructTTrackerTDR::computeOffset(){
  return CLHEP::Hep3Vector(0., 0., _ttracker.z0()-_ttracker.mother().position().z() );
}

void
mu2e::ConstructTTrackerTDR::constructMother(){

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

void
mu2e::ConstructTTrackerTDR::constructMainSupports(){

  SupportStructure const& sup = _ttracker.getSupportStructure();

  for ( auto const& ring : sup.stiffRings() ){

    if ( _verbosityLevel > 0 ) {
      cout << "Support ring position: "
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
                                _doSurfaceCheck
                                );

  }

  for ( auto const& stave : sup.staveBody() ){

    if ( _verbosityLevel > 0 ) {
      cout << "Stave Position: "
           << stave.name()               << " "
           << stave.position()           << " "
           << _motherInfo.centerInWorld  << " "
           << stave.position()-_motherInfo.centerInWorld
           << endl;
    }

    nestTubs( stave.name(),
              stave.tubsParams(),
              findMaterialOrThrow(stave.materialName()),
              &stave.rotation(),
              stave.position()-_motherInfo.centerInWorld,
              _motherInfo,
              0,
              _config.getBool("ttracker.envelopeVisible",false),
              G4Colour::Yellow(),
              _config.getBool("ttracker.envelopeSolid",true),
              _forceAuxEdgeVisible,
              place,
              _doSurfaceCheck
              );


  }

} // end constructMainSupports


// Construct all stations and place them in the mother volume.
// We do not yet represent stations as G4 objects.  Each station is representated
// in G4 as two independent planes (devices).
void
mu2e::ConstructTTrackerTDR::constructStations(){

  int devDraw(_config.getInt("ttracker.devDraw",-1));

  // Build logical volume heirarchy for a single panel (sector).
  VolumeInfo basePanel = preparePanel();

  // Build logical volume heirarchy for all elements of the support structure that live
  // inside each device envelope.
  std::vector<VolumeInfo> baseSupports;
  preparePlaneSupports(baseSupports);

  TubsParams deviceEnvelopeParams = _ttracker.getDeviceEnvelopeParams();
  bool deviceEnvelopeVisible      = _config.getBool("ttracker.deviceEnvelopeVisible",false);
  bool deviceEnvelopeSolid        = _config.getBool("ttracker.deviceEnvelopeSolid",true);

  G4Material* envelopeMaterial = findMaterialOrThrow(_ttracker.envelopeMaterial());

  double dPhiSector = sectorHalfAzimuth();

  for ( int idev=0; idev<_ttracker.nDevices(); ++idev ){

    if ( devDraw > -1 && idev > devDraw ) continue;

    if (_verbosityLevel > 0 ) {
      cout << "Debugging dev: " << idev << " devDraw: " << devDraw << endl;
      cout << __func__ << " working on device:   " << idev << endl;
    }

    const Device& device = _ttracker.getDevice(idev);

    if (!device.exists()) continue;
    if (_verbosityLevel > 0 ) {
      cout << __func__ << " existing   device:   " << idev << endl;
    }

    // device.origin() is in the detector coordinate system.
    // devPosition - is in the coordinate system of the Tracker mother volume.
    CLHEP::Hep3Vector devPosition = device.origin()+ _offset;

    if ( _verbosityLevel > 1 ){
      cout << "Debugging -device.rotation(): " << -device.rotation() << " " << endl;
      cout << "Debugging device.origin():    " << device.origin() << " " << endl;
      cout << "Debugging position in mother: " << devPosition << endl;
    }


    // We need a new logical volume for each device envelope - because the sectors
    // may be placed differently.  We need a distinct name for each logical volume.
    ostringstream os;
    os << "_" << idev;

    VolumeInfo devInfo = nestTubs( trackerEnvelopeName + os.str(),
                                   deviceEnvelopeParams,
                                   envelopeMaterial,
                                   noRotation,
                                   devPosition,
                                   _motherInfo.logical,
                                   idev,
                                   deviceEnvelopeVisible,
                                   G4Colour::Magenta(),
                                   deviceEnvelopeSolid,
                                   _forceAuxEdgeVisible,
                                   place,
                                   noSurfaceCheck
                                   );

    if ( _doSurfaceCheck) {
      checkForOverlaps(devInfo.physical, _config, _verbosityLevel>0);
    }

    addPlaneSupports( baseSupports, idev, devInfo );

    addPanels ( basePanel, idev, devInfo.logical, dPhiSector );

  } // end loop over devices

} // end constructStations

mu2e::VolumeInfo
mu2e::ConstructTTrackerTDR::preparePanel(){

  TubsParams deviceEnvelopeParams = _ttracker.getDeviceEnvelopeParams();
  SupportStructure const& sup     = _ttracker.getSupportStructure();

  // Sectors are identical other than placement - so get required properties from device 0, sector 0.
  Sector const& sec00(_ttracker.getSector(SectorId(0,0)));

  bool sectorEnvelopeVisible = _config.getBool("ttracker.sectorEnvelopeVisible",false);
  bool sectorEnvelopeSolid   = _config.getBool("ttracker.sectorEnvelopeSolid",true);
  bool strawVisible          = _config.getBool("ttracker.strawVisible",false);
  bool strawSolid            = _config.getBool("ttracker.strawSolid",true);
  bool strawLayeringVisible  = _config.getBool("ttracker.strawLayeringVisible",false);
  bool strawLayeringSolid    = _config.getBool("ttracker.strawLayeringSolid",false);
  bool partialStraws         = _config.getBool("ttracker.partialStraws",false);

  // Azimuth of the centerline of a the sector enveleope: sectorCenterPhi
  // Upon construction and before placement, the sector envelope occupies [0,phiMax].
  double sectorCenterPhi = sectorHalfAzimuth();
  double phiMax          = 2.*sectorCenterPhi;

  // Get information about the channel position and depth.
  PlacedTubs const& chanUp(sup.innerChannelUpstream());

  // Internally all sector envelopes are the same.
  // Create one logical sector envelope but, for now, do not place it.
  // Fill it with straws and then place it multiple times.
  TubsParams secEnvParams(deviceEnvelopeParams.innerRadius(),
                          sup.innerChannelUpstream().tubsParams().outerRadius(),
                          chanUp.tubsParams().zHalfLength(),
                          0.,
                          phiMax);

  G4Material* envelopeMaterial = findMaterialOrThrow(_ttracker.envelopeMaterial());

  VolumeInfo sec0Info = nestTubs( "SectorEnvelope",
                                  secEnvParams,
                                  envelopeMaterial,
                                  noRotation,
                                  zeroVector,
                                  nullptr,               // logical volume - not needed since no placement.
                                  0,                     // copyNo
                                  sectorEnvelopeVisible,
                                  G4Colour::Magenta(),
                                  sectorEnvelopeSolid,
                                  true,                  // edge visible
                                  doNotPlace,
                                  _doSurfaceCheck
                                  );


  // The rotation matrix that will place the straw inside the sector envelope.
  // For straws on the upstream side of the support, the sign of the X rotation was chosen to put the
  // z axis of the straw volume to be positive when going clockwise; this is the same convention used
  // within the TTracker class.  For straws on the downstream side of the support, the z axis of the straw
  // volume is positive when going counter-clockwise.  This won't be important since we never work
  // in the straw-local coordiates within G4.
  CLHEP::HepRotationZ secRz(-sectorCenterPhi);
  CLHEP::HepRotationX secRx(M_PI/2.);
  G4RotationMatrix* sec00Rotation = _reg.add(G4RotationMatrix(secRx*secRz));

  // The z of the center of the placed sector envelope, in Mu2e coordinates.
  // This carries a sign, depending on upstream/downstream.
  double zSector(0.);
  for ( int i=0; i<sec00.nLayers(); ++i){
    zSector += sec00.getStraw(StrawId(0,0,i,0)).getMidPoint().z();
  }
  zSector /= sec00.nLayers();

  // A unit vector in the direction from the origin to the wire center within the sector envelope.
  CLHEP::Hep3Vector unit( cos(sectorCenterPhi), sin(sectorCenterPhi), 0.);

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

  // Place the straws into the sector envelope.
  for ( std::vector<Layer>::const_iterator i=sec00.getLayers().begin(); i != sec00.getLayers().end(); ++i ){

    Layer const& lay(*i);

    for ( std::vector<Straw const*>::const_iterator j=lay.getStraws().begin();
          j != lay.getStraws().end(); ++j ){

      Straw const&       straw(**j);
      StrawDetail const& detail(straw.getDetail());

      // To make graphical debugging less busy, create only a subset of the straws.
      if ( partialStraws ){
        if ( straw.id().getStraw()%8 != 0 && straw.id().getStraw() !=47 ) continue;
      }

      // Mid point of the straw in Mu2e coordinates.
      CLHEP::Hep3Vector const& pos(straw.getMidPoint());

      // Mid point of the straw, within the sector envelope.
      double r = (CLHEP::Hep3Vector( pos.x(), pos.y(), 0.)).mag();
      CLHEP::Hep3Vector mid = r*unit;
      mid.setZ(pos.z() - zSector);

      int copyNo=straw.index().asInt();
      bool edgeVisible(true);

      // The enclosing volume for the straw is made of gas.  The walls and the wire will be placed inside.
      VolumeInfo strawVol =  nestTubs( straw.name("TTrackerStrawGas_"),
                                       detail.getOuterTubsParams(),
                                       findMaterialOrThrow(detail.gasMaterialName()),
                                       sec00Rotation,
                                       mid,
                                       sec0Info.logical,
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

      if ( _verbosityLevel > 1 ){
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

    } // end loop over straws within a layer
  } // end loop over layers

  return sec0Info;

} // end preparePanel


// For debugging graphically, add three boxes that are aligned with the coordinate axes,
// on the tracker axis, just downstream of the tracker.
void
mu2e::ConstructTTrackerTDR::constructAxes(){

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


void
mu2e::ConstructTTrackerTDR::addPanels(VolumeInfo& basePanel, int idev, G4LogicalVolume* deviceLogical,
                                      double sectorCenterPhi ){

  Device const& dev = _ttracker.getDevice(idev);

  int secDraw = _config.getInt ("ttracker.secDraw",-1);

  // Assume all devices have the same number of sectors.
  int baseCopyNo = idev * _ttracker.getDevice(0).nSectors();

  // The sector will be placed in the middle of the channel.
  PlacedTubs const& chanUp(_ttracker.getSupportStructure().innerChannelUpstream());

  // Place the sector envelope into the device envelope, one placement per sector.
  for ( int isec=0; isec<dev.nSectors(); ++isec){

    // For debugging, only place the one requested sector.
    if ( secDraw > -1  && isec != secDraw ) continue;
    //if ( secDraw > -1  && isec%2 == 0 ) continue;

    // Choose a representative straw from this this (device,sector).
    Straw const& straw = _ttracker.getStraw( StrawId(idev,isec,0,0) );

    // Azimuth of the midpoint of the wire.
    CLHEP::Hep3Vector const& mid = straw.getMidPoint();
    double phimid = mid.phi();
    if ( phimid < 0 ) phimid += 2.*M_PI;

    // Is this sector on the front or back face of the plane?
    double sign   = ((mid.z() - dev.origin().z())>0. ? 1.0 : -1.0 );
    CLHEP::Hep3Vector sectorPosition(0.,0., sign*std::abs(chanUp.position().z()) );

    // The rotation that will make the straw mid point have the correct azimuth.
    // Sectors on the downstream side are flipped front/back by rotation about Y.
    // so the sense of the z rotation changes sign.
    double phi0 = sectorCenterPhi-phimid;
    if ( sign > 0 ){
      phi0 = sectorCenterPhi + M_PI +phimid;
    }

    CLHEP::HepRotationZ rotZ(phi0);
    CLHEP::HepRotationY rotY(M_PI);
    G4RotationMatrix* rotation  = (sign < 0 )?
      _reg.add(G4RotationMatrix(rotZ)) :
      _reg.add(G4RotationMatrix(rotZ*rotY));

    bool many(false);
    G4VPhysicalVolume* physVol = new G4PVPlacement(rotation,
                                                   sectorPosition,
                                                   basePanel.logical,
                                                   basePanel.name,
                                                   deviceLogical,
                                                   many,
                                                   baseCopyNo + isec,
                                                   false);
    if ( _doSurfaceCheck) {
      checkForOverlaps( physVol, _config, _verbosityLevel>0);
    }

  }

}

double
mu2e::ConstructTTrackerTDR::sectorHalfAzimuth(){

  TubsParams deviceEnvelopeParams = _ttracker.getDeviceEnvelopeParams();
  SupportStructure const& sup     = _ttracker.getSupportStructure();

  // Azimuth of the centerline of a sector enveleope: sectorCenterPhi
  // Upon construction and before placement, the sector envelope occupies [0,2.*sectorCenterPhi].
  double r1              = deviceEnvelopeParams.innerRadius();
  double r2              = sup.innerChannelUpstream().tubsParams().outerRadius();
  double halfChord       = sqrt ( r2*r2 - r1*r1 );
  double sectorCenterPhi = atan2(halfChord,r1);

  return  sectorCenterPhi;

}

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


// Create heirarchy of logical volumes for the elements of the support structure that live inside
// each device envelope.  The argument is passed in empty and is filled by this routine.
// It contains a VolumeInfo object for each logical volume at the top level of the heirarchy;
// these are placed inside a device envelope by the member function addPlaneSupports.
void
mu2e::ConstructTTrackerTDR::preparePlaneSupports( std::vector<VolumeInfo>& supportsToPlace ){

  bool supportVisible = _config.getBool("ttracker.supportVisible",false);
  bool supportSolid   = _config.getBool("ttracker.supportSolid",true);

  // Many parts of the support structure are G4Tubs objects.
  SupportStructure const& sup = _ttracker.getSupportStructure();
  G4Colour  lightBlue (0.0, 0.0, 0.75);
  std::vector<VHelper> vols;
  vols.reserve(11);

  vols.push_back( VHelper(sup.centerPlate(),         lightBlue,           "",                         supportVisible ));
  vols.push_back( VHelper(sup.outerRingUpstream(),   G4Colour::Green(),   "",                         supportVisible ));
  vols.push_back( VHelper(sup.outerRingDownstream(), G4Colour::Green(),   "",                         supportVisible ));
  vols.push_back( VHelper(sup.coverUpstream(),       G4Colour::Cyan(),    "",                         supportVisible ));
  vols.push_back( VHelper(sup.coverDownstream(),     G4Colour::Cyan(),    "",                         supportVisible ));
  vols.push_back( VHelper(sup.gasUpstream(),         G4Colour::Magenta(), "",                         supportVisible ));
  vols.push_back( VHelper(sup.gasDownstream(),       G4Colour::Magenta(), "",                         supportVisible ));
  vols.push_back( VHelper(sup.g10Upstream(),         G4Colour::Blue(),    sup.gasUpstream().name(),   supportVisible ));
  vols.push_back( VHelper(sup.g10Downstream(),       G4Colour::Blue(),    sup.gasDownstream().name(), supportVisible ));
  vols.push_back( VHelper(sup.cuUpstream(),          G4Colour::Yellow(),  sup.gasUpstream().name(),   supportVisible ));
  vols.push_back( VHelper(sup.cuDownstream(),        G4Colour::Yellow(),  sup.gasDownstream().name(), supportVisible ));

  for ( auto& vol : vols ){

    PlacedTubs const& part = vol.part;

    if ( vol.motherName.empty() ){

      // These volumes are top level volumes of the support structure and will be placed inside each device envelope.
      vol.info = nestTubs( part.name(),
                           part.tubsParams(),
                           findMaterialOrThrow(part.materialName()),
                           noRotation,
                           zeroVector,
                           nullptr,
                           0,
                           vol.visible,
                           vol.colour,
                           supportSolid,
                           _forceAuxEdgeVisible,
                           doNotPlace,
                           _doSurfaceCheck
                           );

      vol.info.centerInParent = part.position();

      supportsToPlace.push_back(vol.info);

    } else{

      // These volumes are lower level volumes and are placed, at this time, inside their parents.
      G4LogicalVolume* motherLogical = findMotherLogical( vols, vol.motherName );
      vol.info = nestTubs( part.name(),
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

  }

  {

    PlacedTubs const& ring     = sup.innerRing();
    PlacedTubs const& chanUp   = sup.innerChannelUpstream();
    PlacedTubs const& chanDown = sup.innerChannelDownstream();

    // Parameters of the polycone.
    const int n(10);
    double z[n];
    double rIn[n];
    double rOut[n];

    // I am not sure if I am allowed to have two coincident zplanes so slightly
    // slope the edges of the channel.
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

    VolumeInfo innerRingInfo( ring.name(), ring.position(), zeroVector );
    innerRingInfo.solid = new G4Polycone( ring.name(),
                                          0.,
                                          2.*M_PI,
                                          n,
                                          z,
                                          rIn,
                                          rOut);
    /*
    // I would like to use this c'tor of Polycone but when I do so, ROOT cannot work with the GDML file that it creates.
    const int n(12);
    double z[n];
    double r[n];
    z[ 0] = -ring.tubsParams().zHalfLength();                              r[ 0] = ring.tubsParams().outerRadius();
    z[ 1] = z[0];                                                          r[ 1] = ring.tubsParams().innerRadius();
    z[ 2] = chanUp.position().z()   - chanUp.tubsParams().zHalfLength();   r[ 2] = r[1];
    z[ 3] = z[2];                                                          r[ 3] = chanUp.tubsParams().outerRadius();
    z[ 4] = chanUp.position().z()   + chanUp.tubsParams().zHalfLength();   r[ 4] = r[3];
    z[ 5] = z[4];                                                          r[ 5] = ring.tubsParams().innerRadius();
    z[ 6] = chanDown.position().z() - chanDown.tubsParams().zHalfLength(); r[ 6] = r[5];
    z[ 7] = z[6];                                                          r[ 7] = chanDown.tubsParams().outerRadius();
    z[ 8] = chanDown.position().z() + chanDown.tubsParams().zHalfLength(); r[ 8] = r[7];
    z[ 9] = z[8];                                                          r[ 9] = ring.tubsParams().innerRadius();
    z[10] = ring.tubsParams().zHalfLength();                               r[10] = r[9];
    z[11] = z[10];                                                         r[11] = r[0];

    for ( int i=0; i<n; ++i ){
      cout << "rz: " << i << " " << r[i] << " " << z[i] << endl;
    }

    VolumeInfo innerRingInfo( ring.name(), ring.position(), zeroVector );
    G4Polycone* pcon = new G4Polycone( ring.name(),
                                       0.,
                                       2.*M_PI,
                                       n,
                                       r,
                                       z);
    innerRingInfo.solid = pcon;
    pcon->StreamInfo(cout);
    */

    finishNesting( innerRingInfo,
                   findMaterialOrThrow(ring.materialName()),
                   noRotation,
                   zeroVector,
                   nullptr,
                   0,
                   supportVisible,
                   G4Colour::Red(),
                   supportSolid,
                   _forceAuxEdgeVisible,
                   doNotPlace,
                   false
                   );

    innerRingInfo.centerInParent = ring.position();

    supportsToPlace.push_back(innerRingInfo);
  }


} // end preparePlaneSupports

void mu2e::ConstructTTrackerTDR::addPlaneSupports( std::vector<VolumeInfo>& supportsInfo, int idev, VolumeInfo const& devInfo ){

  for ( auto& info : supportsInfo ){

    if ( _verbosityLevel > 0 ) {
      cout << "Plane Support: "
           << info.name      << " "
           << info.solid     << " "
           << info.logical   << " "
           << info.physical  << " "
           << info.centerInParent << " "
           << info.centerInWorld  << " "
           << endl;
    }

    bool many(false);

    info.physical = new G4PVPlacement(noRotation,
                                      info.centerInParent,
                                      info.logical,
                                      info.name,
                                      devInfo.logical,
                                      many,
                                      idev,
                                      _doSurfaceCheck);

  }

} // end addPlaneSupports
