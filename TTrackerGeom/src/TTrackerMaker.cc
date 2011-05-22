//
// Construct and return an TTracker.
//
//
// $Id: TTrackerMaker.cc,v 1.32 2011/05/22 19:09:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/22 19:09:16 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/for_all.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TTrackerGeom/inc/TTrackerMaker.hh"
#include "TrackerGeom/inc/Sector.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iomanip>
#include <iostream>

using cet::square;

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  TTrackerMaker::TTrackerMaker( SimpleConfig const& config){
    parseConfig(config);
    buildIt( );
  }

  void TTrackerMaker::parseConfig( const SimpleConfig& config ){

    _numDevices         = config.getInt("ttracker.numDevices");
    _sectorsPerDevice   = config.getInt("ttracker.sectorsPerDevice");
    _layersPerSector    = config.getInt("ttracker.layersPerSector");
    _manifoldsPerEnd    = config.getInt("ttracker.manifoldsPerEnd");
    _strawsPerManifold  = config.getInt("ttracker.strawsPerManifold");
    _rotationPattern    = config.getInt("ttracker.rotationPattern");
    _spacingPattern     = config.getInt("ttracker.spacingPattern");

    _zCenter              = config.getDouble("ttracker.z0")*CLHEP::mm;
    _envelopeInnerRadius  = config.getDouble("ttracker.envelopeInnerRadius")*CLHEP::mm;
    _strawOuterRadius     = config.getDouble("ttracker.strawOuterRadius")*CLHEP::mm;
    _strawWallThickness   = config.getDouble("ttracker.strawWallThickness")*CLHEP::mm;
    _strawGap             = config.getDouble("ttracker.strawGap")*CLHEP::mm;
    _deviceSpacing        = config.getDouble("ttracker.deviceSpacing")*CLHEP::mm;
    _deviceHalfSeparation = config.getDouble("ttracker.deviceHalfSeparation")*CLHEP::mm;
    _deviceRotation       = config.getDouble("ttracker.deviceRotation")*CLHEP::degree;

    _outerSupportRadius   = config.getDouble("ttracker.outerSupportRadius")*CLHEP::mm;
    _innerSupportRadius   = config.getDouble("ttracker.innerSupportRadius")*CLHEP::mm;
    _supportHalfThickness = config.getDouble("ttracker.supportHalfThickness")*CLHEP::mm;

    _wireRadius           = config.getDouble("ttracker.wireRadius")*CLHEP::mm;

    _manifoldYOffset      = config.getDouble("ttracker.manifoldYOffset")*CLHEP::mm;
    config.getVectorDouble("ttracker.manifoldHalfLengths", _manifoldHalfLengths, 3);
    for ( size_t i=0; i<_manifoldHalfLengths.size(); ++i ){
      _manifoldHalfLengths.at(i) *= CLHEP::mm;
    }

    config.getVectorString("ttracker.strawMaterials", _strawMaterials, 3);

    _envelopeMaterial = config.getString("ttracker.mat.vacuum");
    _supportMaterial = config.getString("ttracker.mat.support");

    //string ttracker.mat.manifold  = "G4_Al";  // Placeholder.

    // Also define some parameters that may become variable some day.
    _sectorBaseRotations.clear();
    _sectorZSide.clear();

    if (_rotationPattern == 0 && _sectorsPerDevice == 4 ){
      _sectorBaseRotations.push_back(   0.*CLHEP::degree);
      _sectorBaseRotations.push_back(  90.*CLHEP::degree);
      _sectorBaseRotations.push_back( 180.*CLHEP::degree);
      _sectorBaseRotations.push_back( 270.*CLHEP::degree);

      _sectorZSide.push_back(-1.);
      _sectorZSide.push_back(+1.);
      _sectorZSide.push_back(-1.);
      _sectorZSide.push_back(+1.);
    }
    else if (_rotationPattern == 1 && _sectorsPerDevice == 6 ){
      _sectorBaseRotations.push_back(   0.*CLHEP::degree);
      _sectorBaseRotations.push_back(  60.*CLHEP::degree);
      _sectorBaseRotations.push_back( 120.*CLHEP::degree);
      _sectorBaseRotations.push_back( 180.*CLHEP::degree);
      _sectorBaseRotations.push_back( 240.*CLHEP::degree);
      _sectorBaseRotations.push_back( 300.*CLHEP::degree);

      _sectorZSide.push_back(-1.);
      _sectorZSide.push_back(+1.);
      _sectorZSide.push_back(-1.);
      _sectorZSide.push_back(+1.);
      _sectorZSide.push_back(-1.);
      _sectorZSide.push_back(+1.);
    }

    else {
      throw cet::exception("GEOM")
        << "Unrecognized rotation pattern in TTrackerMaker. \n";
    }

    // Parts of the algorithm require that the Aseet style tracker is built
    // of stations with 2 devices each.
    if ( _spacingPattern == 1 ){
      if ( _numDevices%2 == 1 ){
        throw cet::exception("GEOM")
          << "Aseet style tracker requires 2 devices per station.\n"
          << "So ttracker.numDevices must be even.  It was: "
          << _numDevices
          << "\n";
      }
    }

    _sectorBoxXOffsetMag = 0.0;
    _sectorBoxZOffsetMag = 0.0;
    _layerHalfSpacing = 0.0;
    _layerHalfShift = 0.0;
    _manifoldXEdgeExcessSpace = 0.0;
    _manifoldZEdgeExcessSpace = 0.0;

  } // end TTrackerMaker::parseConfig

  void lptest( const Layer& lay){
    cout << lay.id() << " |  "
         << lay.nStraws()  <<  " |  "
         << lay.getStraws().capacity() << " "
         << endl;
  }

  void devtest( const Device& dev){
    cout << "Device: "
         << dev.id() << " "
         << dev.origin() << " "
         << dev.rotation()
         << endl;
  }

  void positionTest( const Layer& lay){
    const Straw& straw = lay.getStraw( 0 );
    cout << "Layer: "
         << lay.id() << " "
         << straw.getMidPoint().z() <<  " "
         << straw.getDirection().z() << " "
         << endl;
  }

  void TTrackerMaker::buildIt(){

    // Make an empty TTracker.
    _tt = auto_ptr<TTracker>(new TTracker());

    // Fill the information about the supports.
    _tt->_supportParams = Support( _innerSupportRadius,
                                   _outerSupportRadius,
                                   _supportHalfThickness,
                                   _supportMaterial);

    _tt->_z0                  = _zCenter;
    _tt->_envelopeInnerRadius = _envelopeInnerRadius;
    _tt->_manifoldHalfLengths = _manifoldHalfLengths;
    _tt->_envelopeMaterial    = _envelopeMaterial;

    // Z location of the first device.
    _z0 = -findFirstDevZ0();
    //    cout << "First device z0: " << _z0 << endl;

    makeDetails();

    // Reserve space for straws so that pointers are valid.
    _nStrawsToReserve = _numDevices * _sectorsPerDevice * _layersPerSector *
      _manifoldsPerEnd * _strawsPerManifold;
    //_tt->_allStraws.reserve(_nStrawsToReserve); // see makeLayer

    _tt->_devices.reserve(_numDevices);
    // Construct the devices and their internals.
    for ( int idev=0; idev<_numDevices; ++idev ){
      makeDevice( DeviceId(idev) );
    }

    // Fill all of the non-persistent information.
    _tt->fillPointers();

    identifyNeighbourStraws();

    //_tt->forAllLayers( lptest);
    //_tt->forAllDevices( devtest);

    //_tt->forAllLayers( positionTest);

  } //end TTrackerMaker::buildIt.

  void TTrackerMaker::makeDevice( DeviceId devId ){

    int idev = devId;

    double devDeltaZ = chooseDeviceSpacing(idev);
    CLHEP::Hep3Vector origin( 0., 0., _z0+devDeltaZ);
    //    cout << "Device z: " << origin.z() << endl;

    //    double phi = chooseDeviceRotation(idev)+_deviceRotation;
    double phi = chooseDeviceRotation(idev);

    _tt->_devices.push_back(Device(devId, origin, phi));
    Device& dev = _tt->_devices.back();
    dev._sectors.reserve(_sectorsPerDevice);

    for ( int isec=0; isec<_sectorsPerDevice; ++isec ){
      makeSector ( SectorId(devId,isec), dev );
    }

  }

  void TTrackerMaker::makeSector( const SectorId& secId, Device& dev ){

    dev._sectors.push_back( Sector(secId) );
    Sector& sector = dev._sectors.back();
    sector._layers.reserve(_layersPerSector);

    // check if the opposite sectors do not overlap
    static double const tolerance = 1.e-6; // this should be in a config file

    if ((2.*_manifoldHalfLengths.at(2)+_supportHalfThickness)>_deviceHalfSeparation + tolerance) {
      cout << "(2.*_manifoldHalfLengths.at(2)+_supportHalfThickness), _deviceHalfSeparation " <<
        (2.*_manifoldHalfLengths.at(2)+_supportHalfThickness) << ", " <<_deviceHalfSeparation << endl;
      throw cet::exception("GEOM")  << "Devices are to close \n";
    }

    makeManifolds( secId );

    double strawSpacing = _strawGap+2.*_strawOuterRadius;

    for ( int ilay=0; ilay<_layersPerSector; ++ilay ){
      makeLayer( LayerId(secId,ilay), sector );

      // checking spacing of the individual layers
      // are the manifolds sized correctly for the straws?

      Layer const & layer = sector.getLayer(ilay);
      //      cout << "Debugging looking at the layer   : " << layer.id() << endl;
      for (int ns = 0; ns!=layer.nStraws()-1; ++ns) {
        double layerDeltaMag =
          (layer.getStraw(ns+1).getMidPoint() - layer.getStraw(ns).getMidPoint()).mag();
        if ( abs(layerDeltaMag-strawSpacing)> tolerance ) {
          cout << "Layer straw spacing is (mm)   : " << layerDeltaMag
               << " for layer " << layer.id() << " straw " << layer.getStraw(ns).id()  << endl;
          cout << "It should be                  : " << strawSpacing << " diff: "
               << (layerDeltaMag-strawSpacing) << endl;

          throw cet::exception("GEOM")  << "Incorrect intralayer straw spacing, check manifold sizes rtc..\n";

        }
      }
    }

    // check spacing between layers/straws

    if (_layersPerSector>1) {

      Layer const & layer0 = sector.getLayer(0);
      Layer const & layer1 = sector.getLayer(1);
      // cout << "Debugging looking at the layers   : " << layer0.id() << ", " << layer1.id() << endl;

      for (int ns = 0; ns!=layer0.nStraws(); ++ns) {
        double xLayerDeltaMag =
          (layer0.getStraw(ns).getMidPoint() - layer1.getStraw(ns).getMidPoint()).mag();
        if ( abs(xLayerDeltaMag-strawSpacing)> tolerance ) {
          cout << "xLayer straw spacing is (mm)   : "
               << xLayerDeltaMag
               << " for straws: "
               << layer0.getStraw(ns).id() << ", " << layer1.getStraw(ns).id()
               << endl;
          cout << "It should be                   : "
               << strawSpacing << " diff: "
               << (xLayerDeltaMag-strawSpacing) << endl;

          throw cet::exception("GEOM")  << "Incorrect interlayer straw spacing \n";

        }
      }


      for (int ns = 1; ns!=layer0.nStraws(); ++ns) {
        double xLayerDeltaMag =
          (layer0.getStraw(ns).getMidPoint() - layer1.getStraw(ns-1).getMidPoint()).mag();
        if ( abs(xLayerDeltaMag-strawSpacing)> tolerance ) {
          cout << "xLayer straw spacing is (mm)   : "
               << xLayerDeltaMag
               << " for straws: "
               << layer0.getStraw(ns).id() << ", " << layer1.getStraw(ns-1).id()
               << endl;
          cout << "It should be                   : "
               << strawSpacing << " diff: "
               << (xLayerDeltaMag-strawSpacing) << endl;

          throw cet::exception("GEOM")  << "Incorrect interlayer straw spacing \n";

        }
      }
    }

    // calculate/make a sector envelope
    computeSectorBoxParams(sector, dev);

  }

  void TTrackerMaker::makeLayer ( const LayerId& layId, Sector& sector ){

    // Make an empty layer object.
    sector._layers.push_back( Layer(layId) );
    Layer& layer = sector._layers.back();

    // Get additional bookkeeping info.
    deque<Straw>& allStraws = _tt->_allStraws;
    const Device& device = _tt->getDevice( layId.getDeviceId() );
    int ilay = layId.getLayer();
    int isec = layId.getSector();

    //    cout << "Debugging TTrackerMaker ilay: " << ilay << endl;

    // Start to populate the layer.
    layer._nStraws = _manifoldsPerEnd*_strawsPerManifold;
    layer._straws.reserve(_manifoldsPerEnd*_strawsPerManifold);

    // |z| of straw center, relative to the center of the device.
    // Sign is taken care of elsewhere.

    // see similar calc in computeSectorBoxParams
    //    double zOffset = _supportHalfThickness + _strawOuterRadius + ilay*2.*_layerHalfSpacing;
    // the above commented out calculation places the straws at the edge of the manifold in Z

    double zOffset = _supportHalfThickness + _manifoldZEdgeExcessSpace +
      _strawOuterRadius + ilay*2.*_layerHalfSpacing;

    // Rotation that puts wire direction and wire mid-point into their
    // correct orientations.
    // CLHEP::HepRotationZ RZ(_sectorBaseRotations.at(isec));
    CLHEP::HepRotationZ RZ(_sectorBaseRotations.at(isec) + device.rotation());

    // Unit vector in the wire direction. (nominal is the sector 0 to the right?)
    CLHEP::Hep3Vector unit = RZ*CLHEP::Hep3Vector(0.,1.,0.);

    // Straw number within the layer; does not reset to zero at each manifold.
    int _istraw(-1);

    // Add all of the straws
    for ( int iman=0; iman<_manifoldsPerEnd; ++iman ){

      // Inner edge of the innermost wire connected to this manifold.
      // this is layer dependent, take care of it at xstraw below

      double xA = _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*iman + _manifoldXEdgeExcessSpace;

      // Add all of the straws connected to this manifold.
      for ( int istr=0; istr<_strawsPerManifold; ++istr ){
        ++_istraw;

        // layers with fewer straws would complicate StrawSD, constructTTrackerv, TTrackerMaker

        // Construct straw midpoint in its base position in the
        // coord system of the device envelope.
        // we will shift the "second" layer from the manifold edge

        double xstraw = (ilay%2==0) ?
          xA + (1 + 2*istr)*_strawOuterRadius + istr*_strawGap :
          xA + (1 + 2*istr)*_strawOuterRadius + istr*_strawGap + 2.0*_layerHalfShift;

        CLHEP::Hep3Vector mid( xstraw, 0., zOffset*_sectorZSide.at(isec) );
        mid += device.origin();

        // Rotate straw midpoint to its actual location.
        CLHEP::Hep3Vector offset = RZ*mid;

        StrawIndex index(allStraws.size());

        allStraws.push_back( Straw( StrawId( layId, _istraw),
                                    index,
                                    offset,
                                    &_tt->_strawDetails.at(iman),
                                    iman,
                                    unit
                                    )
                             );
        layer._straws.push_back(&allStraws.back());
        layer._indices.push_back(index);


//         if ( layId.getDevice() != -1 ){
//           cout << "Position: " << setw(3) <<
//             layId << " | " << setw(3) <<
//             iman << " " << setw(3) <<
//             istr                << " | " << setw(3) <<
//             _istraw << " " << fixed << setprecision(2) << setw(8) <<
//             xstraw << " " << fixed << setprecision(2) << setw(8) <<
//             2.*_strawHalfLengths.at(iman) << " " << fixed << setprecision(2) <<
//             mid << " " << fixed << setprecision(2) << setw(8) <<
//             device.origin() << " | "  << setw(3) <<
//             index <<  " " << setw(3) <<
//             allStraws.size() << " "  << setw(3) <<
//             layer._straws.size() << " | " << setw(5) <<
//             (allStraws.back()).id() << ", " << setw(5) <<
//             (allStraws.back()).index()
//           << endl;
//           }

      }
    }

  }

  void TTrackerMaker::makeManifolds( const SectorId& secId){

    if ( _sectorsPerDevice != 4 && _sectorsPerDevice != 6 ) {
      throw cet::exception("GEOM")
        << "This code only knows how to do 4 or 6 sectors per device.\n";
    }

    double phi = _tt->getDevice(secId.getDevice()).rotation();
    CLHEP::HepRotationZ RZ(phi);


    // manifold objects are not used for now...

    for ( int i=0; i<_manifoldsPerEnd; ++i){

      // First compute everything in their nominal positions: sector 0, right ?
      double x0 = _envelopeInnerRadius +
        _strawsPerManifold*_strawOuterRadius +
        _manifoldHalfLengths.at(0);

      double y0 = _tt->_strawDetails.at(i).halfLength() + _manifoldHalfLengths.at(2);

      double z0 = ( _supportHalfThickness + _manifoldHalfLengths.at(2) );
      //if ( secId.getSector() <= 1 ) z0 = -z0; // is this correct for the 6 sectors?
      // why not *_sectorZSide.at(secId.getSector())

      z0 = z0*_sectorZSide.at(secId.getSector());

      // is the above assuming correct Z? why not secId.getSector()%2?
      // are manifolds ever used?
      // is it correct at all? I mean the origin? It is always the same for each device...
      // x never changes

      CLHEP::Hep3Vector origin(x0,y0,z0);

//       cout << "Manifold device, sector, origin, length[0] :" <<
//         _tt->getDevice(secId.getDevice()).id() << ", " <<
//         secId.getSector() << ", " <<
//         origin << ", " << _manifoldHalfLengths.at(0) <<endl;


      _tt->_allManifolds.push_back( Manifold( origin, _manifoldHalfLengths) );
    }

  }

  // Assumes all devices and all sectors are the same.
  // The straw length depends only on the manifold number.
  // See Mu2e-doc-??? for the algorithm.
  void TTrackerMaker::computeStrawHalfLengths(){

    computeManifoldEdgeExcessSpace();
    _strawHalfLengths.clear();

    for ( int i=0; i<_manifoldsPerEnd; ++i ){

      double xA =
        _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*i + _manifoldXEdgeExcessSpace;

//    double xA = (_layersPerSector==1) ?
//      _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*i + _manifoldXEdgeExcessSpace :
//      _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*i + _manifoldXEdgeExcessSpace +
//       _strawOuterRadius;

      // we ignore the further laying straws in the multi layer case,
      // as this would make the straws shorter than they need to be
      // the wire positioning is not affected by this though

      double yA = sqrt( square(_innerSupportRadius) - square(xA) );
      double yB = yA + _manifoldYOffset;
      // cout << "Straw Length: " << xA << " " << yB*2.0 << endl;
      _strawHalfLengths.push_back(yB);
    }
  }

  void TTrackerMaker::makeDetails(){

    computeConstantSectorBoxParams();

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      _tt->_strawDetails.push_back
        ( StrawDetail
          ( i,
            _strawMaterials,
            _strawOuterRadius,
            _strawWallThickness,
            _strawHalfLengths.at(i),
            _wireRadius
            )
          );
    }

  }

  void  TTrackerMaker::computeSectorBoxParams(Sector& sector, Device& dev){

    // get sector number
    int isec = sector.id().getSector();
    //    int idev = dev.id();
    //    cout << "Debugging TTrackerMaker isec,idev: " << isec << ", " << idev << endl;

    // we copy precalculated _sectorBoxHalfLengths etc.. into the sector

    sector._boxHalfLengths = _sectorBoxHalfLengths;
    double xOffset = _sectorBoxXOffsetMag;
    double zOffset = _sectorBoxZOffsetMag * _sectorZSide.at(isec);

    // Now calculate the rotations and placement of the sector envelope

    sector._boxRxAngle = 0.;
    sector._boxRyAngle = 0.;
    sector._boxRzAngle = _sectorBaseRotations.at(isec) + dev.rotation();

    CLHEP::HepRotationZ RZ(sector._boxRzAngle);

    sector._boxOffset  = RZ*(CLHEP::Hep3Vector( xOffset, 0., zOffset) + dev.origin());

    // we set to 0.0  values smaller than a small number
    static const double max0val = 1.e-06;
    for (int ii=0; ii!=sector._boxOffset.SIZE; ++ii) {
      if (abs(sector._boxOffset[ii])<max0val) sector._boxOffset[ii]=0.0;
    }

//     cout << "Debugging sector box isec, by, bz, bxl, bxs, boxRzAngle, boxOffset: " <<
//       isec << ", " <<
//       by << ", " <<
//       bz << ", " <<
//       bxl << ", " <<
//       bxs << ", " <<
//       sector._boxRzAngle << ", " <<
//       sector._boxOffset <<
//       endl;

//     cout << "Debugging sector box isec, straw lengths: ";
//     for ( int i=0; i<_manifoldsPerEnd; ++i ){
//       cout << i << " " << _strawHalfLengths.at(i);
//     }
//     cout << endl;

    return;

  }

  void TTrackerMaker::computeLayerSpacingAndShift(){

    _layerHalfSpacing = (_layersPerSector<=1) ? 0.0 :
      sqrt(3.0*(square(_strawOuterRadius)+_strawOuterRadius*_strawGap+0.25*square(_strawGap)))*0.5;
    _layerHalfShift   = (_layersPerSector<=1) ? 0.0 : _strawGap*0.25 + _strawOuterRadius*0.5;

  }

  void TTrackerMaker::computeManifoldEdgeExcessSpace(){

    computeLayerSpacingAndShift();

    // Computes space between first/last straw and edge of manifold

    _manifoldXEdgeExcessSpace = _manifoldHalfLengths.at(0) -
      _strawOuterRadius*_strawsPerManifold -
      _strawGap*(_strawsPerManifold-1)*0.5;

    _manifoldZEdgeExcessSpace = _manifoldHalfLengths.at(2) - _strawOuterRadius -
      (_layersPerSector-1)*_layerHalfSpacing;

//     cout << "Debugging,  _manifoldXEdgeExcessSpace, _manifoldZEdgeExcessSpace: " <<
//       _manifoldXEdgeExcessSpace << ", " << _manifoldZEdgeExcessSpace << endl;

    if ( _manifoldXEdgeExcessSpace < 0.0 || _manifoldZEdgeExcessSpace < 0.0){
      throw cet::exception("GEOM")
        << "Manifolds are too small to hold straws!\n";
    }

  }

  void TTrackerMaker::computeConstantSectorBoxParams() {

    computeStrawHalfLengths();

    // the box is a trapezoid ;

    // note that G4 has it own coordinate convention for each solid
    // (see mu2e<->G4 translation below):

    // trapezoid y dimensions (or x in G4)
    // trapezoid x dimension  (or z in G4)
    // trapezoid z dimension  (or y in G4)

    // shorter y             is the length of the straws in the top/last (shortest) manifold
    // longer  y is longer than the length of the straws in the longest manifold
    // the other dimentions are "x", the combined manifold width + _manifoldXEdgeExcessSpace
    // the "thickness" of the trpezoid z, ~ the layer thickness * number of layers

    // x

    // if there are more than one layer per sector, manifolds have a
    // straw "sticking out" by a straw radius (and half of the straw
    // gap), but only the last one extends beyond the entire structure
    // in the z direction; remember that those are "halfLengths", so
    // the sticking out straw radius has to be divided by 2 and the gap by 4

    // _layerHalfShift is 0 in one layer case
    double bx = _manifoldHalfLengths.at(0)*double(_manifoldsPerEnd) + _layerHalfShift;

    // calculating "longer" x;
    // starting from the tng of the slope

    // calculate the largest slope starting from the longest straws manifold
    double maxtg = 0.0;

    // the code below looks at the slope "seen" from the longest set of straws
    for (int i=1; i!=_manifoldsPerEnd; ++i) {

      double ttg = ( _manifoldHalfLengths.at(0) + _layerHalfShift )*double(i) /
        ( _strawHalfLengths.at(0) - _strawHalfLengths.at(i) ) ;

      if (maxtg < ttg ) {
        maxtg = ttg;
      }

    }

    // finally y (long, short)
    double byl = (_manifoldHalfLengths.at(0) + _layerHalfShift)/maxtg + _strawHalfLengths.at(0);

    double bys = byl - bx/maxtg;

    // z
    // manifold better be thicker than the straws:

    double bz = _manifoldHalfLengths.at(2);

    // now push it all back into the vector
    // std::vector<double> _sectorBoxHalfLengths;

    static const size_t sectorBoxHalfLengthsSize= 5;

    _sectorBoxHalfLengths.reserve(sectorBoxHalfLengthsSize);

    // Pad the trapezoid to be slightly larger than it needs to be
    static const double pad = 0.0; // this needs to be in the geom file... if to be non-zero

    // the order is forced by the LTracker and nestTrp/G4Trd and Sector data

    _sectorBoxHalfLengths.push_back(0.0); //dummy to be compatible with LTracker
    _sectorBoxHalfLengths.push_back(bx+pad);
    _sectorBoxHalfLengths.push_back(bz+pad);

    _sectorBoxHalfLengths.push_back(bys+pad);
    _sectorBoxHalfLengths.push_back(byl+pad);

    if (_sectorBoxHalfLengths.size()!=sectorBoxHalfLengthsSize) {
      cout << " _sectorBoxHalfLengths.size() sould be " << sectorBoxHalfLengthsSize <<
        ", but is : " << _sectorBoxHalfLengths.size() << endl;
      throw cet::exception("GEOM")
        << "something is wrong with sector _sectorBoxHalfLengths calculations \n";
    }

    _sectorBoxXOffsetMag = _envelopeInnerRadius  + _sectorBoxHalfLengths.at(1); // bx + pad
    _sectorBoxZOffsetMag = _supportHalfThickness + _sectorBoxHalfLengths.at(2); // bz + pad

    // we need to make sure the trapezoids do not extend beyond the device envelope...
    // we will check if the dev envelope radius acomodates the newly created box

    double outerSupportRadiusRequireds =
      sqrt(square(_envelopeInnerRadius + 2.0*_sectorBoxHalfLengths.at(1))+
           square(_sectorBoxHalfLengths.at(3)));
    double outerSupportRadiusRequiredl =
      max ( outerSupportRadiusRequireds,
            sqrt(square(_envelopeInnerRadius-pad)+ square(_sectorBoxHalfLengths.at(4))) );

//     if (true) {
//       cout << "Debugging _strawHalfLengths: ";
//       for (size_t i=0; i!=_manifoldsPerEnd; ++i) {
//         cout << _strawHalfLengths.at(i)  << ", ";
//       }
//       cout << endl;
//       cout << "Debugging _supportParams.innerRadius   :   " << _tt->_supportParams.innerRadius() << endl;
//       cout << "Debugging _supportParams.outerRadius   :   " << _tt->_supportParams.outerRadius() << endl;
//       cout << "Debugging _supportParams.outerRadius rs:   " << outerSupportRadiusRequireds << endl;
//       cout << "Debugging _supportParams.outerRadius rl:   " << outerSupportRadiusRequiredl << endl;
//     }

    if (_tt->_supportParams.outerRadius() < outerSupportRadiusRequiredl) {
      cout << " _supportParams.outerRadius         :   " << _tt->_supportParams.outerRadius() << endl;
      cout << " _supportParams.outerRadius required:   " << outerSupportRadiusRequiredl << endl;
      throw cet::exception("GEOM")
        << "outerSupportRadius is to small given other paramters \n";
    }

    return;

  }


  // Compute the rotation for the given device.
  double TTrackerMaker::chooseDeviceRotation( int idev ) const{


    if ( _rotationPattern == 0 ){
      int k = idev%3;
      if ( k == 0 ) {
        return 0.;
      } else if ( k == 1 ){
        return _deviceRotation;
      } else{
        return -_deviceRotation;
      }
    }
    else if ( _rotationPattern == 1 ){
      int k = idev%2;
      if ( k == 0 ) {
        return 0.;
      } else {
        return _deviceRotation/2.;
      }
    }
    else {
      throw cet::exception("GEOM")
        << "Unrecognized rotation pattern in TTrackerMaker. \n";
    }

  } //end TTrackerMaker::chooseDeviceRotation


  // Compute the spacing for the given device.
  double TTrackerMaker::chooseDeviceSpacing( int idev ) const{


    if ( _spacingPattern == 0 ) {
      return idev * _deviceSpacing;
    }

    else if ( _spacingPattern == 1 ) {
      int station = idev/2;
      int k = idev%2;
      if (k == 0 ) {
        return  station*_deviceSpacing - _deviceHalfSeparation;
      } else if (k == 1 ) {
        return  station*_deviceSpacing + _deviceHalfSeparation;
      }
    }

    throw cet::exception("GEOM")
      << "Unrecognized separation pattern in TTrackerMaker. \n";

  }

  double TTrackerMaker::findFirstDevZ0() const{

    if ( _spacingPattern == 0 ) {
      return _deviceSpacing*double(_numDevices-1)/2.0;
    }

    else if ( _spacingPattern == 1 ) {
      int nStations = _numDevices/2;
      return double(nStations-1)/2.0 * _deviceSpacing;
    }
    else {
      throw cet::exception("GEOM")
        << "Unrecognized separation pattern in TTrackerMaker. \n";
    }
  }


  // Identify the neighbour straws for all straws in the tracker
  void TTrackerMaker::identifyNeighbourStraws() {

    deque<Straw>& allStraws = _tt->_allStraws;

    for (deque<Straw>::iterator i = allStraws.begin();
         i != allStraws.end(); ++i) {
      // throw exception if more than 2 layers per sector

      if (_tt->getSector(i->id().getSectorId()).nLayers() > 2 ) {
        throw cet::exception("GEOM")
          << "The code works with no more than 2 layers per sector. \n";
      }

      LayerId lId = i->id().getLayerId();
      int layer = lId.getLayer();
      int nStrawLayer = _tt->getLayer(lId)._nStraws;

      //  cout << lId << " has " << nStrawLayer << " straws" << endl;
      //  cout << "Analyzed straw: " << i->id() << '\t' << i->index() << endl;

      // add the "same layer" n-1 neighbours straw (if exist)

      if ( i->id().getStraw() ) {
        const StrawId nsId(lId, (i->id().getStraw())-1 );
        i->_nearestById.push_back( nsId );
        // Straw temp = _tt->getStraw( nsId );
        // cout << "Neighbour left straw: " << temp.id() << '\t' << temp.index() << endl;
        i->_nearestByIndex.push_back( _tt->getStraw(nsId).index() );
      }

      // add the "same layer" n+1 neighbours straw (if exist)

      if ( i->id().getStraw() < (nStrawLayer-1) ) {
        const StrawId nsId(lId, (i->id().getStraw())+ 1 );
        i->_nearestById.push_back( nsId );
        // Straw temp = _tt->getStraw( nsId );
        // cout << "Neighbour right straw: " << temp.id() << '\t' << temp.index() << endl;
        i->_nearestByIndex.push_back( _tt->getStraw(nsId).index() );
      }

      // add the "opposite layer" n neighbours straw (if more than 1 layer)

      if (_layersPerSector == 2) {
        const StrawId nsId( i->id().getSectorId(), (layer+1)%2, (i->id().getStraw()) );

        // throw exception if the two layer of the same sector have different
        // number of straws
        if (_tt->getLayer(lId)._nStraws != nStrawLayer) {
          throw cet::exception("GEOM")
            << "The code works only with the same number of straws "
            << "per layer in the same sector. \n";
        }

        i->_nearestById.push_back( nsId );
        // Straw temp = _tt->getStraw( nsId );
        // cout << "Neighbour opposite straw: " << temp.id() << '\t' << temp.index() << endl;
        i->_nearestByIndex.push_back( _tt->getStraw( nsId ).index() );

        // add the "opposite layer" n+-1 neighbours straw (if exist)

        if ( i->id().getStraw() > 0 && i->id().getStraw() < nStrawLayer-1 ) {
          const StrawId nsId( i->id().getSectorId(), (layer+1)%2,
                              (i->id().getStraw()) + (layer?1:-1));
          i->_nearestById.push_back( nsId );
          // Straw temp = _tt->getStraw( nsId );
          // cout << "Neighbour opposite +- 1 straw: " << temp.id() << '\t' << temp.index() << endl;
          i->_nearestByIndex.push_back( _tt->getStraw( nsId ).index() );
        }

      }
    }
  }


} // namespace mu2e
