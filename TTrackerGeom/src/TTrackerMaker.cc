//
// Construct and return a TTracker.
//
//
// $Id: TTrackerMaker.cc,v 1.57 2014/04/11 04:40:46 genser Exp $
// $Author: genser $
// $Date: 2014/04/11 04:40:46 $
//
// Original author Rob Kutschke
//
// Significant change mf 5/30/11:
//    Insertion into the Straw objects data useful for early (hit/miss
//    based) pattern recognition.  See identifyDirectionalNeighborStraws().

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "ConfigTools/inc/SimpleConfig.hh"
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
using cet::diff_of_squares;
using cet::sum_of_squares;

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  TTrackerMaker::TTrackerMaker( SimpleConfig const& config){
    parseConfig(config);
    buildIt( );

    // print straw layout for debbugging pourposes

    if (_verbosityLevel>2) {


      int idev = -1;
      int isec = -1;
      int ilay = -1;
      double iang = -36000;

      size_t nstraws = (_tt->_allStraws).size();
      for (size_t istr=0; istr!=nstraws; ++istr) {

        const Straw& straw = _tt->getStraw(StrawIndex(istr));

        int cdev = straw.id().getDevice();
        int csec = straw.id().getSector();
        int clay = straw.id().getLayer();

        const Device& device = _tt->getDevice(cdev);
        const Sector& sector = device.getSector(csec);
        const Layer&  layer  = sector.getLayer(clay);

        size_t nStrawsPerSector = sector.nLayers()  * layer.nStraws();
        size_t nStrawsPerDevice = device.nSectors() * nStrawsPerSector;

        double cang = sector.boxRzAngle()/M_PI*180.;
        double dang = device.rotation()/M_PI*180.;
        double sroz = (sector.boxOffset() - device.origin()).z();

        size_t isecf = nStrawsPerSector*csec + nStrawsPerDevice*cdev;
 
        cout << __func__ << " Straw "
             << fixed << setw(6) << istr 
             << " secfloor " << setw(6) << isecf << " "
             << straw.id() 
             << " sector rotation: " << cang
             << " origin " << sector.boxOffset()
             << " device rotation: " << dang
             << " origin " << device.origin()
             << " sec rel origin z " << sroz
	     << " straw exists " << _tt->strawExists(StrawIndex(istr))
	     << " device exists " << device.exists();

        if (isec>csec && idev==cdev) cout << " <--S";
        if (iang>cang && idev==cdev) cout << " <--A";
        if (ilay>clay && isec==csec) cout << " <--L";
        if ((csec%2 == 0 && sroz>0.) || (csec%2 != 0 && sroz<0.)) cout << " <--Z";
        if (idev!=cdev) idev=cdev;
        if (isec!=csec) isec=csec;
        if (ilay!=clay) ilay=clay;

        cout << endl;

      }

    }

  }

  void TTrackerMaker::parseConfig( const SimpleConfig& config ){

    _verbosityLevel     = config.getInt("ttracker.verbosityLevel",0);


    _motherRIn        = config.getDouble("ttracker.mother.rIn"        )*CLHEP::mm;
    _motherROut       = config.getDouble("ttracker.mother.rOut"       )*CLHEP::mm;
    _motherHalfLength = config.getDouble("ttracker.mother.halfLength" )*CLHEP::mm;
    _motherZ0         = config.getDouble("ttracker.mother.z0"         )*CLHEP::mm;

    _numDevices         = config.getInt("ttracker.numDevices");
    _sectorsPerDevice   = config.getInt("ttracker.sectorsPerDevice");
    _layersPerSector    = config.getInt("ttracker.layersPerSector");
    _manifoldsPerEnd    = config.getInt("ttracker.manifoldsPerEnd");
    _strawsPerManifold  = config.getInt("ttracker.strawsPerManifold");
    _rotationPattern    = config.getInt("ttracker.rotationPattern");
    _spacingPattern     = config.getInt("ttracker.spacingPattern");

    _oddStationRotation   =  config.getDouble("ttracker.oddStationRotation")*CLHEP::degree;
    _zCenter              =  config.getDouble("ttracker.z0")*CLHEP::mm;
    _xCenter              = -config.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
    _envelopeInnerRadius  =  config.getDouble("ttracker.envelopeInnerRadius")*CLHEP::mm;
    _strawOuterRadius     =  config.getDouble("ttracker.strawOuterRadius")*CLHEP::mm;
    _strawWallThickness   =  config.getDouble("ttracker.strawWallThickness")*CLHEP::mm;
    _strawGap             =  config.getDouble("ttracker.strawGap")*CLHEP::mm;
    _deviceSpacing        =  config.getDouble("ttracker.deviceSpacing")*CLHEP::mm;
    _deviceHalfSeparation =  config.getDouble("ttracker.deviceHalfSeparation")*CLHEP::mm;

    _outerSupportRadius   =  config.getDouble("ttracker.outerSupportRadius")*CLHEP::mm;
    _innerSupportRadius   =  config.getDouble("ttracker.innerSupportRadius")*CLHEP::mm;
    _supportHalfThickness =  config.getDouble("ttracker.supportHalfThickness")*CLHEP::mm;
    _wireRadius           =  config.getDouble("ttracker.wireRadius")*CLHEP::mm;
    _manifoldYOffset      =  config.getDouble("ttracker.manifoldYOffset")*CLHEP::mm;
    _virtualDetectorHalfLength = config.getDouble("vd.halfLength")*CLHEP::mm;

    config.getVectorInt("ttracker.nonExistingDevices", _nonExistingDevices,  vector<int>() );

    _verbosityLevel > 0 && _nonExistingDevices.size()>0 &&
      cout << __func__ << " inactive devices : f/l   " 
	   << _nonExistingDevices.front() << " / "
	   << _nonExistingDevices.back()
	   << endl;

    // station
    // TODO Maybe -- These might eventually want to be config parameter driven
    _planesPerStation     = 2;
    _panelsPerFace        = 3;                    // hardwired 3 sectors/face
    _facesPerPlane        = _sectorsPerDevice / _panelsPerFace;
    _facesPerStation      = _planesPerStation*_facesPerPlane;
    _zLayersPerPanel      = _layersPerSector;
    _strawsPerZLayer      = _manifoldsPerEnd*_strawsPerManifold;

    config.getVectorDouble("ttracker.manifoldHalfLengths", _manifoldHalfLengths, 3);
    for ( size_t i=0; i<_manifoldHalfLengths.size(); ++i ){
      _manifoldHalfLengths.at(i) *= CLHEP::mm;
    }

    config.getVectorString("ttracker.strawMaterials", _strawMaterials, 3);

    _envelopeMaterial = config.getString("ttracker.mat.vacuum");
    _supportMaterial = config.getString("ttracker.mat.support");

    _passivationMargin        =  config.getDouble("ttracker.passivationMargin")*CLHEP::mm;

    // For the detailv0 support model there are lots of extra parameters.
    _supportModel = SupportModel( config.getString("ttrackerSupport.model","simple"));
    if ( _supportModel == SupportModel::detailedv0 ) {
      _endRingOuterRadius      = config.getDouble( "ttrackerSupport.endRing.outerRadius" );
      _endRingInnerRadius      = config.getDouble( "ttrackerSupport.endRing.innerRadius" );
      _endRingHalfLength       = config.getDouble( "ttrackerSupport.endRing.halfLength"  );
      _endRingZOffset          = config.getDouble( "ttrackerSupport.endRing.zOffset"     );
      _endRingMaterial         = config.getString( "ttrackerSupport.endRing.material"    );

      config.getVectorInt( "ttrackerSupport.midRing.slot", _midRingSlot );
      _midRingHalfLength       = config.getDouble(    "ttrackerSupport.midRing.halfLength" );


      _nStaves                 = config.getInt   ( "ttrackerSupport.stave.nstaves"       );
      _stavePhi0               = config.getDouble( "ttrackerSupport.stave.phi0"          );
      _stavePhiWidth           = config.getDouble( "ttrackerSupport.stave.phiWidth"      );
      _staveMaterial           = config.getString( "ttrackerSupport.stave.material"      );

      _innerRingInnerRadius    = config.getDouble( "ttrackerSupport.innerRing.innerRadius" );
      _innerRingOuterRadius    = config.getDouble( "ttrackerSupport.innerRing.outerRadius" );
      _innerRingHalfLength     = config.getDouble( "ttrackerSupport.innerRing.halfLength"  );
      _innerRingMaterial       = config.getString( "ttrackerSupport.innerRing.material"    );

      _centerPlateHalfLength   = config.getDouble( "ttrackerSupport.centerPlate.halfLength" );
      _centerPlateMaterial     = config.getString( "ttrackerSupport.centerPlate.material"   );

      _outerRingInnerRadius    = config.getDouble( "ttrackerSupport.outerRing.innerRadius" );
      _outerRingOuterRadius    = config.getDouble( "ttrackerSupport.outerRing.outerRadius" );
      _outerRingMaterial       = config.getString( "ttrackerSupport.outerRing.material"    );

      _coverHalfLength          = config.getDouble( "ttrackerSupport.cover.halfLength"           );
      _coverMaterial            = config.getString( "ttrackerSupport.cover.material"             );
      _electronicsG10HalfLength = config.getDouble( "ttrackerSupport.electronics.g10.halfLength" );
      _electronicsG10Material   = config.getString( "ttrackerSupport.electronics.g10.material"   );
      _electronicsCuHhalfLength = config.getDouble( "ttrackerSupport.electronics.cu.halfLength"  );
      _electronicsCuMaterial    = config.getString( "ttrackerSupport.electronics.cu.material"    );
      _channelZOffset           = config.getDouble( "ttrackerSupport.channel.zOffset"            );
      _channelDepth             = config.getDouble( "ttrackerSupport.channel.depth"              );
      _channelMaterial          = config.getString( "ttrackerSupport.channel.material"           );
      _electronicsSpaceMaterial = config.getString( "ttrackerSupport.electronicsSpace.material"  );

      _wallOuterMetalThickness  = config.getDouble("ttracker.straw.wallOuterMetal.thickness")*CLHEP::mm;
      _wallInnerMetal1Thickness = config.getDouble("ttracker.straw.wallInnerMetal1.thickness")*CLHEP::mm;
      _wallInnerMetal2Thickness = config.getDouble("ttracker.straw.wallInnerMetal2.thickness")*CLHEP::mm;
      _wirePlateThickness       = config.getDouble("ttracker.straw.wirePlate.thickness")*CLHEP::mm;
      _wallOuterMetalMaterial   = config.getString("ttracker.straw.wallOuterMetal.material");
      _wallInnerMetal1Material  = config.getString("ttracker.straw.wallInnerMetal1.material");
      _wallInnerMetal2Material  = config.getString("ttracker.straw.wallInnerMetal2.material");
      _wirePlateMaterial        = config.getString("ttracker.straw.wirePlate.material");

    }

    //string ttracker.mat.manifold  = "G4_Al";  // Placeholder.

    // Also define some parameters that may become variable some day.
    _sectorBaseRotations.clear();
    _sectorZSide.clear();

    if (_sectorsPerDevice == 6 ){
  // the Z pattern here is forced by the 'alternating sector id' convention
  // that permeates the rest of the code.
	_sectorZSide.push_back(-1.);
	_sectorZSide.push_back(+1.);
	_sectorZSide.push_back(-1.);
	_sectorZSide.push_back(+1.);
	_sectorZSide.push_back(-1.);
	_sectorZSide.push_back(+1.);
      if(_rotationPattern == 1){
// cdr geometry, taken from DOC 888, also alternatives 1 and 3 from doc 2799
	// faces overlap by 60 degrees
// Implicitly define the rotations for the even and odd (sequentially) sectors.
	_sectorBaseRotations.push_back(	 45.*CLHEP::degree);
	_sectorBaseRotations.push_back( 105.*CLHEP::degree);
	_sectorBaseRotations.push_back( 165.*CLHEP::degree);
	_sectorBaseRotations.push_back( 225.*CLHEP::degree);
	_sectorBaseRotations.push_back( 285.*CLHEP::degree);
	_sectorBaseRotations.push_back( 345.*CLHEP::degree);
	_sectorBaseRotations.push_back(  75.*CLHEP::degree);
	_sectorBaseRotations.push_back(  15.*CLHEP::degree);
	_sectorBaseRotations.push_back( 195.*CLHEP::degree);
	_sectorBaseRotations.push_back( 135.*CLHEP::degree);
	_sectorBaseRotations.push_back( 315.*CLHEP::degree);
	_sectorBaseRotations.push_back( 255.*CLHEP::degree);
		
      } else if(_rotationPattern==2){
      	// alternative 2 from DOC 2799
	// faces overlap by 60 degrees
	_sectorBaseRotations.push_back(  45.*CLHEP::degree);
	_sectorBaseRotations.push_back(  75.*CLHEP::degree);
	_sectorBaseRotations.push_back( 165.*CLHEP::degree);
	_sectorBaseRotations.push_back( 195.*CLHEP::degree);
	_sectorBaseRotations.push_back( 285.*CLHEP::degree);
	_sectorBaseRotations.push_back( 315.*CLHEP::degree);
	_sectorBaseRotations.push_back( 105.*CLHEP::degree);
	_sectorBaseRotations.push_back(  15.*CLHEP::degree);
	_sectorBaseRotations.push_back( 225.*CLHEP::degree);
	_sectorBaseRotations.push_back( 135.*CLHEP::degree);
	_sectorBaseRotations.push_back( 345.*CLHEP::degree);
     	_sectorBaseRotations.push_back( 255.*CLHEP::degree);

      } else if(_rotationPattern==3){
      // faces overlap by 60 degrees, second device 'flipped'
	_sectorBaseRotations.push_back(   0.*CLHEP::degree);
	_sectorBaseRotations.push_back(  90.*CLHEP::degree);
	_sectorBaseRotations.push_back(  120.*CLHEP::degree);
	_sectorBaseRotations.push_back(  210.*CLHEP::degree);
	_sectorBaseRotations.push_back(  240.*CLHEP::degree);
	_sectorBaseRotations.push_back(  330.*CLHEP::degree);
	_sectorBaseRotations.push_back(  30.*CLHEP::degree);
	_sectorBaseRotations.push_back(  60.*CLHEP::degree);
	_sectorBaseRotations.push_back(  150.*CLHEP::degree);
	_sectorBaseRotations.push_back(  180.*CLHEP::degree);
	_sectorBaseRotations.push_back(  270.*CLHEP::degree);
	_sectorBaseRotations.push_back(  300.*CLHEP::degree);
      } else {
	throw cet::exception("GEOM")
	  << "Unrecognized rotation pattern in TTrackerMaker. \n";
      }
    } else {
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
    _tt = unique_ptr<TTracker>(new TTracker());

    makeMother();

    computeLayerSpacingAndShift();
    computeManifoldEdgeExcessSpace();

    _tt->_supportModel = _supportModel;

    // Fill information about the new style, fully detailed support structure.
    if ( _supportModel == SupportModel::detailedv0 ) {
      makeSupportStructure();
    }

    // Fill the information about the, old style minimal supports.
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

    makeDetails();

    // Reserve space for straws so that pointers are valid.
    _nStrawsToReserve = _numDevices * _sectorsPerDevice * _layersPerSector *
      _manifoldsPerEnd * _strawsPerManifold;
    //_tt->_allStraws.reserve(_nStrawsToReserve); // see makeLayer

    // we need to be able to assign true/false for a given straw index
    _tt->_strawExists.resize(_nStrawsToReserve);

    _tt->_devices.reserve(_numDevices);
    // Construct the devices and their internals.
    for ( int idev=0; idev<_numDevices; ++idev ){
      makeDevice( DeviceId(idev) );
    }

    // Fill all of the non-persistent information.
    _tt->fillPointers();

    identifyNeighbourStraws();
    identifyDirectionalNeighbourStraws();

   // Stations
   // Construct the stations and their internals based on devices internals
   if ( _numDevices%2 != 0 ) {
         throw cet::exception("GEOM")  << "_numDevices = " << _numDevices
         << ": Current TTracker geometry assumes even number of devices  \n";
   }
   _numStations = _numDevices/2;
    _tt->_stations.reserve(_numStations);
    // Construct the devices and their internals.
    for ( int istation=0; istation<_numStations; ++istation ){
      makeStation( StationId(istation) );
    }

    // Order is important here.
    computeDeviceEnvelope();
    computeTrackerEnvelope();
    recomputeHalfLengths();
    makeStrawTubs();

    // This uses information from the devices
    makeThinSupportRings();


    finalCheck();

    if ( _verbosityLevel > 0 ) {
      cout << "TTracker Support Structure: \n" << _tt->_supportStructure << endl;
    }

    // Test the forAll methods.
    //_tt->forAllLayers( lptest);
    //_tt->forAllDevices( devtest);
    //_tt->forAllLayers( positionTest);

  } //end TTrackerMaker::buildIt.

  void TTrackerMaker::makeMother(){

    _tt->_mother = PlacedTubs ( "TrackerMother",
                                TubsParams( _motherRIn, _motherROut, _motherHalfLength),
                                CLHEP::Hep3Vector( _xCenter, 0., _motherZ0),
                                _envelopeMaterial );

  } //end TTrackerMaker::makeMother

  // In the present code the straw positions are computed using the manifold information.
  // The channel position is computed using the SupportStructure information.  These two
  // must be self consistent.  Relatively soon the manifold information will be removed
  // and they will be replaced with the support structure information.
  // For now, check for self-consistency.
  void TTrackerMaker::finalCheck( ){

    // Only do this test for the new model(s).
    if ( _supportModel == SupportModel::simple ) return;

    double ztest = 0.5*(_tt->getStraw ( StrawId(0,1,0,0) ).getMidPoint().z()+_tt->getStraw ( StrawId(0,1,1,0) ).getMidPoint().z())
      -  _tt->getDevice( DeviceId(0) ).origin().z();

    double tolerance = 1.e-6;
    if ( std::abs(ztest-_channelZOffset) > tolerance ){
      throw cet::exception("GEOM")  << "Inconsistent channel center and wire z location. \n"
                                    << "channel Z offset: " <<  _channelZOffset << "\n"
                                    << "device z center:  " << _tt->getDevice( DeviceId(0) ).origin().z() << "\n "
                                    << "Straw z layer 0:  " << _tt->getStraw ( StrawId(0,1,0,0) ).getMidPoint().z() << "\n"
                                    << "Straw z layer 1:  " << _tt->getStraw ( StrawId(0,1,1,0) ).getMidPoint().z() << "\n"
                                    << "z test:           " << ztest << " delta " << std::abs(ztest-_channelZOffset) << "\n";
    }

  }


  void TTrackerMaker::makeDevice( DeviceId devId ){

//std::cout << "->->-> makeDevice\n";
    int idev = devId;

    double devDeltaZ = chooseDeviceSpacing(idev);
    CLHEP::Hep3Vector origin( 0., 0., _z0+devDeltaZ);
    
    // device rotation is no longer used.
    double phi = 0.0;
    _tt->_devices.push_back(Device(devId, origin, phi));
    Device& dev = _tt->_devices.back();
    dev._exists = ( find ( _nonExistingDevices.begin(), _nonExistingDevices.end(), idev) == 
		      _nonExistingDevices.end() );

    dev._sectors.reserve(_sectorsPerDevice);

    for ( int isec=0; isec<_sectorsPerDevice; ++isec ){
      makeSector ( SectorId(devId,isec), dev );
    }

//std::cout << "<-<-<- makeDevice\n";
  }

  void TTrackerMaker::makeSector( const SectorId& secId, Device& dev ){
//std::cout << "->->-> makeSector\n";

    dev._sectors.push_back( Sector(secId) );
    Sector& sector = dev._sectors.back();
    sector._layers.reserve(_layersPerSector);

    // check if the opposite sectors do not overlap
    static double const tolerance = 1.e-6; // this should be in a config file

    if ((2.*_manifoldHalfLengths.at(2)+_supportHalfThickness)>_deviceHalfSeparation + tolerance) {
      cout << "(2.*_manifoldHalfLengths.at(2)+_supportHalfThickness), _deviceHalfSeparation " <<
        (2.*_manifoldHalfLengths.at(2)+_supportHalfThickness) << ", " <<_deviceHalfSeparation << endl;
      throw cet::exception("GEOM")  << "Devices are too close \n";
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

//std::cout << "<-<-<- makeSector\n";
  }  // makeSector

  void TTrackerMaker::makeLayer ( const LayerId& layId, Sector& sector ){
//std::cout << "->->-> makeLayer\n";

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
    CLHEP::HepRotationZ RZ(sectorRotation(isec,layId.getDeviceId()));

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

	_tt->_strawExists[index.asInt()] = device.exists();

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

//std::cout << "<-<-<- makeLayer\n";
  }

  void TTrackerMaker::makeManifolds( const SectorId& secId){
//std::cout << "->->-> makeManifolds\n";

    if ( _sectorsPerDevice != 4 && _sectorsPerDevice != 6 ) {
      throw cet::exception("GEOM")
        << "This code only knows how to do 4 or 6 sectors per device.\n";
    }

//    double phi = _tt->getDevice(secId.getDevice()).rotation();
//    CLHEP::HepRotationZ RZ(phi);


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

//std::cout << "<-<-<- makeManifolds\n";
  }

// ======= Station view makers ============

  void TTrackerMaker::makeStation( StationId stationId ){
//std::cout << "->->-> makeStation\n";

    int ist = stationId;
    int idev1 = 2*ist;
    int idev2 = idev1 + 1;
    double stationZ = 0.5 *
        ( _tt->_devices.at(idev1).origin().z() +
          _tt->_devices.at(idev2).origin().z() );
    _tt->_stations.push_back(Station(stationId, stationZ));
    Station & st = _tt->_stations.back();
    st._planes.reserve (_planesPerStation);
    st._faces.reserve (_facesPerStation);

    for ( int iplane = 0; iplane < _planesPerStation; ++iplane ) {
      makePlane ( PlaneId ( stationId, iplane ), st );
    }

//std::cout << "<-<-<- makeStation\n";
  }

  // station view
  void TTrackerMaker::makePlane( const PlaneId& planeId, Station & station ){
//std::cout << "->->-> makePlane\n";

    int idev = 2*planeId.getStation() + planeId.getPlane();
    const Device & device = _tt->_devices.at(idev);
    double planeZ = device.origin().z();
    station._planes.push_back(Plane(planeId, planeZ));
    Plane & plane = station._planes.back();
    plane._faces.reserve(_facesPerPlane);
    for ( int iface = 0; iface < _facesPerPlane; ++iface ) {
      int faceNum = 2*planeId.getPlane() + iface;
      makeFace ( FaceId ( station.id(), faceNum ), plane, device );
      station._faces.push_back(&(plane._faces.back()));  // quelle haque
    }

//std::cout << "<-<-<- makePlane\n";
  }  // makePlane

// station view
  void TTrackerMaker::makeFace( const FaceId & faceId,
                                      Plane  & plane,
                                const Device & device )
  {
//std::cout << "->->-> makeFace " << faceId << "\n";
    bool faceZisKnown = false;
    double faceZ = 0.0;
    double faceZtolerance = .0001;
    int faceParity = faceId.getFace()%2;

    for (int isector = faceParity; isector < _sectorsPerDevice; isector += 2) {
      const Sector & sector = device.getSector(isector);
      if (faceZisKnown) {
        if (abs ( sector.straw0MidPoint().z() - faceZ ) > faceZtolerance ) {
          throw cet::exception ( "GEOM" )
          << "Inconsistent sector positions within a face for faceId = "
          << faceId << "\nExpected Z position = " << faceZ
          << "  Sector base position Z = " << sector.straw0MidPoint().z()
          << " in sector " << sector.id() << "\n";
        }
      } else {
        faceZ = sector.straw0MidPoint().z();
        // TODO -- If this is absolutely correct, I will be surprised.
        //         We need to look at where this places the faces, and
        //         that will clue us in to how to do it right.
        faceZisKnown = true;
      }
    }

    plane._faces.push_back(Face(faceId, faceZ));
    Face & f =  plane._faces.back();
    f._panels.reserve  (_panelsPerFace);
    for ( int ipanel = 0; ipanel < _panelsPerFace; ++ipanel ) {
      makePanel ( PanelId ( faceId, ipanel ), f,
                  device.getSector(2*ipanel + faceParity) );
    }

//std::cout << "<-<-<- makeFace\n";
  }  // makeFace

  // station view
  void TTrackerMaker::makePanel( const PanelId & panelId,
                                       Face  & face,
                                 const Sector & sector )
  {
//std::cout << "->->-> makePanel " << panelId << "\n";
    double panelZ = face.midZ();
    face._panels.push_back(Panel(panelId, panelZ));
    Panel & pnl =  face._panels.back();
    pnl._zlayers.reserve  (_zLayersPerPanel);
    for ( int izlayer = 0; izlayer < _zLayersPerPanel; ++izlayer ) {
      int faceParity = face.id().getFace()%2;
      int ilayer = 1 - faceParity + (2*faceParity -1)*izlayer;
      makeZLayer ( ZLayerId ( panelId, izlayer ), pnl,
                                                  sector.getLayer(ilayer) );
    }
    // Determine phi (and view) by looking at just one sample straw;
    // check that the remaining straws are consisdtent direction-wise
    Straw const * sampleStraw = pnl.getZLayer(0).getStrawptr(0);
    CLHEP::Hep3Vector dir = sampleStraw->getDirection();
    const double directionTolerance = .01e-8;
    typedef  std::vector<ZLayer>::const_iterator ZLayersIt;
    ZLayersIt izlend =  pnl.getZLayers().end();
    for ( ZLayersIt izl = pnl.getZLayers().begin(); izl != izlend; ++izl )
    {
      typedef std::vector<Straw const *>::const_iterator StrawptrsIt;
      StrawptrsIt isend = izl->getStraws().end();
      for ( StrawptrsIt is = izl->getStraws().begin(); is != isend; ++is )
      {
        CLHEP::Hep3Vector dir2 = (*is)->getDirection();
        double directionDiscrepancy =
          (dir.x() - dir2.x())*(dir.x() - dir2.x()) +
          (dir.y() - dir2.y())*(dir.y() - dir2.y()) +
          dir2.z()*dir2.z();
        // Note that a significantly non-zero z component of any direction
        // is also clasified as problematic
        if ( directionDiscrepancy > directionTolerance ) {
          throw cet::exception("GEOM")
            << "makePanel: Straw Direction Discrepancy \n"
            << "Straw " << sampleStraw->id() << " direction " << dir  << "\n"
            << "Straw " << (*is)->id()       << " direction " << dir2 << "\n";
        }
      }
    }
    double phi = dir.phi() - M_PI/2;
    if (phi < 0) phi += 2*M_PI;  // phi is now in range [0,2 pi)
    if ( (phi < 0) || (phi >= 2*M_PI) ) {
//       throw cet::exception("GEOM")
//	std::cout
//         << "makePanel: An ill-understood phi - "
//         << dir.phi() << phi-2*M_PI << std::endl;
    }
    pnl._phi = phi;
    double phi0 = 0.0;  // TODO -- cope with non-zero station rotation offsets
    double clockAngle = phi - phi0 + M_PI/36;
    if (clockAngle < 0)       clockAngle += 2*M_PI;
    if (clockAngle > 2*M_PI ) clockAngle -= 2*M_PI;
    int hour = static_cast<int> (std::floor(clockAngle/(M_PI/6)));
    if ( (hour < 0 || hour >= 12) ) {
       throw cet::exception("GEOM")
         << "makePanel: An ill-understood conversion of phi to hour - "
         << dir.phi() << " --> " << hour;
    }
    pnl._view = hour;
//std::cout << "<-<-<- makePanel\n";
  }  // makePanel

  // station view
  void TTrackerMaker::makeZLayer( const ZLayerId & zlayerId,
                                        Panel  & panel,
                                  const Layer & layer )
  {
    double layerZ = layer.straw0MidPoint().z();
    panel._zlayers.push_back(ZLayer(zlayerId, layerZ));
    ZLayer & zlay =  panel._zlayers.back();
    zlay._straws.reserve  (_strawsPerZLayer);
    zlay._straws = layer.getStraws();
  }  // makeZLayer

// ======= end of Station view makers ================

  // Assumes all devices and all sectors are the same.
  // The straw length depends only on the manifold number.
  // See Mu2e-doc-??? for the algorithm.
  void TTrackerMaker::computeStrawHalfLengths(){

    _strawHalfLengths.clear();
    _strawActiveHalfLengths.clear();

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

      double yA = sqrt( diff_of_squares(_innerSupportRadius, xA) );
      double yB = yA + _manifoldYOffset;
      _strawHalfLengths.push_back(yB);

      // This variable is only used if SupportModel==simple.
      double activeHLen=sqrt( diff_of_squares(_innerSupportRadius,xA+2.5))-_passivationMargin;
      activeHLen = std::max( activeHLen, 1.0);
      _strawActiveHalfLengths.push_back(activeHLen);

    }

  } // end TTrackerMaker::computeStrawHalfLengths

  void TTrackerMaker::makeDetails(){

    computeConstantSectorBoxParams();

    if ( _supportModel == SupportModel::simple ){
      _tt->_strawDetails.reserve(_manifoldsPerEnd);
    } else{
      // This will be extended in recomputeStrawLengths - reserve enough space for the extention.
      _tt->_strawDetails.reserve(_manifoldsPerEnd*_layersPerSector);
    }

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      _tt->_strawDetails.push_back
        ( StrawDetail
          ( i,
            _strawMaterials,
            _strawOuterRadius,
            _strawWallThickness,
            _strawHalfLengths.at(i),
            _strawActiveHalfLengths.at(i),
            _wireRadius
            )
          );
    }

  } // end TTrackerMaker::makeDetails

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
    sector._boxRzAngle = sectorRotation(isec,dev.id() );

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

  } // end TTrackerMaker::computeSectorBoxParams

  void TTrackerMaker::computeLayerSpacingAndShift(){

    _layerHalfSpacing = (_layersPerSector<=1) ? 0.0 :
      sqrt(3.0*(square(_strawOuterRadius)+_strawOuterRadius*_strawGap+0.25*square(_strawGap)))*0.5;
    _layerHalfShift   = (_layersPerSector<=1) ? 0.0 : _strawGap*0.25 + _strawOuterRadius*0.5;

  }

  void TTrackerMaker::computeManifoldEdgeExcessSpace(){

    // Computes space between first/last straw and edge of manifold

    _manifoldXEdgeExcessSpace = _manifoldHalfLengths.at(0) -
      _strawOuterRadius*_strawsPerManifold -
      _strawGap*(_strawsPerManifold-1)*0.5;

    _manifoldZEdgeExcessSpace = _manifoldHalfLengths.at(2) - _strawOuterRadius -
      (_layersPerSector-1)*_layerHalfSpacing;

    //cout << "Debugging,  _manifoldXEdgeExcessSpace, _manifoldZEdgeExcessSpace: " <<
    //     _manifoldXEdgeExcessSpace << ", " << _manifoldZEdgeExcessSpace << endl;

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

    // the order is forced by the nestTrp/G4Trd and Sector data

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
      sqrt( sum_of_squares(_envelopeInnerRadius + 2.0*_sectorBoxHalfLengths.at(1),
                           _sectorBoxHalfLengths.at(3)
                           )
            );
    double outerSupportRadiusRequiredl =
      max( outerSupportRadiusRequireds,
           sqrt( sum_of_squares(_envelopeInnerRadius-pad,
                                _sectorBoxHalfLengths.at(4))
                 )
           );

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

  } // end TTrackerMaker::computeConstantSectorBoxParams


  // Compute the spacing for the given device.
  double TTrackerMaker::chooseDeviceSpacing( int idev ) const {


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
    return 0.0;
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

        // TODO -- CORRECT A LOGIC ERROR (or check this reasoning is amiss):
        //
        // The block below adds straw n +/- 1 in the opposite layer.  But
        // when n is 0 and layer is 1 (or when n is nStrawLayer-1 and layer
        // is zero) it is supposed to add straw 1 of layer 0 (or straw
        // nStrawLayer-2 of layer 1).  The computation is in fact correct,
        // but the check will cause the insertioin to be skipped.  Thus
        // one neighbor of each of two straws in the panel will be omitted.

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
  } // identifyNeighborStraws

  // Identify the neighbour straws in inner/outer same-layer or zig-zag order
  void TTrackerMaker::identifyDirectionalNeighbourStraws() {

    // TODO:  The following algorithm relies on a few more geometry assumptions
    //        than would strictly be necessary.  For example, it relies on
    //        the layer labelling such that leayer 0 lies half a straw to the
    //        inner side of layer 1, in each panel.  Some of these assumptions
    //        can be lifted, and others at least checked, so that if the
    //        geometry changes, the code will still produce the right results,
    //        or will at least throw to indicate a serious problem.

    deque<Straw>& allStraws = _tt->_allStraws;

    for (deque<Straw>::iterator straw = allStraws.begin();
         straw != allStraws.end(); ++straw) {

      // throw exception if more than 2 layers in the sector of this straw
      if (_tt->getSector(straw->id().getSectorId()).nLayers() != 2 ) {
        throw cet::exception("GEOM")
          << "The code expects exactly 2 layers per sector. \n";
      }

      LayerId layerId = straw->id().getLayerId();
      int layerNumber = layerId.getLayer();
      bool layerStaggeredToInside = (layerNumber == 0);
      // In all cases, layer 0 is staggered to the inside,
      // layer 1 is staggered to the outside.  We will now check this:
      // TODO -- Do the check
      LayerId otherLayerId ( layerId.getSectorId(), 1-layerNumber );

      int nStrawLayer = _tt->getLayer(layerId)._nStraws;

      int strawNumberWithinLayer = straw->id().getStraw();
      int incrementedStrawNumber =
        ( strawNumberWithinLayer + 1 < nStrawLayer )
        ? strawNumberWithinLayer + 1
        : StrawIndex::NO_STRAW;
      int decrementedStrawNumber =
        ( strawNumberWithinLayer - 1 >=  0)
        ? strawNumberWithinLayer - 1
        : StrawIndex::NO_STRAW;

      straw->_nextOuterL = ttStrawIndex (layerId, incrementedStrawNumber);
      straw->_nextInnerL = ttStrawIndex (layerId, decrementedStrawNumber);
      if (layerStaggeredToInside) {
        straw->_nextOuterP = ttStrawIndex(otherLayerId, strawNumberWithinLayer);
        straw->_nextInnerP = ttStrawIndex(otherLayerId, decrementedStrawNumber);
      } else {
        straw->_nextOuterP = ttStrawIndex(otherLayerId, incrementedStrawNumber);
        straw->_nextInnerP = ttStrawIndex(otherLayerId, strawNumberWithinLayer);
      }

#ifdef CHECK_STRAW_NAVIGATION_ASSIGNMENTS_BY_PRINTING
      cout << "Straw " << straw->id() << ":\n"
           << "_nextOuterL: "   << straw->_nextOuterL.asInt()
           << "  _nextInnerL: " << straw->_nextInnerL.asInt()
           << "\n_nextOuterP: " << straw->_nextOuterP.asInt()
           << "  _nextInnerP: " << straw->_nextInnerP.asInt() << "\n";
#endif

      // TODO -- Insert logic to check that the radius of the purported
      // next straw differs by the right amount and sign, in each of these
      // four cases.

    } // end of loop over all straws

  } // identifyDirectionalNeighbourStraws

  StrawIndex TTrackerMaker::ttStrawIndex(LayerId const & layerId, int snum)
  {
    if ( snum == StrawIndex::NO_STRAW ) {
      return  StrawIndex(StrawIndex::NO_STRAW);
    }
    StrawId sid ( layerId, snum );
    return _tt->getStraw(sid).index();
  } //; ttStrawIndex

  void TTrackerMaker::makeSupportStructure(){

    SupportStructure& sup  = _tt->_supportStructure;

    // Positions for the next few objects in Mu2e coordinates.
    TubsParams endRingTubs( _endRingInnerRadius, _endRingOuterRadius, _endRingHalfLength);

    TubsParams midRingTubs ( _endRingInnerRadius, _endRingOuterRadius, _midRingHalfLength);
    sup._stiffRings.push_back(PlacedTubs ( "TTrackerEndRingUpstream",   endRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter-_endRingZOffset), _endRingMaterial ));

    {
      if ( _numDevices%2 !=0 ){
        throw cet::exception("GEOM")
          << "TTrackerMaker::makeSupportStructure expected and even number of devices. Saw " << _numDevices << " devices.\n";
      }

      // From upsgream end of most upstream station to the downstream end of the most downstream station.
      // Including all materials.
      double overallLength = (_numDevices/2-1)*_deviceSpacing + 2.*_deviceHalfSeparation + 2.* _innerRingHalfLength;

      // Staves touch the big ring at the upstream end; they end flush with the downstream edge of the most downstream station.
      double z1 = -(_endRingZOffset-_endRingHalfLength);
      double z2 = overallLength/2.;
      double zoff = (z1+z2)/2.;
      double zHalf = ( z2-z1)/2.;

      double dphi       = 360*CLHEP::degree/_nStaves;
      double phi0       = _stavePhi0*CLHEP::degree;
      double pad        = 2.*_virtualDetectorHalfLength;
      TubsParams staveTubsParams( _outerRingOuterRadius+pad, _endRingOuterRadius, zHalf, 0., _stavePhiWidth*CLHEP::degree);
      for ( int i=0; i<_nStaves; ++i ){
        double phi = phi0 + i*dphi;
        CLHEP::HepRotationZ rot(phi);
        ostringstream os;
        os << "TTrackerStave_" << i;
        sup._staveBody.push_back( PlacedTubs ( os.str(), staveTubsParams, CLHEP::Hep3Vector(_xCenter,0., _zCenter+zoff), rot, _staveMaterial) );
      }
    }

    // Positions from here onward are in the coordinates of the device envelope.
    {
      TubsParams centerPlateTubs( _innerRingOuterRadius, _outerRingOuterRadius, _centerPlateHalfLength);
      sup._centerPlate = PlacedTubs( "TTrackerSupportCenterPlate", centerPlateTubs, CLHEP::Hep3Vector(0.,0.,0.), _centerPlateMaterial);
    }

    {
      TubsParams innerRingTubs( _innerRingInnerRadius, _innerRingOuterRadius, _innerRingHalfLength);
      sup._innerRing = PlacedTubs( "TTrackerSupportInnerRing", innerRingTubs, CLHEP::Hep3Vector(0.,0.,0), _innerRingMaterial );
    }

    {
      double outerRadius        = _innerRingInnerRadius + _channelDepth;
      double channelHalfLength  = (_layerHalfSpacing-_strawOuterRadius) + 2.*_strawOuterRadius;
      TubsParams innerChannelTubs( _innerRingInnerRadius, outerRadius, channelHalfLength );
      sup._innerChannelUpstream   = PlacedTubs( "TTrackerSupportInnerChannelUpstream",   innerChannelTubs, CLHEP::Hep3Vector(0.,0.,-_channelZOffset), _envelopeMaterial );
      sup._innerChannelDownstream = PlacedTubs( "TTrackerSupportInnerChannelDownstream", innerChannelTubs, CLHEP::Hep3Vector(0.,0., _channelZOffset), _envelopeMaterial );
    }

    {
      double halfLength = (_innerRingHalfLength - _centerPlateHalfLength)/2.;
      double dz         = _centerPlateHalfLength + halfLength;
      TubsParams outerRingTubs( _outerRingInnerRadius, _outerRingOuterRadius, halfLength);
      sup._outerRingUpstream   = PlacedTubs ( "TTrackerSupportOuterRingUpstream",   outerRingTubs, CLHEP::Hep3Vector(0.,0.,-dz), _outerRingMaterial );
      sup._outerRingDownstream = PlacedTubs ( "TTrackerSupportOuterRingDownstream", outerRingTubs, CLHEP::Hep3Vector(0.,0., dz), _outerRingMaterial );
    }

    {
      double dz = _innerRingHalfLength-_coverHalfLength;
      TubsParams coverTubs( _innerRingOuterRadius, _outerRingInnerRadius, _coverHalfLength);
      sup._coverUpstream   = PlacedTubs( "TTrackerSupportCoverUpstream",   coverTubs, CLHEP::Hep3Vector(0.,0.,-dz), _coverMaterial );
      sup._coverDownstream = PlacedTubs( "TTrackerSupportCoverDownstream", coverTubs, CLHEP::Hep3Vector(0.,0., dz), _coverMaterial );
    }

    {
      double halfLength = (_innerRingHalfLength - _centerPlateHalfLength - 2.*_coverHalfLength)/2.;
      double dz         = _centerPlateHalfLength + halfLength;
      TubsParams gasTubs( _innerRingOuterRadius, _outerRingInnerRadius, halfLength );
      sup._gasUpstream   = PlacedTubs ( "TTrackerSupportGasUpstream",  gasTubs, CLHEP::Hep3Vector(0.,0.,-dz), _electronicsSpaceMaterial );
      sup._gasDownstream = PlacedTubs ( "TTrackerSupportGasDownstream", gasTubs, CLHEP::Hep3Vector(0.,0., dz), _electronicsSpaceMaterial );
    }

    // Positions for the next two are in the coordinates of the electronics space.
    {
      TubsParams g10Tubs( _innerRingOuterRadius, _outerRingInnerRadius, _electronicsG10HalfLength );
      sup._g10Upstream   = PlacedTubs ( "TTrackerSupportElecG10Upstream",   g10Tubs, CLHEP::Hep3Vector(0.,0.,-_electronicsG10HalfLength), _electronicsG10Material );
      sup._g10Downstream = PlacedTubs ( "TTrackerSupportElecG10Downstream", g10Tubs, CLHEP::Hep3Vector(0.,0.,-_electronicsG10HalfLength), _electronicsG10Material );
    }

    {
      TubsParams cuTubs( _innerRingOuterRadius, _outerRingInnerRadius, _electronicsCuHhalfLength);
      sup._cuUpstream   = PlacedTubs ( "TTrackerSupportElecCuUpstream",   cuTubs, CLHEP::Hep3Vector(0.,0.,_electronicsCuHhalfLength), _electronicsCuMaterial);
      sup._cuDownstream = PlacedTubs ( "TTrackerSupportElecCuDownstream", cuTubs, CLHEP::Hep3Vector(0.,0.,_electronicsCuHhalfLength), _electronicsCuMaterial);
    }

  }

  // This needs to know the z positions of the
  void TTrackerMaker::makeThinSupportRings(){
    SupportStructure& sup  = _tt->_supportStructure;

    TubsParams thinRingTubs ( _endRingInnerRadius, _outerRingOuterRadius, _midRingHalfLength);

    for ( size_t i=0; i< _midRingSlot.size(); ++i){
      ostringstream name;
      int station = _midRingSlot.at(i);
      int idev1 = station*2+1;
      int idev2 = idev1+1;
      if ( idev2 >= _tt->nDevices() ){
        throw cet::exception("GEOM")
          << "Requested a thin support after station: "
          << station
          << " This is between devices: "
          << idev1 << " and "
          << idev2 << "\n"
          << "But there are only "
          << _tt->nDevices()
          << " devices\n";
      }
      name << "ThinSupportRing_" << i;

      // Center the support in the gap between two stations.
      double z = 0.5*( _tt->getDevice(idev1).origin().z() +  _tt->getDevice(idev2).origin().z());
      sup._stiffRings.push_back(PlacedTubs ( name.str(),  thinRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter+z), _endRingMaterial ));
    }

  }

  // Create all of the Tubs objects needed for G4 to describe a straw.
  void TTrackerMaker::makeStrawTubs(){
    std::vector<StrawDetail>& details = _tt->_strawDetails;

    for ( std::vector<StrawDetail>::iterator i=details.begin();
          i != details.end(); ++i ){

      StrawDetail& detail(*i);

      detail._outerMetalThickness  = _wallOuterMetalThickness;
      detail._innerMetal1Thickness = _wallInnerMetal1Thickness;
      detail._innerMetal2Thickness = _wallInnerMetal2Thickness;
      detail._wirePlateThickness   = _wirePlateThickness;
      detail._outerMetalMaterial   = _wallOuterMetalMaterial;
      detail._innerMetal1Material  = _wallInnerMetal1Material ;
      detail._innerMetal2Material  = _wallInnerMetal2Material;
      detail._wirePlateMaterial    = _wirePlateMaterial;

    }
  }

  // Envelope that holds one device ("TTrackerDeviceEnvelope")
  void TTrackerMaker::computeDeviceEnvelope(){

    if ( _supportModel == SupportModel::simple ){
      double halfThick = _tt->_supportParams.halfThickness() + 2.*_tt->_manifoldHalfLengths[2];
      _tt->_deviceEnvelopeParams = TubsParams( _tt->_envelopeInnerRadius,
                                               _tt->_supportParams.outerRadius(),
                                               halfThick);
    } else if ( _supportModel == SupportModel::detailedv0 ){
      _tt->_deviceEnvelopeParams = TubsParams( _envelopeInnerRadius,
                                               _outerRingOuterRadius,
                                               _innerRingHalfLength );
    }else{
      throw cet::exception("GEOM")
        << "Unknown value of _supportModel in TTrackerMaker::computeDeviceEnvelopeParams "
        << _supportModel
        << "\n";
    }
  }

  // Envelope that holds the full TTracker ("TrackerMother")
  void TTrackerMaker::computeTrackerEnvelope(){

    if ( _supportModel == SupportModel::simple ){

      // Envelope of a single device.
      TubsParams deviceEnvelope = _tt->getDeviceEnvelopeParams();

      // Full length from center to center of the first and last devices.
      double fullLength = _tt->_devices.back().origin().z()-_tt->_devices.front().origin().z();

      // Remember the thickness of the devices.
      double halfLength = fullLength/2. + deviceEnvelope.zHalfLength();

      _tt->_innerTrackerEnvelopeParams = TubsParams( deviceEnvelope.innerRadius(),
                                                     deviceEnvelope.outerRadius(),
                                                     halfLength);

    } else if ( _supportModel == SupportModel::detailedv0 ){

      double fullLength = _tt->_devices.back().origin().z()-_tt->_devices.front().origin().z();
      double halfLength = fullLength/2. + _tt->getDeviceEnvelopeParams().zHalfLength();

      TubsParams val( _envelopeInnerRadius, _outerRingOuterRadius, halfLength);
      _tt->_innerTrackerEnvelopeParams = val;

    } else{

      throw cet::exception("GEOM")
        << "Unknown value of _supportModel in TTrackerMaker::computeTrackerEnvelopeParams "
        << _supportModel
        << "\n";
    }

  } //end TTrackerMaker::computeTrackerEvnvelopeParams


  void TTrackerMaker::recomputeHalfLengths(){

    // This code is only valid for the model detailedv0.
    if ( _supportModel != SupportModel::detailedv0 ) {
      return;
    }

    size_t idev0(0);
    size_t isec0(0);

    deque<Straw>& allStraws              = _tt->_allStraws;
    std::vector<StrawDetail>& allDetails = _tt->_strawDetails;
    size_t originalSize                  = allDetails.size();

    // Step 1: Check that the pattern of _detailIndex is as expected.
    int nBad(0);
    for ( deque<Straw>::const_iterator i=allStraws.begin();
          i != allStraws.end(); ++i ){
      Straw const & straw(*i);
      if ( straw._detailIndex != straw._id.getStraw() ){
        ++nBad;
        cout << "Unexecpted value of detailIndex. Straw "
             << straw.id()
             << " has detailIndex of: "
             << straw._detailIndex
             << endl;
      }
    }
    if ( nBad > 0 ){
      throw cet::exception("GEOM")
        << "TTRackerMaker::recomputeHalfLengths: patterm of _detailIndex is not as expected."
        << "\n";
    }

    // Inner and outer radii of the channel in inner ring.
    double rmin = _innerRingInnerRadius;
    double rmax = rmin + _channelDepth;

    // Number of straws that are too short.
    int nShort(0);

    // Number of straws that have an illegal active length.
    int nIllegalActive(0);

    // Step 2: For all layers in SectorId(0,0) recompute the straw lengths.
    //         For layers > 0:
    //            - create a new StrawDetail object to hold the new length.
    //            - Reseat the _detail and _detailIndex objects in the Straw object.
    vector<Layer>& lays = _tt->_devices.at(idev0)._sectors.at(isec0)._layers;
    for ( size_t ilay=0; ilay<lays.size(); ++ilay){
      Layer& lay(lays.at(ilay));

      for (int ist=0; ist<lay.nStraws(); ++ist ){
        int idx             = lay._indices.at(ist).asInt();
        Straw& straw        = allStraws.at(idx);
        StrawDetail& detail = allDetails.at(straw._detailIndex);
        double r0 = straw.getMidPoint().perp();
        double r1 = r0 - _strawOuterRadius;
        double r2 = r0 + _strawOuterRadius;

        // Choose half length so that the outer edge of straw just touches the outer
        // limit of the channel.
        double hlen = sqrt(diff_of_squares( rmax, r2));

        // Active half-length of the straw.
        double activeHalfLen = sqrt( diff_of_squares(rmin,r0) )-_passivationMargin;
        activeHalfLen = std::max( activeHalfLen, 1.0);

        if ( activeHalfLen > hlen ){
          ++nIllegalActive;
          cout << "TTRackerMaker::recomputeHalfLengths: illegal passivation Margin.\n"
               << "The active length of straw " << ist << " is longer than its physical length\n";
        }

        // Check that the inner edge of the straw reaches the support
        double r3 = sqrt(sum_of_squares(hlen,r1));

        if ( r3 < rmin ){
          ++nShort;
          cout << "Straw is too short to reach the inner edge of thesupport.\n"
               << "Straw; " << straw.id()
               << " Radius at inner corner: " << r3 << " mm\n"
               << "Radius at inner edge of the support ring: " << rmin
               << endl;
        }

        if ( ilay == 0 ){
          detail._halfLength       = hlen;
          detail._activeHalfLength = activeHalfLen;
        } else {
          StrawDetail newDetail       = detail;
          newDetail._halfLength       = hlen;
          newDetail._activeHalfLength = activeHalfLen;
          newDetail._id               = allDetails.size();

          // Enough space has been reserved so that push_back does not invalidate iterators.
          allDetails.push_back(newDetail);
          straw._detail      = &allDetails.back();
          straw._detailIndex = newDetail._id;
        }
      }
    }
    if ( nShort > 0 ){
      throw cet::exception("GEOM")
        << "TTRackerMaker::recomputeHalfLengths: some straaws are too short.\n"
        << "Probably the answer is to deepen the channel."
        << "\n";
    }

    if ( nIllegalActive > 0 ){
      throw cet::exception("GEOM")
        << "TTRackerMaker::recomputeHalfLengths: some straaws have an active length longer than their physical length.\n"
        << "Probably the passivationMargin is too negative."
        << "\n";
    }

    // Step 4: reseat _detail and _detailIndex for all other straws.
    for ( deque<Straw>::iterator i=allStraws.begin();
          i != allStraws.end(); ++i ){
      Straw& straw(*i);

      // These are already done:
      if ( straw.id().getSectorId() == SectorId(0,0) ) continue;
      if ( straw.id().getLayer()    ==            0  ) continue;

      // Get the new detail object for this straw.
      int idx = straw.id().getLayer()*originalSize + straw.id().getStraw();
      StrawDetail const& detail = allDetails.at(idx);

      // Update the info about the detail object.
      straw._detail      = &detail;
      straw._detailIndex = idx;
    }

  } //end TTrackerMaker::recomputeHalfLengths
  double
  TTrackerMaker::sectorRotation(int isec,int idev) const {
    int jdev = idev%2;
    int ista = (idev/2)%2;
    int jsec = isec + jdev*_sectorsPerDevice;
    double phi = _sectorBaseRotations.at(jsec);
    if(ista==1)phi += _oddStationRotation;
    return phi;
  }

} // namespace mu2e
