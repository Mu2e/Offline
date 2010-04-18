//
// Construct and return an TTracker.
//
//
// $Id: TTrackerMaker.cc,v 1.1 2010/04/18 00:37:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:37:16 $
//
// Original author Rob Kutschke
//

#include <iostream>
#include <iomanip>
#include <cmath>


// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "TTrackerGeom/inc/TTrackerMaker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Sector.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "Mu2eUtilities/inc/for_all.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"
#include "GeneralUtilities/inc/pow.hh"

#ifndef __CINT__ 

using namespace CLHEP;

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  TTrackerMaker::TTrackerMaker( SimpleConfig const& config){
    
    parseConfig(config);

    buildIt( );
  }
  
  TTrackerMaker::~TTrackerMaker (){}

  void TTrackerMaker::parseConfig( const SimpleConfig& config ){

    _numDevices         = config.getInt("ttracker.numDevices");
    _sectorsPerDevice   = config.getInt("ttracker.sectorsPerDevice");
    _layersPerSector    = config.getInt("ttracker.layersPerSector");
    _manifoldsPerEnd    = config.getInt("ttracker.manifoldsPerEnd");
    _strawsPerManifold  = config.getInt("ttracker.strawsPerManifold");
    _rotationPattern    = config.getInt("ttracker.rotationPattern");

    _zCenter              = config.getDouble("ttracker.z0")*mm;
    _envelopeInnerRadius  = config.getDouble("ttracker.envelopeInnerRadius")*mm;
    _strawOuterRadius     = config.getDouble("ttracker.strawOuterRadius")*mm;
    _strawWallThickness   = config.getDouble("ttracker.strawWallThickness")*mm;
    _deviceSeparation     = config.getDouble("ttracker.deviceSeparation")*mm;
    _deviceRotation       = config.getDouble("ttracker.deviceRotation")*degree;

    _outerSupportRadius   = config.getDouble("ttracker.outerSupportRadius")*mm;
    _innerSupportRadius   = config.getDouble("ttracker.innerSupportRadius")*mm;
    _supportHalfThickness = config.getDouble("ttracker.supportHalfThickness")*mm;

    _wireRadius           = config.getDouble("ttracker.wireRadius")*mm;
    
    _manifoldYOffset      = config.getDouble("ttracker.manifoldYOffset")*mm;
    config.getVectorDouble("ttracker.manifoldHalfLengths", _manifoldHalfLengths, 3);
    for ( int i=0; i<_manifoldHalfLengths.size(); ++i ){
      _manifoldHalfLengths[i] *= mm;
    }
    
    config.getVectorString("ttracker.strawMaterials", _strawMaterials, 3);

    _envelopeMaterial = config.getString("ttracker.mat.vacuum");
    
    //string ttracker.mat.manifold  = "G4_Al";  // Placeholder.
    //string ttracker.mat.support   = "G4_Al";  // Placeholder.

    // Also define some parameters that may become variable some day.
    _sectorBaseRotations.clear();
    _sectorBaseRotations.push_back(   0.*degree);
    _sectorBaseRotations.push_back(  90.*degree);
    _sectorBaseRotations.push_back( 180.*degree);
    _sectorBaseRotations.push_back( 270.*degree);
    _sectorZSide.clear();
    _sectorZSide.push_back(-1.);
    _sectorZSide.push_back(+1.);
    _sectorZSide.push_back(-1.);
    _sectorZSide.push_back(+1.);

  }

  void lptest( const Layer& lay){
    cout << lay.Id() << " |  " 
         << lay.nStraws()  <<  " |  "
         << lay.getStraws().capacity() << " "
         << endl;
  }
  
  void devtest( const Device& dev){
    cout << "Device: "
         << dev.Id() << " "
         << dev.origin() << " "
         << dev.rotation()
         << endl;
  }

  void TTrackerMaker::buildIt(){

    // Make an empty TTracker.
    _tt = auto_ptr<TTracker>(new TTracker());

    makeDetails();

    // Fill the information about the supports.
    _tt->_supportParams = Support( _innerSupportRadius,
                                   _outerSupportRadius,
                                   _supportHalfThickness,
                                   "CarbonFiber");

    _tt->_z0                  = _zCenter;
    _tt->_envelopeInnerRadius = _envelopeInnerRadius;
    _tt->_manifoldHalfLengths = _manifoldHalfLengths;
    _tt->_envelopeMaterial    = _envelopeMaterial;

    // Z location of the first device.
    _z0 = -_deviceSeparation*double(_numDevices-1)/2.0;

    // Reserve space for straws so that pointers are valid.
    _nStrawsToReserve = _numDevices * _sectorsPerDevice * _layersPerSector * 
      _manifoldsPerEnd * _strawsPerManifold;
    //_tt->_allStraws.reserve(_nStrawsToReserve);

    _tt->_devices.reserve(_numDevices);
    // Construct the devices and their internals.
    for ( int idev=0; idev<_numDevices; ++idev ){
      makeDevice( DeviceId(idev) );
    }

    // Fill all of the non-persistent information.
    _tt->fillPointers();

    //_tt->forAllLayers( lptest);
    //_tt->forAllDevices( devtest);

  }

  void TTrackerMaker::makeDevice( DeviceId devId ){

    int idev = devId;

    Hep3Vector origin( 0., 0., _z0+_deviceSeparation*idev);
    double phi = chooseDeviceRotation(idev)+_deviceRotation;

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

    makeManifolds( secId );

    for ( int ilay=0; ilay<_layersPerSector; ++ilay ){
      makeLayer( LayerId(secId,ilay), sector );
    }

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

    // Start to populate the layer.
    layer._nStraws = _manifoldsPerEnd*_strawsPerManifold;
    layer._straws.reserve(_manifoldsPerEnd*_strawsPerManifold);

    // Space between first/last straw and edge of manifold.
    double dx = _manifoldHalfLengths[0] - _strawOuterRadius*_strawsPerManifold;

    // |z| of straw center, relative to the center of the device.
    // Sign is taken care of elsewhere.
    double zOffset = _supportHalfThickness + _manifoldHalfLengths[2];

    // Rotation that puts wire direction and wire mid-point into their
    // correct orientations.
    HepRotationZ RZ(_sectorBaseRotations.at(isec));

    // Unit vector in the wire direction.
    Hep3Vector unit = RZ*Hep3Vector(0.,1.,0.);

    // Straw number within the layer; does not reset to zero at each manifold.
    int istraw(-1);

    // Add all of the straws
    for ( int iman=0; iman<_manifoldsPerEnd; ++iman ){

      // Inner edge of the innermost wire connected to this manifold.
      double xA = _envelopeInnerRadius + 2.*_manifoldHalfLengths[0]*iman - dx;

      // Add all of the straws connected to this manifold.
      for ( int istr=0; istr<_strawsPerManifold; ++istr ){
        ++istraw;

        // Construct straw midpoint in its base position in the 
        // coord system of the device envelope.
        double xstraw = xA + (1. + 2.*istr)*_strawOuterRadius ;
        Hep3Vector mid( xstraw, 0., zOffset*_sectorZSide.at(isec) );
        mid += device.origin();
        
        // Rotate straw midpoint to its actual location.
        Hep3Vector offset = RZ*mid;

        StrawIndex index(allStraws.size());

        allStraws.push_back( Straw( StrawId( layId, istraw),
                                    index,
                                    offset,
                                    &_tt->_strawDetails.at(iman),
                                    iman,
                                    unit
                                    )
                             );
        layer._straws.push_back(&allStraws.back());
        layer._indices.push_back(index);
        

        /*
        if ( layId.getDevice() == 0 ){
          cout << "Position: "
               << layId << " | "
               << iman << " "
               << istr                << " | "
               << istraw << " "
               << xstraw  << " "
               << 2.*_strawHalfLengths.at(iman) << " "
               << mid << " "
               << device.origin() << " | " 
               << index <<  " "
               << allStraws.size() << " " 
               << layer._straws.size() << " "
               << endl;
        }
        */
      }
    }
    
  }

  void TTrackerMaker::makeManifolds( const SectorId& secId){

    if ( _sectorsPerDevice != 4 ) {
      throw cms::Exception("GEOM")
        << "This code only knows how to do 4 sectors per device.\n";
    }

    double phi = _tt->getDevice(secId.getDevice()).rotation();
    HepRotationZ RZ(phi);


    for ( int i=0; i<_manifoldsPerEnd; ++i){
    
      // First compute everything in their nominal positions: sector 2, top
      double x0 = _envelopeInnerRadius + 
	_strawsPerManifold*_strawOuterRadius +
	_manifoldHalfLengths[0];
    
      double y0 = _tt->_strawDetails[i].halfLength() + _manifoldHalfLengths[2];

      double z0 = ( _supportHalfThickness + _manifoldHalfLengths[2] );
      if ( secId.getSector() <= 1 ) z0 = -z0;
      
      Hep3Vector origin(x0,y0,z0);


      _tt->_allManifolds.push_back( Manifold( origin, _manifoldHalfLengths) );
    }

  }

  // Assumes all devices and all sectors are the same.
  // The straw length depends only on the manifold number.
  // See Mu2e-doc-??? for the algorithm.
  void TTrackerMaker::computeStrawHalfLengths(){
    _strawHalfLengths.clear();

    // Should get this from somewhere else.
    static const double tolerance = 1.e-4*mm;

    // Space between first/last straw and edge of manifold.
    double dx = _manifoldHalfLengths[0] - _strawOuterRadius*_strawsPerManifold;

    // Check for a legal size.
    if ( dx < 0.0 ){
      if ( dx < -tolerance ){
        throw cms::Exception("GEOM")
          << "Manifolds are too small to hold straws!\n";
      }
      dx = 0;
    }

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      double xA = _envelopeInnerRadius + 2.*_manifoldHalfLengths[0]*i - dx;
      double yA = sqrt( square(_innerSupportRadius) - square(xA) );
      double yB = yA + _manifoldYOffset;
      //      cout << "Straw Length: " << xA << " " << yB*2.0 << endl;
      _strawHalfLengths.push_back(yB);
    }
  }
  
  void TTrackerMaker::makeDetails(){

    computeStrawHalfLengths();

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      _tt->_strawDetails.push_back
	( StrawDetail
	  ( i,
	    _strawMaterials,
	    _strawOuterRadius,
	    _strawWallThickness,
	    _strawHalfLengths[i],
	    _wireRadius
	    )
	  );
    }
    
  }

  // Compute the rotation for the given device.
  // For now, only one pattern is known.
  double TTrackerMaker::chooseDeviceRotation( int idev ) const{
    if ( _rotationPattern != 0 ){
      throw cms::Exception("GEOM")
        << "Unrecognized rotation pattern in TTrackerMaker. \n";
    }
    int k = idev%3;
    if ( k == 0 ) {
      return 0.;
    } else if ( k == 1 ){
      return _deviceRotation;
    } else{
      return -_deviceRotation;
    }
  }

} // namespace mu2e

#endif
