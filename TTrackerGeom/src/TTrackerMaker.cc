//
// Construct and return an TTracker.
//
//
// $Id: TTrackerMaker.cc,v 1.5 2010/07/29 21:02:42 genser Exp $
// $Author: genser $
// $Date: 2010/07/29 21:02:42 $
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
    _spacingPattern     = config.getInt("ttracker.spacingPattern");

    _zCenter              = config.getDouble("ttracker.z0")*CLHEP::mm;
    _envelopeInnerRadius  = config.getDouble("ttracker.envelopeInnerRadius")*CLHEP::mm;
    _strawOuterRadius     = config.getDouble("ttracker.strawOuterRadius")*CLHEP::mm;
    _strawWallThickness   = config.getDouble("ttracker.strawWallThickness")*CLHEP::mm;
    _deviceSeparation     = config.getDouble("ttracker.deviceSeparation")*CLHEP::mm;
    _deviceSpacing        = config.getDouble("ttracker.deviceSpacing")*CLHEP::mm;
    _deviceRotation       = config.getDouble("ttracker.deviceRotation")*CLHEP::degree;

    _outerSupportRadius   = config.getDouble("ttracker.outerSupportRadius")*CLHEP::mm;
    _innerSupportRadius   = config.getDouble("ttracker.innerSupportRadius")*CLHEP::mm;
    _supportHalfThickness = config.getDouble("ttracker.supportHalfThickness")*CLHEP::mm;

    _wireRadius           = config.getDouble("ttracker.wireRadius")*CLHEP::mm;
    
    _manifoldYOffset      = config.getDouble("ttracker.manifoldYOffset")*CLHEP::mm;
    config.getVectorDouble("ttracker.manifoldHalfLengths", _manifoldHalfLengths, 3);
    for ( int i=0; i<_manifoldHalfLengths.size(); ++i ){
      _manifoldHalfLengths.at(i) *= CLHEP::mm;
    }
    
    config.getVectorString("ttracker.strawMaterials", _strawMaterials, 3);

    _envelopeMaterial = config.getString("ttracker.mat.vacuum");
    
    //string ttracker.mat.manifold  = "G4_Al";  // Placeholder.
    //string ttracker.mat.support   = "G4_Al";  // Placeholder.

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

    else
      throw cms::Exception("GEOM")
        << "Unrecognized rotation pattern in TTrackerMaker. \n";



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
    _z0 = -findFirstDevZ0();

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

    //_tt->forAllLayers( lptest);
    //_tt->forAllDevices( devtest);

  }

  void TTrackerMaker::makeDevice( DeviceId devId ){

    int idev = devId;

    double devDeltaZ = chooseDeviceSpacing(idev);
    CLHEP::Hep3Vector origin( 0., 0., _z0+devDeltaZ);
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

    makeManifolds( secId );

    for ( int ilay=0; ilay<_layersPerSector; ++ilay ){
      makeLayer( LayerId(secId,ilay), sector );
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

    // Start to populate the layer.
    layer._nStraws = _manifoldsPerEnd*_strawsPerManifold;
    layer._straws.reserve(_manifoldsPerEnd*_strawsPerManifold);

    // Space between first/last straw and edge of manifold.
    double dx = _manifoldHalfLengths.at(0) - _strawOuterRadius*_strawsPerManifold;

    // |z| of straw center, relative to the center of the device.
    // Sign is taken care of elsewhere.
    double zOffset = _supportHalfThickness + _manifoldHalfLengths.at(2);

    // Rotation that puts wire direction and wire mid-point into their
    // correct orientations.
    // CLHEP::HepRotationZ RZ(_sectorBaseRotations.at(isec));
    CLHEP::HepRotationZ RZ(_sectorBaseRotations.at(isec) + device.rotation());

    // Unit vector in the wire direction. (nominal is the sector 0 to the right?)
    CLHEP::Hep3Vector unit = RZ*CLHEP::Hep3Vector(0.,1.,0.);

    // Straw number within the layer; does not reset to zero at each manifold.
    int istraw(-1);

    // Add all of the straws
    for ( int iman=0; iman<_manifoldsPerEnd; ++iman ){

      // Inner edge of the innermost wire connected to this manifold.
      double xA = _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*iman - dx;

      // Add all of the straws connected to this manifold.
      for ( int istr=0; istr<_strawsPerManifold; ++istr ){
        ++istraw;

        // Construct straw midpoint in its base position in the 
        // coord system of the device envelope.
        double xstraw = xA + (1. + 2.*istr)*_strawOuterRadius ;
        CLHEP::Hep3Vector mid( xstraw, 0., zOffset*_sectorZSide.at(isec) );
        mid += device.origin();
        
        // Rotate straw midpoint to its actual location.
        CLHEP::Hep3Vector offset = RZ*mid;

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
        

//           if ( layId.getDevice() != -1 ){
//           cout << "Position: " << setw(3) <<
//           layId << " | " << setw(3) <<
//           iman << " " << setw(3) <<
//           istr                << " | " << setw(3) <<
//           istraw << " " << fixed << setprecision(2) << setw(8) <<
//           xstraw << " " << fixed << setprecision(2) << setw(8) <<
//           2.*_strawHalfLengths.at(iman) << " " << fixed << setprecision(2) <<
//           mid << " " << fixed << setprecision(2) << setw(8) <<
//           device.origin() << " | "  << setw(3) <<
//           index <<  " " << setw(3) <<
//           allStraws.size() << " "  << setw(3) <<
//           layer._straws.size() << " "
//           << endl;
//           }

      }
    }
    
  }

  void TTrackerMaker::makeManifolds( const SectorId& secId){

    if ( _sectorsPerDevice != 4 && _sectorsPerDevice != 6 ) {
      throw cms::Exception("GEOM")
        << "This code only knows how to do 4 or 6 sectors per device.\n";
    }

    double phi = _tt->getDevice(secId.getDevice()).rotation();
    CLHEP::HepRotationZ RZ(phi);


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
//         _tt->getDevice(secId.getDevice()).Id() << ", " <<
//         secId.getSector() << ", " <<
//         origin << ", " << _manifoldHalfLengths.at(0) <<endl;


      _tt->_allManifolds.push_back( Manifold( origin, _manifoldHalfLengths) );
    }

  }

  // Assumes all devices and all sectors are the same.
  // The straw length depends only on the manifold number.
  // See Mu2e-doc-??? for the algorithm.
  void TTrackerMaker::computeStrawHalfLengths(){
    _strawHalfLengths.clear();

    // Should get this from somewhere else.
    static const double tolerance = 1.e-4*CLHEP::mm;

    // Space between first/last straw and edge of manifold.
    double dx = _manifoldHalfLengths.at(0) - _strawOuterRadius*_strawsPerManifold;

    // Check for a legal size.
    if ( dx < 0.0 ){
      if ( dx < -tolerance ){
        throw cms::Exception("GEOM")
          << "Manifolds are too small to hold straws!\n";
      }
      dx = 0;
    }

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      double xA = _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*i - dx;
      double yA = sqrt( square(_innerSupportRadius) - square(xA) );
      double yB = yA + _manifoldYOffset;
      // cout << "Straw Length: " << xA << " " << yB*2.0 << endl;
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
            _strawHalfLengths.at(i),
            _wireRadius
            )
          );
    }
    
  }

  void  TTrackerMaker::computeSectorBoxParams(Sector& sector, Device& dev){

    // get sector number

    int isec = sector.Id().getSector();

    //  void  TTrackerMaker::computeSectorBoxParams(vector<double>& boxHalfLengths){

    // the box is a trapezoid 
    // shorter x             is the length of the straws in the top/last (shortest) manifold
    // longer  x is longer than the length of the straws in the longest manifold 
    // the other dimentions are "z", the combined manifold width
    // the "thickness" of the trpezoid y, the layer thickness * number of layers 


    // z
    double bz = _manifoldHalfLengths.at(0)*double(_manifoldsPerEnd);

    // calculating "longer" x;
    // starting from the tng of the slope 

    // calculate the largest slope starting from the longest straws manifold
    double maxtg = 0.0;

    for (size_t i=1; i!=_manifoldsPerEnd; ++i) {
      double ttg = (_manifoldHalfLengths.at(0)*double(i) ) /
        ( _strawHalfLengths.at(0) - _strawHalfLengths.at(i) );
      if (maxtg < ttg ) {
        maxtg = ttg;
      }
    }

    double slopeAlpha = atan(maxtg);

    // finally x (long, short)
    double bxl = _manifoldHalfLengths.at(0)/tan(slopeAlpha)+_strawHalfLengths.at(0);
    double bxs = bxl - bz/tan(slopeAlpha);

    // we need to make sure the trapezoids do not extend beyond the device envelope...
    // we may need to do "volume subtraction"

    // y
    // double by = _layersPerSector*_tt->_strawDetails.at(0).outerRadius();
    // manifold better be thicker than the straws:
    double by = _layersPerSector* max(_manifoldHalfLengths.at(2),_tt->_strawDetails.at(0).outerRadius());

    // now push it all back into the vector

    // Pad the box to be slightly larger than it needs to be
    const double pad = 0.001; // this needs to be in the geom file...

    sector._boxHalfLengths.reserve(5);

    // the order is forced by the LTracker and nestTrp/G4Trd and Sector data

    sector._boxHalfLengths.push_back(0.0); //dummy to be compatible with LTracker
    sector._boxHalfLengths.push_back(bz+pad);
    sector._boxHalfLengths.push_back(by+pad);

    sector._boxHalfLengths.push_back(bxs+pad);
    sector._boxHalfLengths.push_back(bxl+pad);

    // Now calculate the rotations and placement of the sector envelope

    sector._boxRxAngle = 0.;
    sector._boxRyAngle = 0.;
    sector._boxRzAngle = _sectorBaseRotations.at(isec) + dev.rotation();

    double zOffset = _layersPerSector*(_supportHalfThickness + _manifoldHalfLengths.at(2))*
      _sectorZSide.at(isec);
    
    double xOffset = _envelopeInnerRadius + _manifoldHalfLengths.at(0)*_manifoldsPerEnd;

    CLHEP::HepRotationZ RZ(sector._boxRzAngle);

    sector._boxOffset  = RZ*(CLHEP::Hep3Vector( xOffset, 0., zOffset) + dev.origin());

    // we may want to set to 0.0 some values smaller than say 10-6...
    const double max0val = 1.e-06;
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
      } else 
	return _deviceRotation;
    }
    else 
      throw cms::Exception("GEOM")
	<< "Unrecognized rotation pattern in TTrackerMaker. \n";
  
  }


  // Compute the spacing for the given device.
  double TTrackerMaker::chooseDeviceSpacing( int idev ) const{


    if ( _spacingPattern == 0 ) {
      return idev * _deviceSeparation;
    }

    else if ( _spacingPattern == 1 ) {
      int k = idev%2;
      if (k == 0 ) {
	return  ( idev/2 )  * ( _deviceSeparation + _deviceSpacing);
      } else if (k == 1 ) {
	return ( (idev-1)/2 ) * ( _deviceSeparation + _deviceSpacing ) + _deviceSpacing;
      }
    }
    else 
      throw cms::Exception("GEOM")
	<< "Unrecognized separation pattern in TTrackerMaker. \n";
  
  }

  double TTrackerMaker::findFirstDevZ0() const{

    if ( _spacingPattern == 0 ) {
      return _deviceSeparation*double(_numDevices-1)/2.0;
    }

    else if ( _spacingPattern == 1 ) {
      return (_deviceSeparation+_deviceSpacing) *double(_numDevices-1)/4.0;
    }
    else 
      throw cms::Exception("GEOM")
	<< "Unrecognized separation pattern in TTrackerMaker. \n";
  }



} // namespace mu2e
