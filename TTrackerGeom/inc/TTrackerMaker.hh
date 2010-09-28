#ifndef TTRACKERMAKER_HH
#define TTRACKERMAKER_HH
//
// Construct and return a TTracker.
//
// $Id: TTrackerMaker.hh,v 1.6 2010/09/28 21:43:04 genser Exp $
// $Author: genser $
// $Date: 2010/09/28 21:43:04 $
//
// Original author Rob Kutschke
//

#include <vector>
#include <memory>
#include <string>

//#include "TTrackerGeom/inc/TLayerInfo.hh"

#include "TrackerGeom/inc/Device.hh"
#include "TrackerGeom/inc/Sector.hh"
#include "TrackerGeom/inc/Layer.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class TTracker;
class SimpleConfig;

class TTrackerMaker{

public:

  TTrackerMaker( SimpleConfig const& config );  
  ~TTrackerMaker ();

  // This is the accessor that will remain.
  std::auto_ptr<TTracker> getTTrackerPtr() { return _tt; }

private:

  // Extract info from the config file.
  void parseConfig( const SimpleConfig& config );

  // Fill the details of the different straw types.
  void makeDetails();

  void makeDevice( DeviceId devId );
  void makeSector( const SectorId& secId, Device& dev );
  void makeLayer ( const LayerId& layId,  Sector& sec );
  void makeManifolds( const SectorId& secId);

  void computeStrawHalfLengths();
  void computeSectorBoxParams(Sector& sector, Device& dev);
  void computeConstantSectorBoxParams();

  // Do the work of constructing it.
  void buildIt();

  double chooseDeviceRotation( int idev ) const;
  double chooseDeviceSpacing( int idev ) const;
  double findFirstDevZ0() const;

  // Basic parameters needed to describe the TTracker.
  int    _numDevices;                  // Number of devices.
  int    _sectorsPerDevice;            // Number of sectors in one device.
  int    _layersPerSector;             // Number of layers in one sector.
  int    _manifoldsPerEnd;             // Number of manifolds along one end of the wires in a layer.
  int    _strawsPerManifold;           // Number of straws connected to each manifold.
  int    _rotationPattern;             // Pattern of rotations from device to device.
  int    _spacingPattern;              // Pattern of spacing from device to device.
  double _zCenter;                     // Position of the center of the tracker, in the Mu2e coord system.
  double _envelopeInnerRadius;         // Inner radius of inside of innermost straw.
  double _strawOuterRadius;            // Radius of each straw.
  double _strawWallThickness;          // Thickness of each straw.
  double _strawGap;                    // Gap between straws.
  double _deviceSpacing;               // Z-separation between adjacent stations.
  double _deviceHalfSeparation;        // Z-separation between adjacent devices.
  double _deviceRotation;              // Relative rotation of each succesive device.
  double _innerSupportRadius;          // Inner radius of support frame.
  double _outerSupportRadius;          // Outer radius of support frame.
  double _supportHalfThickness;        // Thickness of support frame.
  double _wireRadius;                  // Wire radius
  double _manifoldYOffset;             // Offset of manifold from inner edge of support.
  std::string _envelopeMaterial;       // Material for the envelope volume.

  std::vector<double> _manifoldHalfLengths; // Dimensions of each manifold.
  std::vector<std::string> _strawMaterials; // Names of the materials.

  // Base rotations of a sector; does not include device rotation.
  std::vector<double> _sectorBaseRotations;
  std::vector<double> _sectorZSide;

  std::auto_ptr<TTracker> _tt;

  // Derived parameters.

  // Compute at start and compare against final result.
  int _nStrawsToReserve;

  // Lengths of straws indexed by manifold, from innermost radius, outwards.
  std::vector<double> _strawHalfLengths;

  // sector box half lengths
  std::vector<double> _sectorBoxHalfLengths;

  // Z Location of the first device.
  double _z0;

  // global straw counter (will be used per sector)
  int _istraw;

};

}  //namespace mu2e

#endif
