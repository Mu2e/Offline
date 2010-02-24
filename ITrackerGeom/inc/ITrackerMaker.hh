#ifndef ITRACKERMAKER_HH
#define ITRACKERMAKER_HH

#include <vector>
#include <memory>
#include <string>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "CLHEP/Vector/ThreeVector.h"
#include "ITrackerGeom/inc/SuperLayerInfo.hh"
#include "ITrackerGeom/inc/SuperLayer.hh"
#include "ITrackerGeom/inc/Wall.hh"

namespace mu2e {

class ITracker;
class SimpleConfig;
//class Wall;

class ITrackerMaker{

public:

  ITrackerMaker( SimpleConfig const& config );  
  ~ITrackerMaker ();

  // This is the accessor that will remain.
  std::auto_ptr<ITracker> getITrackerPtr() { return _ltt; }

private:

  void Build();
  void ITFldWireLocater ( boost::shared_ptr<WireDetail> &wdetail, boost::shared_ptr<ITLayer> &itl/*ITLayer *itl*/, int NofWire, double PosRadius, double Theta, double ThetaOffset, double Stereo, double halfAlpha );
  void ITWireLocater ( boost::shared_ptr<WireDetail> &wdetaill, Wire::Wtype wireType, boost::shared_ptr<ITLayer> &itl/*ITLayer *itl*/, int NofWire, double PosRadius, double Theta, double ThetaOffset, double Stereo, double halfAlpha, int copyNunOffset=0, boost::shared_ptr<CellDetail> *celldetail = NULL );

  // Name of external gdml geometry file description.
  std::string _extFile;
  std::string _extWireFile;
  bool _isExternal;

  // Basic geometry element parameter
  int _nSWire;				//Number of sense wire at the first Ring of the first SuperLayer
  int _nSDeltaWire;			//Increment of the sense wire number for each SuperLayer
  int _nSuperLayer;			//Number of SuperLayer
  int _nRing;				//Number of Ring in each Superlayer
  int _nVerticalFWire;		//Number of additional Vertical Field Wire in the case of Square cells
  int _StoFWireRatio;		//Sense to Field wire ratio per cell in the case of Square cells
  double _cellDimension;	//Cell dimension in the case of Square cells
  double _FWireStep;		//Field Wire step distance in the case of Square cells
  double _r0;				//Nominal Inner radius of the tracker
  double _halfLength;		//Nominal Half-Length of the tracker in the barrel region
  double _rOut;				//Nominal Outer radius of the tracker
  double _drop;				//Drop distance of the wires (the wires stereo angles depend by it)

  int _geomType;			//Cell Geometry type: 2 hexagonal, 3 square
  int _endCapType;			//EndCap shape type: 0 plane, 1 spherical
  double _voxFactor;		//voxelization optimization factor
  bool _notExtVoxel;
  bool _displayGasLayer;	//Allow to display the gas inside the chamber
  bool _displayWires;		//Allow to display every wires inside gas inside the chamber.

  double _z0;				//Shift along z of the center of the tracker

  // Number of layers and of cells per layer in superlayer.
  std::vector<SuperLayerInfo> _slayersInfo;

  //Active gas material
  std::string _fillMaterial;

  //Wire dimensions and materials composition
  double _fWireDiameter;
  double _sWireDiameter;
  // Names of the materials of the wires.
  std::vector<std::string> _fwMaterialsName;
  std::vector<std::string> _swMaterialsName;

  std::vector<double> _fwShellsThicknesses;
  std::vector<double> _swShellsThicknesses;

  //Walls descriptions
  std::multimap<Wall::Walltype,Wall* > _walls;

  // Center of the tracker.
  CLHEP::Hep3Vector _center;

  std::auto_ptr<ITracker> _ltt;

};

}  //namespace mu2e

#endif /*ITRACKERMAKER_HH*/
