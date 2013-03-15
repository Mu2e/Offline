// ITracker geometry maker
//
// $Id: ITrackerMaker.hh,v 1.13 2013/03/15 16:20:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 16:20:00 $
//
// Original author G. Tassielli
//

#ifndef ITrackerGeom_ITrackerMaker_hh
#define ITrackerGeom_ITrackerMaker_hh

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>

#include "CLHEP/Vector/ThreeVector.h"
#include "ITrackerGeom/inc/SuperLayer.hh"
#include "ITrackerGeom/inc/SuperLayerInfo.hh"
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
  std::unique_ptr<ITracker> getITrackerPtr() { return std::move(_ltt); }

private:

  void Build();
  void ITFldWireLocater ( boost::shared_ptr<WireDetail> &wdetail, boost::shared_ptr<ITLayer> &itl/*ITLayer *itl*/, int NofWire, double PosRadius, double Theta, double ThetaOffset, double Stereo, double halfAlpha );
  void ITWireLocater ( boost::shared_ptr<WireDetail> &wdetaill, Wire::Wtype wireType, boost::shared_ptr<ITLayer> &itl/*ITLayer *itl*/, int NofWire, double PosRadius, double Theta, double ThetaOffset, double Stereo, double halfAlpha, int copyNunOffset=0, boost::shared_ptr<CellDetail> *celldetail = nullptr );

  void AssignFieldWireToCell();

  // Name of external gdml geometry file description.
  std::string _extFile;
  std::string _extWireFile;
  bool _isExternal;

  // Basic geometry element parameter
  int _nSWire;                  //Number of sense wire at the first Ring of the first SuperLayer
  int _nSDeltaWire;             //Increment of the sense wire number for each SuperLayer
  int _nSuperLayer;             //Number of SuperLayer
  int _nRing;                   //Number of Ring in each Superlayer
  int _nVerticalFWire;          //Number of additional Vertical Field Wire in the case of Square cells
  int _StoFWireRatio;           //Sense to Field wire ratio per cell in the case of Square cells
  double _cellDimension;        //Cell dimension in the case of Square cells
  double _FWireStep;            //Field Wire step distance in the case of Square cells
  double _r0;                   //Nominal Inner radius of the tracker
  double _halfLength;           //Nominal Half-Length of the tracker in the barrel region
  double _rOut;                 //Nominal Outer radius of the tracker
  double _drop;                 //Drop distance of the wires (in case of constant drop version (v41) the wires stereo angles depend by it)
  double _alpha;                //alpha of the wires (in case of constant alpha version (v42) the wires stereo angles and their drops depend by it)
  bool   _isDumbbell;           //true=dumbbell option, false=the wires are unique
  std::vector<double> _zZones;  //limits in z of the active wire zones in the case of dumbbell option

  int _geomType;                //Cell Geometry type: 2 hexagonal, 3 square
  int _endCapType;              //EndCap shape type: 0 plane, 1 spherical
  double _voxFactor;            //voxelization optimization factor
  bool _notExtVoxel;
  bool _displayGasLayer;        //Allow to display the gas inside the chamber
  bool _displayWires;           //Allow to display every wires inside gas inside the chamber.

  double _z0;                   //Shift along z of the center of the tracker

  bool _detailedWireSupport;
  bool _isElectCont;
  double _elctContWallThick;

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

  //Guard Wire dimensions and materials composition
  int _nInGuardWires;
  int _nOutGuardWires;
  double _inGuardRad;
  double _outGuardRad;

  double _inGWireDiameter;
  double _outGWireDiameter;
  // Names of the materials of the wires.
  std::vector<std::string> _inGwMaterialsName;
  std::vector<std::string> _outGwMaterialsName;

  std::vector<double> _inGwShellsThicknesses;
  std::vector<double> _outGwShellsThicknesses;

  //Walls descriptions
  std::multimap<Wall::Walltype,Wall* > _walls;

  // Center of the tracker.
  CLHEP::Hep3Vector _center;

  std::unique_ptr<ITracker> _ltt;

};

}  //namespace mu2e

#endif /* ITrackerGeom_ITrackerMaker_hh */
