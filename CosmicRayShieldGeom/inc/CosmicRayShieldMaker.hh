#ifndef CosmicRayShieldGeom_CosmicRayShieldMaker_hh
#define CosmicRayShieldGeom_CosmicRayShieldMaker_hh
//
// Class to construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.hh,v 1.10 2012/06/07 04:55:00 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2012/06/07 04:55:00 $
//
// Original author KLG
//
 
#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CosmicRayShield;
  class SimpleConfig;

  class CRSScintillatorShield;
  class CRSScintillatorModule;
  class CRSScintillatorLayer;

class CosmicRayShieldMaker {

public:

  CosmicRayShieldMaker( SimpleConfig const & config, double solenoidOffset );

  void parseConfig( SimpleConfig const & _config );

  void makeCRSSteelShield(SimpleConfig const & _config);

  void calculateCommonCRSScintillatorParameters();
  void makeDetails();
  void makeShields();
  void makeModules(CRSScintillatorShield& shield);
  void makeLayers(CRSScintillatorModule& module);
  void makeBars(CRSScintillatorLayer& layer);


  // This is deprecated and will go away soon.
  // Still needed for root graphics version.
  const CosmicRayShield& getCosmicRayShield() const { return *_crs;}

  // This is the accessor that will remain.
  std::auto_ptr<CosmicRayShield> getCosmicRayShieldPtr() { return _crs; }

private:

  std::auto_ptr<CosmicRayShield> _crs;

  int _diagLevel;

  bool _hasPassiveShield;
  bool _hasActiveShield;

  double _HallSteelHalfThick;
  double _HallSteelHalfLengthX;
  double _HallSteelHalfLengthY;
  double _HallSteelHalfLengthZ;
  double _HallSteelHoleRadius;

  CLHEP::Hep3Vector _HallSteelOffset; // aka CRPassiveShield offset

  CLHEP::Hep3Vector _TopHallSteelOffset;
  CLHEP::Hep3Vector _LeftHallSteelOffset;
  CLHEP::Hep3Vector _RightHallSteelOffset;
  CLHEP::Hep3Vector _DownstreamHallSteelOffset;

  std::string _HallSteelMaterialName;

  int                 _scintillatorLayersPerModule;
  int                 _scintillatorBarsPerFullLayer;
  std::vector<double> _scintillatorBarHalfLengths;
  std::string         _scintillatorBarMaterialName;
  double              _scintillatorLayerShift;
  double              _scintillatorLayerGap;
  double              _scintillatorBarPigmentationHalfThickness;
  std::string         _scintillatorBarPigmentationMaterialName;
  std::vector<double> _scintillatorModuleOuterSheetHalfLengths;
  std::string         _scintillatorModuleOuterSheetMaterialName;
  std::string         _scintillatorModuleInterLayerSheetMaterialName;
  double              _scintillatorModuleInterLayerSheetHalfThickness;
  double              _scintillatorOverlap;
  std::vector<double> _moduleUnistrutHalfLengths;
  //  double              _wallUnistrutHalfThickness; // was used to displace the CRV

  std::vector<int>    _shieldR_NumberOfModules;
  std::vector<int>    _shieldL_NumberOfModules;
  std::vector<int>    _shieldD_NumberOfModules;
  std::vector<int>    _shieldT_NumberOfModules;
  std::vector<int>    _shieldTSR_NumberOfModules;
  std::vector<int>    _shieldTSL_NumberOfModules;
  std::vector<int>    _shieldTST_NumberOfModules;

  CLHEP::Hep3Vector   _scintillatorShieldOffset; // aka CRActiveShield offset

  std::vector<double> _shieldR_Offset;
  std::vector<double> _shieldL_Offset;
  std::vector<double> _shieldD_Offset;
  std::vector<double> _shieldT_Offset;
  std::vector<double> _shieldTSR_Offset;
  std::vector<double> _shieldTSL_Offset;
  std::vector<double> _shieldTST_Offset;

  // derived quantities etc...

  int    _totalNumberOfBars;

  double _scintillatorFullLayerHalfWidth;
  double _scintillatorHalfLayerHalfWidth;
  double _scintillatorLayerHalfLength;

  double _scintillatorModuleCoreHalfThickness; // no unistruts

  double _scintillatorFullModuleHalfWidth;
  double _scintillatorHalfModuleHalfWidth;

  double _scintillatorModuleHalfThickness;

  double _scintillatorModuleHalfLength;

  double _scintillatorModuleOverlap;

  double _scintillatorShieldHalfThickness;

// double _scintillatorShieldOffsetToTheSideOfHallSteel; // was calculated based on _wallUnistrutHalfThickness

};

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CosmicRayShieldMaker_hh */
