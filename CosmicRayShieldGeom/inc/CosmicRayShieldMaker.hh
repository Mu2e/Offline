#ifndef CosmicRayShieldGeom_CosmicRayShieldMaker_hh
#define CosmicRayShieldGeom_CosmicRayShieldMaker_hh
//
// Class to construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.hh,v 1.5 2011/05/19 21:36:23 wb Exp $
// $Author: wb $
// $Date: 2011/05/19 21:36:23 $
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

  CosmicRayShieldMaker( SimpleConfig const & config );

  void parseConfig( SimpleConfig const & _config );

  void makeCRSSteelShield(SimpleConfig const & _config);

  void calculateCRSOffsets(SimpleConfig const & _config);

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

  double _HallSteelHalfThick;
  double _HallSteelHalfLengthXY;
  double _HallSteelHalfLengthZ;
  double _HallSteelHoleRadius;

  std::vector<double> _HallSteelOffset;

  CLHEP::Hep3Vector _TopHallSteelOffset;
  CLHEP::Hep3Vector _BottomHallSteelOffset;
  CLHEP::Hep3Vector _LeftHallSteelOffset;
  CLHEP::Hep3Vector _RightHallSteelOffset;
  CLHEP::Hep3Vector _DownstreamHallSteelOffset;
  CLHEP::Hep3Vector _UpstreamHallSteelOffset;

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
  double              _wallUnistrutHalfThickness;

  std::vector<int>    _shieldR_NumberOfModules;
  std::vector<int>    _shieldL_NumberOfModules;
  std::vector<int>    _shieldD_NumberOfModules;
  std::vector<int>    _shieldU_NumberOfModules;
  std::vector<int>    _shieldT_NumberOfModules;
  std::vector<int>    _shieldB_NumberOfModules;
  std::vector<int>    _shieldTS_NumberOfModules;

  std::vector<double> _shieldR_Offset;
  std::vector<double> _shieldL_Offset;
  std::vector<double> _shieldD_Offset;
  std::vector<double> _shieldU_Offset;
  std::vector<double> _shieldT_Offset;
  std::vector<double> _shieldB_Offset;
  std::vector<double> _shieldTS_Offset;

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

  double _scintillatorShieldOffsetToTheSideOfHallSteel;

};

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CosmicRayShieldMaker_hh */
