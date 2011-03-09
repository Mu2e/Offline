//
// Construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.cc,v 1.2 2011/03/09 19:46:51 genser Exp $
// $Author: genser $ 
// $Date: 2011/03/09 19:46:51 $
//
// Original author KLG based on Rob Kutschke's ...Maker classes
//

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "CosmicRayShieldGeom/inc/CosmicRayShieldMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BeamlineGeom/inc/Beamline.hh"


using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  CosmicRayShieldMaker::CosmicRayShieldMaker(SimpleConfig const & _config)
  {

    _diagLevel = 0;

    _crs = auto_ptr<CosmicRayShield>(new CosmicRayShield());

    if( ! _config.getBool("hasCosmicRayShield",false) ) return;

    parseConfig(_config);

    calculateCRSOffsets(_config);

    if ( _diagLevel > 0) {
      cout << __func__ << " _crs->_localOffset : " << _crs->_localOffset << endl;
      cout << __func__ << " _crs->_globalOffset: " << _crs->_globalOffset << endl;
    }

    makeCRSSteelShield(_config);

    makeShields();


  }

  void CosmicRayShieldMaker::parseConfig( SimpleConfig const & _config ){

    // we readin/store crs parameters needed in more than one function

    _HallSteelHalfThick     = _config.getDouble("fluxcrv.HallSteelHalfThick");
    _HallSteelHalfLengthXY  = _config.getDouble("fluxcrv.HallSteelHalfLengthXY");
    _HallSteelHalfLengthZ   = _config.getDouble("fluxcrv.HallSteelHalfLengthZ");
    _HallSteelMaterialName  = _config.getString("fluxcrv.HallSteelMaterialName");
    _HallSteelHoleRadius    = _config.getDouble("fluxcrv.HallSteelHoleRadius");
    _config.getVectorDouble("fluxcrv.HallSteelOffset",_HallSteelOffset,3);

    _scintillatorLayersPerModule  = _config.getInt("crs.scintillatorLayersPerModule");
    _scintillatorBarsPerFullLayer = _config.getInt("crs.scintillatorBarsPerFullLayer");

    if ( _scintillatorBarsPerFullLayer%2!=0) {
      throw cms::Exception("GEOM")
        << "crs.scintillatorBarsPerFullLayer number should be even\n";
    }

    _config.getVectorDouble("crs.scintillatorBarHalfLengths",_scintillatorBarHalfLengths,3);
    _scintillatorBarMaterialName  = _config.getString("crs.scintillatorBarMaterialName");
    _scintillatorLayerShift       = _config.getDouble("crs.scintillatorLayerShift");
    _scintillatorLayerGap       = _config.getDouble("crs.scintillatorLayerGap");

    _scintillatorBarPigmentationHalfThickness       =
      _config.getDouble("crs.scintillatorBarPigmentationHalfThickness");
    _scintillatorBarPigmentationMaterialName        = 
      _config.getString("crs.scintillatorBarPigmentationMaterialName");
    _config.getVectorDouble("crs.scintillatorModuleOuterSheetHalfLengths",
                            _scintillatorModuleOuterSheetHalfLengths,3);
    _scintillatorModuleOuterSheetMaterialName       =
      _config.getString("crs.scintillatorModuleOuterSheetMaterial");
    _scintillatorModuleInterLayerSheetMaterialName  =
      _config.getString("crs.scintillatorModuleInterLayerSheetMaterialName");
    _scintillatorModuleInterLayerSheetHalfThickness =
      _config.getDouble("crs.scintillatorModuleInterLayerSheetHalfThickness");
    _scintillatorOverlap = _config.getDouble("crs.scintillatorOverlap");

    _config.getVectorInt("crs.shieldR_NumberOfModules",_shieldR_NumberOfModules,2);
    _config.getVectorInt("crs.shieldL_NumberOfModules",_shieldL_NumberOfModules,2);
    _config.getVectorInt("crs.shieldD_NumberOfModules",_shieldD_NumberOfModules,2);
    _config.getVectorInt("crs.shieldU_NumberOfModules",_shieldU_NumberOfModules,2);
    _config.getVectorInt("crs.shieldT_NumberOfModules",_shieldT_NumberOfModules,2);
    _config.getVectorInt("crs.shieldB_NumberOfModules",_shieldB_NumberOfModules,2);
    _config.getVectorInt("crs.shieldTS_NumberOfModules",_shieldTS_NumberOfModules,2);

    _config.getVectorDouble("crs.shieldR_Offset",_shieldR_Offset,3);
    _config.getVectorDouble("crs.shieldL_Offset",_shieldL_Offset,3);
    _config.getVectorDouble("crs.shieldD_Offset",_shieldD_Offset,3);
    _config.getVectorDouble("crs.shieldU_Offset",_shieldU_Offset,3);
    _config.getVectorDouble("crs.shieldT_Offset",_shieldT_Offset,3);
    _config.getVectorDouble("crs.shieldB_Offset",_shieldB_Offset,3);
    _config.getVectorDouble("crs.shieldTS_Offset",_shieldTS_Offset,3);

    _config.getVectorDouble("crs.moduleUnistrutHalfLengths",_moduleUnistrutHalfLengths,3);
    _wallUnistrutHalfThickness  = _config.getDouble("crs.wallUnistrutHalfThickness");

  }

  void CosmicRayShieldMaker::makeDetails() {

    CRSScintillatorBarDetail& detail = _crs->_barDetails;
    
    // there is only one detail for now; we'll stick module materials there as well
    // this should be split once more details will be needed

    detail._id=0;

    detail._materialNames.push_back(_scintillatorBarMaterialName);
    detail._materialNames.push_back(_scintillatorBarPigmentationMaterialName);
    detail._materialNames.push_back(_scintillatorModuleOuterSheetMaterialName);
    detail._materialNames.push_back(_scintillatorModuleInterLayerSheetMaterialName);

    detail._halfLengths = _scintillatorBarHalfLengths;

  }

  void CosmicRayShieldMaker::makeShields() {

    // the defining factors for a shield is the number of Full/Half modules
    // the rest follows from that, the module size/overlaps define the size of one shield
    
    // like Straws and ScintillatorBars the ScintillatorModules have "universal" and
    // specific features, universal ones are their dimensions, specific are
    // their placements and the fact if they are full or half modules

    // a half module has different dimentions from the full one... 

    // the "universal" parameters are kept in  "details"

    // calculate  "global" parameters
    calculateCommonCRSScintillatorParameters();

    makeDetails();

    // shield counter
    int ishield=0;
    
    // we generate/assign the names here:

    // R Shield is the "Right" shield

    std::string name = "CRSScintillatorRShield";

    std::vector<int> numberOfModules = _shieldR_NumberOfModules;

    // the shields are placed in a "nominal" postion first

    CLHEP::Hep3Vector CRSScintillatorShieldOffset = 
      CLHEP::Hep3Vector(_scintillatorShieldOffsetToTheSideOfHallSteel + _HallSteelHalfLengthXY, 0., 0.) + 
      CLHEP::Hep3Vector(_shieldR_Offset[0],_shieldR_Offset[1],_shieldR_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " << 
        CRSScintillatorShieldOffset << endl;
    }

    // we rotate the shields so that they "start" from the downstream end

    std::vector<double> CRSScintillatorShieldRotationAngles; // x,y,z
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);

    // we may want to have an enumeration/function translating from 
    // R, L, T, B, D, U, TS to numbers, 
    // as of now the "translation" is done by using ishield:
    // 0  1  2  3  4  5  6 

    // the constructors are "simple", most work is done in the maker

    _crs->_scintillatorShields[name] = 
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldOffset + _crs->_localOffset,  // in Hall Air
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _crs->_globalOffset, // in Mu2e
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;
    
    //

    name = "CRSScintillatorLShield";

    numberOfModules = _shieldL_NumberOfModules;

    CRSScintillatorShieldOffset = 
      CLHEP::Hep3Vector(_scintillatorShieldOffsetToTheSideOfHallSteel + _HallSteelHalfLengthXY, 0., 0.) + 
      CLHEP::Hep3Vector(_shieldL_Offset[0],_shieldL_Offset[1],_shieldL_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " << 
        CRSScintillatorShieldOffset << endl;

    }

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(0.);

    _crs->_scintillatorShields[name] = 
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldOffset + _crs->_localOffset,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _crs->_globalOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorTShield";

    numberOfModules = _shieldT_NumberOfModules;

    CRSScintillatorShieldOffset = 
      CLHEP::Hep3Vector(_scintillatorShieldOffsetToTheSideOfHallSteel + _HallSteelHalfLengthXY, 0., 0.) + 
      CLHEP::Hep3Vector(_shieldT_Offset[0],_shieldT_Offset[1],_shieldT_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " << 
        CRSScintillatorShieldOffset << endl;
    }

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::halfpi);

    _crs->_scintillatorShields[name] = 
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldOffset + _crs->_localOffset,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _crs->_globalOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;


    //

    name = "CRSScintillatorBShield";

    numberOfModules = _shieldB_NumberOfModules;

    CRSScintillatorShieldOffset = 
      CLHEP::Hep3Vector(_scintillatorShieldOffsetToTheSideOfHallSteel + _HallSteelHalfLengthXY, 0., 0.) + 
      CLHEP::Hep3Vector(_shieldB_Offset[0],_shieldB_Offset[1],_shieldB_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " << 
        CRSScintillatorShieldOffset << endl;
    }

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);
    CRSScintillatorShieldRotationAngles.push_back(-CLHEP::halfpi);

    _crs->_scintillatorShields[name] = 
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldOffset + _crs->_localOffset,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _crs->_globalOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorDShield";

    numberOfModules = _shieldD_NumberOfModules;

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_scintillatorShieldOffsetToTheSideOfHallSteel + 
                        _HallSteelHalfLengthZ - _HallSteelHalfThick, 0., 0.) + 
      CLHEP::Hep3Vector(_shieldD_Offset[0],_shieldD_Offset[1],_shieldD_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " << 
        CRSScintillatorShieldOffset << endl;
    }

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(-CLHEP::halfpi);
    CRSScintillatorShieldRotationAngles.push_back(0.);

    _crs->_scintillatorShields[name] = 
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldOffset + _crs->_localOffset,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _crs->_globalOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;


    for (std::map<std::string,CRSScintillatorShield>::iterator itshield=_crs->_scintillatorShields.begin();
         itshield!=_crs->_scintillatorShields.end(); ++itshield) {

      if ( _diagLevel > 0) {
        cout << __func__ << " shield._name        : " << (itshield->second)._name << endl;
        cout << __func__ << " shield._localOffset : " << (itshield->second)._localOffset << endl;
        cout << __func__ << " shield._globalOffset: " << (itshield->second)._globalOffset << endl;
      }

      makeModules(itshield->second);

    }
    
  }

  void CosmicRayShieldMaker::makeModules(CRSScintillatorShield& shield) {

    // calculating total (half) width taken by all the modules

    int numberOfFullModules = shield._numberOfFullModules;
    int numberOfHalfModules = shield._numberOfHalfModules;
    int numberOfModules     = numberOfFullModules + numberOfHalfModules;

    double shieldHalfWidth = (2*numberOfFullModules*_scintillatorFullModuleHalfWidth +
                              2*numberOfHalfModules*_scintillatorHalfModuleHalfWidth - 
                              (numberOfModules-1)*_scintillatorModuleOverlap)*.5;

    if ( _diagLevel > 0) {
      cout << __func__ << " Shield              : " << shield._name << endl;
      cout << __func__ << " numberOfModules     : " << numberOfModules << endl;
      cout << __func__ << " numberOfFullModules : " << numberOfFullModules << endl;
      cout << __func__ << " numberOfHalfModules : " << numberOfHalfModules << endl;
      cout << __func__ << " shieldHalfWidth     : " << shieldHalfWidth << endl;
      cout << __func__ << " _scintillatorFullModuleHalfWidth : " << _scintillatorFullModuleHalfWidth << endl;
      cout << __func__ << " _scintillatorModuleOverlap       : " << _scintillatorModuleOverlap << endl;
    }

    // module counter for the given shield
    int imodule = 0;
    // we place/make the full modules first

    double firstFullModuleLocalOffsetL = -shieldHalfWidth + _scintillatorFullModuleHalfWidth;

    for (int ii = 0; ii<numberOfFullModules; ++ii) {

      if ( _diagLevel > 0) {
        cout << __func__ << " creating module          : " << imodule << endl;
      }

      // module local offsets distance from the center of the shield for a given module

      double moduleLocalOffsetL = firstFullModuleLocalOffsetL +
        ii*(2.*_scintillatorFullModuleHalfWidth-_scintillatorModuleOverlap);

      if ( _diagLevel > 0) {
        cout << __func__ << " moduleLocalOffsetL               : " << moduleLocalOffsetL << endl;
      }

      // creating an empty module (need to know if full/half)
      // than modifying it

      shield._modules.push_back(CRSScintillatorModule( CRSScintillatorModuleId(shield._id,imodule),
                                                       _scintillatorBarsPerFullLayer));

      CRSScintillatorModule& module = shield._modules.back();

      module._layers.reserve(_scintillatorLayersPerModule);

      // we need to shift the modules transversly as well
      // we shall assume that the first (0th) one will be touching the steel (%2 below)

      module._localOffset = CLHEP::Hep3Vector( imodule%2 != 0 ? 
                                               _scintillatorModuleCoreHalfThickness : 
                                               -_scintillatorModuleCoreHalfThickness,
                                               0., moduleLocalOffsetL);

      module._globalRotationAngles = shield._globalRotationAngles;

      ++imodule;

    }

    // same for half modules (may create a function for the two)

    double firstHalfModuleLocalOffsetL = -shieldHalfWidth + 
      numberOfFullModules*(2.*_scintillatorFullModuleHalfWidth - _scintillatorModuleOverlap)
      + _scintillatorHalfModuleHalfWidth;

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorModuleOverlap       : " << _scintillatorModuleOverlap << endl;
    }

    for (int ii = 0; ii<numberOfHalfModules; ++ii) {

      if ( _diagLevel > 0) {
        cout << __func__ << " creating module          : " << imodule << endl;
      }

      double moduleLocalOffsetL = firstHalfModuleLocalOffsetL +
        ii*(2.0*_scintillatorHalfModuleHalfWidth-_scintillatorModuleOverlap);

      if ( _diagLevel > 0) {
        cout << __func__ << " moduleLocalOffsetL               : " << moduleLocalOffsetL << endl;
      }

      // creating an empty module than modifying it

      shield._modules.push_back(CRSScintillatorModule( CRSScintillatorModuleId(shield._id,imodule),
                                                       _scintillatorBarsPerFullLayer/2 ));

      CRSScintillatorModule& module = shield._modules.back();

      module._layers.reserve(_scintillatorLayersPerModule);

      module._localOffset = CLHEP::Hep3Vector( imodule%2 != 0 ? 
                                               _scintillatorModuleCoreHalfThickness : 
                                               -_scintillatorModuleCoreHalfThickness,
                                               0., moduleLocalOffsetL);

      module._globalRotationAngles = shield._globalRotationAngles;

      ++imodule;

    }

    // set global offsets, make layers

    for (int ii = 0; ii<numberOfModules; ++ii) {

      if ( _diagLevel > 0) {
        cout << __func__ << " making layers for module : " << ii << endl;
      }

      CRSScintillatorModule& module = shield._modules.at(ii);
      module._globalOffset = shield._globalOffset + module._localOffset;
      if ( _diagLevel > 0) {
        cout << __func__ << " module._localOffset : " << module._localOffset << endl;
        cout << __func__ << " module._globalOffset: " << module._globalOffset << endl;
      }
      makeLayers(module);

    }

  }
  
  void CosmicRayShieldMaker::makeLayers(CRSScintillatorModule& module) {

    // each module has the same number of layers, but not the same number of bars per layer

    // each module knows its number of bars per layer

    int numberOfBarsPerLayer = module._nBarsPerLayer;
    int numberOfLayers       = _scintillatorLayersPerModule;

    // calculating longitudinal offset

    double firtstLayerLocalOffsetL = _scintillatorLayerShift*0.5*(1-_scintillatorLayersPerModule);

    // calculating transverse offset

    double scintillatorLayerShiftT = 2.0*(_scintillatorBarHalfLengths[0]+
                                          _scintillatorModuleInterLayerSheetHalfThickness);

    double firtstLayerLocalOffsetT = 
      (_scintillatorBarHalfLengths[0]+
       _scintillatorModuleInterLayerSheetHalfThickness)*
      (1-_scintillatorLayersPerModule);


    for (int ii = 0; ii<numberOfLayers; ++ii) {

      if ( _diagLevel > 0) {
        cout << __func__ << " making layer : " << ii << endl;
      }

      module._layers.push_back(CRSScintillatorLayer( CRSScintillatorLayerId(module._id,ii),
                                                     numberOfBarsPerLayer));
      CRSScintillatorLayer& layer = module._layers.back();

      // fill in the layer offsets etc...

      // we will place innermost layers shifted most to the left
      layer._localOffset = CLHEP::Hep3Vector(firtstLayerLocalOffsetT + ii*scintillatorLayerShiftT,
                                             0., 
                                             firtstLayerLocalOffsetL + ii*_scintillatorLayerShift);
      layer._globalRotationAngles = module._globalRotationAngles;
      layer._globalOffset = module._globalOffset + layer._localOffset;
      layer._nBars = module._nBarsPerLayer;

      if ( _diagLevel > 0) {
        cout << __func__ << " layer._localOffset : " << layer._localOffset << endl;
        cout << __func__ << " layer._globalOffset: " << layer._globalOffset << endl;
      }

      makeBars(layer);

    }

  }

  void CosmicRayShieldMaker::makeBars(CRSScintillatorLayer& layer) {

    // put the bar in its global "bar registry/place" 
    int numberOfBars = layer._nBars;

    layer._bars.reserve(numberOfBars);

    double barSpaceHalfWidth = (numberOfBars == _scintillatorBarsPerFullLayer) ?
      _scintillatorFullLayerHalfWidth : _scintillatorHalfLayerHalfWidth;

    double firtstBarLocalOffsetL = -barSpaceHalfWidth + _scintillatorBarHalfLengths[2];

    for (int ii = 0; ii<numberOfBars; ++ii) {

      if ( _diagLevel > 1) {
        cout << __func__ << " making bar   : " << ii << endl;
      }

      // calculate the bar offsets, ids, indeces and enter it into the container

      CRSScintillatorBarIndex index(_crs->_allCRSScintillatorBars.size());

      // local longitudinal offset wrt the center of the layer

      _crs->_allCRSScintillatorBars.push_back(CRSScintillatorBar(CRSScintillatorBarId(layer._id,ii),
                                                                 index));

      CRSScintillatorBar& bar = _crs->_allCRSScintillatorBars.back();

      bar._localOffset = CLHEP::Hep3Vector(0.,0.,
                                           firtstBarLocalOffsetL + ii*(2.*_scintillatorBarHalfLengths[2]+
                                                                       _scintillatorLayerGap));

      bar._globalRotationAngles = layer._globalRotationAngles;

      // untill now the object positions were calculated locally
      // wrt to the the module etc... in a "nominal" position

      // we will now calculate position of the bar wrt the CRS origin and rotate it wrt that origin

      CLHEP::Hep3Vector barCRSOffset = layer._globalOffset + bar._localOffset - _crs->_globalOffset;

      CLHEP::HepRotationX RX(bar._globalRotationAngles[0]);
      CLHEP::HepRotationY RY(bar._globalRotationAngles[1]);
      CLHEP::HepRotationZ RZ(bar._globalRotationAngles[2]);

      CLHEP::HepRotation barRotation(RX*RY*RZ);

      CLHEP::Hep3Vector barCRSOffsetRotated = barRotation * barCRSOffset;

      //  bar._globalOffset = layer._globalOffset + bar._localOffset;

      //  note that difference, see the comments above, rotation done for the bars only so far

      bar._globalOffset = barCRSOffsetRotated + _crs->_globalOffset;

      layer._bars.push_back(&bar);
      bar._detail = &_crs->_barDetails;
      layer._indices.push_back(index);

      if ( _diagLevel > 1) {
        cout << __func__ << " barCRSOffset        : " << barCRSOffset        << endl;
        cout << __func__ << " barCRSOffsetRotated : " << barCRSOffsetRotated << endl;
        cout << __func__ << " bar._localOffset    : " << bar._localOffset    << endl;
        cout << __func__ << " bar._globalOffset   : " << bar._globalOffset   << endl;
      }

    }

  }

  void CosmicRayShieldMaker::calculateCommonCRSScintillatorParameters() {

    // the modules have struts, outerSheet, scintillator layers, interLayerSheets

    // this is the combined width of the scintillator bars only
    _scintillatorFullLayerHalfWidth = 
      _scintillatorBarsPerFullLayer*_scintillatorBarHalfLengths[2] + 
      (_scintillatorBarsPerFullLayer-1)*_scintillatorLayerGap*0.5;

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorFullLayerHalfWidth  : " <<
        _scintillatorFullLayerHalfWidth << endl;
    }

//     _scintillatorHalfLayerHalfWidth = _scintillatorFullLayerHalfWidth -
//       (_scintillatorBarsPerFullLayer/2)*(_scintillatorBarHalfLengths[2]+_scintillatorLayerGap);

    _scintillatorHalfLayerHalfWidth = 
      _scintillatorBarsPerFullLayer/2*_scintillatorBarHalfLengths[2] + 
      (_scintillatorBarsPerFullLayer/2-1)*_scintillatorLayerGap*0.5;

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorHalfLayerHalfWidth   : " <<
        _scintillatorHalfLayerHalfWidth << endl;
    }

    // module is longer than the layers due to the side bracket
    _scintillatorFullModuleHalfWidth = _scintillatorModuleOuterSheetHalfLengths[2]; 
    _scintillatorHalfModuleHalfWidth = 
      _scintillatorFullModuleHalfWidth - _scintillatorHalfLayerHalfWidth;

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorFullModuleHalfWidth : " <<
        _scintillatorFullModuleHalfWidth << endl;
    }

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorHalfModuleHalfWidth : " <<
        _scintillatorHalfModuleHalfWidth << endl;
    }

    if (_scintillatorFullModuleHalfWidth<
        (_scintillatorFullLayerHalfWidth + (_scintillatorLayersPerModule-1)/2*_scintillatorLayerShift)) {
      throw cms::Exception("GEOM")
        << "inconsistent data crs.scintillatorModuleOuterSheetHalfLengths to small?\n";
    }

    _scintillatorLayerHalfLength  = _scintillatorBarHalfLengths[1];
    _scintillatorModuleHalfLength = _scintillatorBarHalfLengths[1];

    // given the required _scintillatorOverlap
    _scintillatorModuleOverlap = _scintillatorOverlap + 
      2.*(_scintillatorFullModuleHalfWidth - 
          _scintillatorFullLayerHalfWidth - 
          (_scintillatorLayersPerModule-1)/2*_scintillatorLayerShift);

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorModuleOverlap       : " <<
        _scintillatorModuleOverlap << endl;
    }

    // outer thickness should be the same for all the shields/modules/layers
    _scintillatorModuleHalfThickness = 
      2.*(_moduleUnistrutHalfLengths[0]+
          _scintillatorModuleOuterSheetHalfLengths[0]) +
      (_scintillatorLayersPerModule-1)*_scintillatorModuleInterLayerSheetHalfThickness + 
      _scintillatorLayersPerModule*_scintillatorBarHalfLengths[0];

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorModuleHalfThickness : " << 
        _scintillatorModuleHalfThickness << endl;
    }

    _scintillatorModuleCoreHalfThickness = 
      2.*(_scintillatorModuleOuterSheetHalfLengths[0]) +
      (_scintillatorLayersPerModule-1)*_scintillatorModuleInterLayerSheetHalfThickness + 
      _scintillatorLayersPerModule*_scintillatorBarHalfLengths[0];

    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorModuleCoreHalfThickness : " << 
        _scintillatorModuleCoreHalfThickness << endl;
    }

    // same calc
    double scintillatorShieldHalfThickness1 = _scintillatorModuleHalfThickness + 
      2.*_scintillatorModuleOuterSheetHalfLengths[0] +
      (_scintillatorLayersPerModule-1)*_scintillatorModuleInterLayerSheetHalfThickness + 
      _scintillatorLayersPerModule*_scintillatorBarHalfLengths[0];
    
    if ( _diagLevel > 0) {
      cout << __func__ << " scintillatorShieldHalfThickness1 : " << 
        scintillatorShieldHalfThickness1 << endl;
    }

    _scintillatorShieldHalfThickness = 
      2.*(_scintillatorModuleHalfThickness - _moduleUnistrutHalfLengths[0]);
    
    if ( _diagLevel > 0) {
      cout << __func__ << " _scintillatorShieldHalfThickness : " << 
        _scintillatorShieldHalfThickness << endl;
    }

    // this is offset outward from the plane centered in the center of the steel wall
    //
    _scintillatorShieldOffsetToTheSideOfHallSteel = _HallSteelHalfThick + 
      2.*_wallUnistrutHalfThickness + _scintillatorShieldHalfThickness;

    if ( _diagLevel > 0) {
      cout << __func__ << " ShieldOutwardOffset  : " << 
        _scintillatorShieldOffsetToTheSideOfHallSteel << endl;
    }

    // the shield has modules overlaping on their edges, hence a global 2.* below
    // see fig. 3.4 in meco-crs-05-001

    // should be the same as above
    double ShieldOutwardOffset1 = 
      _HallSteelHalfThick +
      2.*(_wallUnistrutHalfThickness    +
          _moduleUnistrutHalfLengths[0] +
          2.*_scintillatorModuleOuterSheetHalfLengths[0] +
          _scintillatorLayersPerModule*_scintillatorBarHalfLengths[0]+
          (_scintillatorLayersPerModule-1)*_scintillatorModuleInterLayerSheetHalfThickness);

    if ( _diagLevel > 0) {
      cout << __func__ << " ShieldOutwardOffset1 : " << ShieldOutwardOffset1 << endl;
    }

    // we need to allocate space for all the bars at once to avoid pointer invalidatation

    _totalNumberOfBars = 
      (_shieldR_NumberOfModules[0] +
       _shieldL_NumberOfModules[0] +
       _shieldD_NumberOfModules[0] +
       _shieldU_NumberOfModules[0] +
       _shieldT_NumberOfModules[0] +
       _shieldB_NumberOfModules[0] +
       _shieldTS_NumberOfModules[0])*_scintillatorLayersPerModule*_scintillatorBarsPerFullLayer +     
      (_shieldR_NumberOfModules[1] + 
       _shieldL_NumberOfModules[1] + 
       _shieldD_NumberOfModules[1] + 
       _shieldU_NumberOfModules[1] + 
       _shieldT_NumberOfModules[1] + 
       _shieldB_NumberOfModules[1] + 
       _shieldTS_NumberOfModules[1])*_scintillatorLayersPerModule*(_scintillatorBarsPerFullLayer/2);

    if ( _diagLevel > 0) {
      cout << __func__ << " _totalNumberOfBars : " << _totalNumberOfBars << endl;
    }
    _crs->_allCRSScintillatorBars.reserve(_totalNumberOfBars);

  }

  void CosmicRayShieldMaker::makeCRSSteelShield(SimpleConfig const & _config) {
    // first make the steel (fluxreturn)

    // Compute dimensions of 5 sides in Mu2e coordinates
    double CRSSteelShieldTopHalfX   = _HallSteelHalfLengthXY + _HallSteelHalfThick;
    double CRSSteelShieldTopHalfY   = _HallSteelHalfThick;
    double CRSSteelShieldTopHalfZ   = _HallSteelHalfLengthZ;
    double CRSSteelShieldSideHalfX  = _HallSteelHalfThick;
    double CRSSteelShieldSideHalfY  = _HallSteelHalfLengthXY - _HallSteelHalfThick;
    double CRSSteelShieldSideHalfZ  = _HallSteelHalfLengthZ; 
    double CRSSteelShieldUpstreamHalfX = _HallSteelHalfLengthXY - _HallSteelHalfThick;
    double CRSSteelShieldUpstreamHalfY = _HallSteelHalfLengthXY - _HallSteelHalfThick;
    double CRSSteelShieldUpstreamHalfZ = _HallSteelHalfThick;

    double CRSSteelShieldTopDims[3] ={
      CRSSteelShieldTopHalfX,
      CRSSteelShieldTopHalfY,
      CRSSteelShieldTopHalfZ
    };

    double CRSSteelShieldSideDims[3] ={
      CRSSteelShieldSideHalfX,
      CRSSteelShieldSideHalfY,
      CRSSteelShieldSideHalfZ
    };

    double CRSSteelShieldUpstreamDims[3] ={
      CRSSteelShieldUpstreamHalfX,
      CRSSteelShieldUpstreamHalfY,           
      CRSSteelShieldUpstreamHalfZ                        
    };

    if ( _diagLevel > 0) {

      CLHEP::Hep3Vector CRSSteelShieldTopDimsV(CRSSteelShieldTopHalfX,
                                               CRSSteelShieldTopHalfY,
                                               CRSSteelShieldTopHalfZ);

      CLHEP::Hep3Vector CRSSteelShieldSideDimsV(CRSSteelShieldSideHalfX,
                                                CRSSteelShieldSideHalfY,
                                                CRSSteelShieldSideHalfZ);

      CLHEP::Hep3Vector CRSSteelShieldUpstreamDimsV(CRSSteelShieldUpstreamHalfX,
                                                 CRSSteelShieldUpstreamHalfY,           
                                                 CRSSteelShieldUpstreamHalfZ);

      cout << __func__ << " CRSSteelShieldTopDimsV   : " << CRSSteelShieldTopDimsV   << endl;
      cout << __func__ << " CRSSteelShieldSideDimsV  : " << CRSSteelShieldSideDimsV  << endl;
      cout << __func__ << " CRSSteelShieldUpstreamDimsV : " << CRSSteelShieldUpstreamDimsV << endl;

    }

    _TopHallSteelOffset   = CLHEP::Hep3Vector(0.,      _HallSteelHalfLengthXY,  0.);
    _BottomHallSteelOffset= CLHEP::Hep3Vector(0.,    -(_HallSteelHalfLengthXY), 0.);
    _LeftHallSteelOffset  = CLHEP::Hep3Vector(         _HallSteelHalfLengthXY,  0., 0.);
    _RightHallSteelOffset = CLHEP::Hep3Vector(       -(_HallSteelHalfLengthXY), 0., 0.);
    _DownstreamHallSteelOffset  = CLHEP::Hep3Vector(0., 0.,  _HallSteelHalfLengthZ - _HallSteelHalfThick);
    _UpstreamHallSteelOffset = CLHEP::Hep3Vector(0., 0.,-(_HallSteelHalfLengthZ - _HallSteelHalfThick));

    if ( _diagLevel > 0) {

      cout << __func__ << " _TopHallSteelOffset     : " << _TopHallSteelOffset    << endl;
      cout << __func__ << " _BottomHallSteelOffset  : " << _BottomHallSteelOffset << endl;
      cout << __func__ << " _LeftHallSteelOffset    : " << _LeftHallSteelOffset   << endl;
      cout << __func__ << " _RightHallSteelOffset   : " << _RightHallSteelOffset  << endl;
      cout << __func__ << " _DownstreamHallSteelOffset    : " << _DownstreamHallSteelOffset   << endl;
      cout << __func__ << " _UpstreamHallSteelOffset   : " << _UpstreamHallSteelOffset  << endl;

    }

    // finaly create the steel shield objects

    std::string name = "CRSSteelTopShield";
    _crs->_steelShields[name] = 
      CRSSteelShield(name,
                     _TopHallSteelOffset + _crs->_localOffset,  // offset in the hall Air
                     0,                                           // HepRotation*
                     _TopHallSteelOffset + _crs->_globalOffset, // global offset in Mu2e
                     CRSSteelShieldTopDims);

    name = "CRSSteelBottomShield";
    _crs->_steelShields[name] = 
      CRSSteelShield(name,
                     _BottomHallSteelOffset + _crs->_localOffset,
                     0,
                     _BottomHallSteelOffset + _crs->_globalOffset,
                     CRSSteelShieldTopDims);
                     
    name = "CRSSteelLeftShield";
    _crs->_steelShields[name] = 
      CRSSteelShield(name,
                     _LeftHallSteelOffset + _crs->_localOffset,
                     0,
                     _LeftHallSteelOffset + _crs->_globalOffset,
                     CRSSteelShieldSideDims);

    name = "CRSSteelRightShield";
    _crs->_steelShields[name] = 
      CRSSteelShield(name,
                     _RightHallSteelOffset + _crs->_localOffset,
                     0,
                     _RightHallSteelOffset + _crs->_globalOffset,
                     CRSSteelShieldSideDims);

    name = "CRSSteelDownstreamShield";
    _crs->_steelShields[name] = 
      CRSSteelShield(name,
                     _DownstreamHallSteelOffset + _crs->_localOffset,
                     0,
                     _DownstreamHallSteelOffset + _crs->_globalOffset,
                     CRSSteelShieldUpstreamDims);

    name = "CRSSteelUpstreamShield";
    _crs->_steelShields[name] = 
      CRSSteelShield(name,
                     _UpstreamHallSteelOffset + _crs->_localOffset,
                     0,
                     _UpstreamHallSteelOffset + _crs->_globalOffset,
                     CRSSteelShieldUpstreamDims,
                     _config.getDouble("fluxcrv.HallSteelHoleRadius")
                     );
    
  }


  void CosmicRayShieldMaker::calculateCRSOffsets(SimpleConfig const & _config) {

    double dsCoilZ0          = _config.getDouble("toyDS.z0");

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();

    CLHEP::Hep3Vector detSolCoilPosition(solenoidOffset, 0., -dsCoilZ0);

    if ( _diagLevel > 0) {
      cout << __func__ << " detSolCoilPosition  : " << detSolCoilPosition  << endl;
    }

    // reconstructing the position of the hall air in the world etc...

    vector<double> worldHLen;
    _config.getVectorDouble("world.halfLengths", worldHLen, 3);

    // Get parameters related to the overall dimensions of the hall and to
    // the earthen overburden.

    double floorThick           =  _config.getDouble("hall.floorThick");
    double ceilingThick         =  _config.getDouble("hall.ceilingThick");
    double overburdenDepth      =  _config.getDouble("dirt.overburdenDepth");
    vector<double> hallInHLen;
    _config.getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Top of the floor in G4 world coordinates.
    double yFloor   = -worldHLen[1] + floorThick;

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*hallInHLen[1];
    
    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + ceilingThick;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + overburdenDepth;
    
    // Half length and y origin of the dirt box.
    double yLDirt = ( ySurface + worldHLen[1] )/2.;
    double y0Dirt = -worldHLen[1] + yLDirt;
    
    // Center of the dirt box, in the G4 world system.
    CLHEP::Hep3Vector dirtOffset(0.,y0Dirt,0.);

    if ( _diagLevel > 0) {
      cout << __func__ << " dirtOffset          : " << dirtOffset    << endl;
    }

    // Position of the center of the hall in the world volume.
    vector<double> hallPosition;
    _config.getVectorDouble("hall.offset", hallPosition,3);

    if ( _diagLevel > 0) {
      cout << __func__ << " hallPosition/Offset : " << CLHEP::Hep3Vector(hallPosition[0],
                                                                         hallPosition[1],
                                                                         hallPosition[2])
           << endl;
    }

    double hallY0 = yFloor + hallInHLen[1] + hallPosition[1];

    // Center of the concrete volume in the coordinate system of the dirt.
    // It is a "local" offset

    CLHEP::Hep3Vector wallOffset = 
      CLHEP::Hep3Vector(hallPosition[0], hallY0, hallPosition[2]) - dirtOffset; //

    if ( _diagLevel > 0) {
      cout << __func__ << " wallOffset          : " << wallOffset    << endl;
    }

    if ( _diagLevel > 0) {
      cout << __func__ << " wallOffsetInWorld   : " << wallOffset + dirtOffset  << endl;
    }

    // Origin of the hall air volume in the system of the hall concrete volume.
    CLHEP::Hep3Vector hallAirOffset = CLHEP::Hep3Vector( 0., (floorThick-ceilingThick)/2., 0.);

    if ( _diagLevel > 0) {
      cout << __func__ << " hallAirOffset          : " << hallAirOffset << endl;
    }

    if ( _diagLevel > 0) {
      cout << __func__ << " hallAirOffsetInWorld   : " << hallAirOffset + wallOffset + dirtOffset << endl;
    }

    // Position of the origin of the mu2e coordinate system
    CLHEP::Hep3Vector mu2eOrigin = CLHEP::Hep3Vector(
                                                     _config.getDouble("world.mu2eOrigin.xoffset"),
                                                     _config.getDouble("world.mu2eOrigin.height") + yFloor,
                                                     _config.getDouble("world.mu2eOrigin.zoffset")
                                                     );
    if ( _diagLevel > 0) {
      cout << __func__ << " mu2eOrigin          : " << mu2eOrigin    << endl;
    }

    std::vector<double> CRSSteelShieldOffsetSTDV;
    _config.getVectorDouble("fluxcrv.HallSteelOffset", CRSSteelShieldOffsetSTDV, 3);

    CLHEP::Hep3Vector CRSSteelShieldOffset(CRSSteelShieldOffsetSTDV[0],
                                           CRSSteelShieldOffsetSTDV[1],
                                           CRSSteelShieldOffsetSTDV[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSSteelShieldOffset     : " << CRSSteelShieldOffset  << endl;
    }

    // Imagine a box that exactly contains the flux return steel.
    // This is the center of that box, in the coordinate system of the mother volume (the hall air).

    _crs->_localOffset =  - (hallAirOffset + wallOffset + dirtOffset - mu2eOrigin + detSolCoilPosition + CRSSteelShieldOffset);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSOffset (in hallAir) : " << _crs->_localOffset    << endl;
    }

    _crs->_globalOffset = _crs->_localOffset + hallAirOffset + wallOffset + dirtOffset - mu2eOrigin;
    // _crs->_globalOffset = _crs->_localOffset + hallAirOffset + wallOffset + dirtOffset;

    // the global offset will be wrt mu2e

    if ( _diagLevel > 0) {
      cout << __func__ << " _crs->_globalOffset    : " << _crs->_globalOffset << endl;
    }

  }

} // namespace mu2e
