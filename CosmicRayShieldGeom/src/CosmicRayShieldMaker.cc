//
// Construct and return CosmicRayShield
//
// $Id: CosmicRayShieldMaker.cc,v 1.26 2013/08/27 17:10:28 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/08/27 17:10:28 $
//
// Original author KLG based on Rob Kutschke's ...Maker classes
//
// Notes
//
// Right now the hole in CRSSteelDownstreamShield is large enough to
// accommodate the DS3Vacuum; it should eventually be only large
// enough to accommodate MBS

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
 
// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "CosmicRayShieldGeom/inc/CosmicRayShieldMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "MBSGeom/inc/MBS.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  CosmicRayShieldMaker::CosmicRayShieldMaker(SimpleConfig const & _config, double solenoidOffset)
  {

    _crs = unique_ptr<CosmicRayShield>(new CosmicRayShield());

    if( ! _config.getBool("hasCosmicRayShield",false) ) return;

    parseConfig(_config);

    if ( _diagLevel > 0) {
      cout << __func__ << " _HallSteelOffset:    "       << _HallSteelOffset    << endl;
      cout << __func__ << " _scintillatorShieldOffset: " << _scintillatorShieldOffset << endl;
    }

    if ( _hasPassiveShield ) {
      makeCRSSteelShield(_config);
      _crs->_hasPassiveShield = true;
    }
    // this order of calculations is important
    if ( _hasActiveShield ) {
      makeShields();
      _crs->_hasActiveShield = true;
    }
  }

  void CosmicRayShieldMaker::parseConfig( SimpleConfig const & _config ){

    // we readin/store crs parameters needed in more than one function
    _diagLevel = _config.getInt("crs.verbosityLevel",0);

    _hasPassiveShield =  _config.getBool("crs.hasPassiveShield");
    _hasActiveShield  =  _config.getBool("crs.hasActiveShield");

    if(_hasPassiveShield == true)
    {
      std::cout<<"The passive (\"steel\") shield will be replaced by a new class."<<std::endl;
      std::cout<<"The code will remain in CosmicRayShieldMaker.cc until the new class is done."<<std::endl;
      std::cout<<"However, it will be disabled from now on."<<std::endl;
    }
    _hasPassiveShield = false; 

    if ( !_hasPassiveShield && !_hasActiveShield ){
      throw cet::exception("GEOM")
        << " CosmicRayShield requested but none of the parts selected \n";
    }

    if (_hasPassiveShield) {
      _HallSteelHalfThick           = _config.getDouble("fluxcrv.HallSteelHalfThick");
      _HallSteelHalfSideShieldHeight= _config.getDouble("fluxcrv.HallSteelHalfSideShieldHeight");
      _HallSteelHalfRShieldLength   = _config.getDouble("fluxcrv.HallSteelHalfRShieldLength");
      _HallSteelHalfLShieldLength   = _config.getDouble("fluxcrv.HallSteelHalfLShieldLength");
      _HallSteelHalfTShieldLength   = _config.getDouble("fluxcrv.HallSteelHalfTShieldLength");
      _HallSteelHalfTSRShieldLength = _config.getDouble("fluxcrv.HallSteelHalfTSRShieldLength");
      _HallSteelHalfTSLShieldLength = _config.getDouble("fluxcrv.HallSteelHalfTSLShieldLength");
      _HallSteelHalfTSTShieldLength = _config.getDouble("fluxcrv.HallSteelHalfTSTShieldLength");
      _HallSteelMaterialName    = _config.getString("fluxcrv.HallSteelMaterialName");
      _HallSteelHoleRadius      = _config.getDouble("fluxcrv.HallSteelHoleRadius");
      _HallSteelOffset          = _config.getHep3Vector("fluxcrv.HallSteelOffset");
      _HallSteelRShieldCenter   = _config.getHep3Vector("fluxcrv.HallSteelRShieldCenter");
      _HallSteelLShieldCenter   = _config.getHep3Vector("fluxcrv.HallSteelLShieldCenter");
      _HallSteelTShieldCenter   = _config.getHep3Vector("fluxcrv.HallSteelTShieldCenter");
      _HallSteelDShieldCenter   = _config.getHep3Vector("fluxcrv.HallSteelDShieldCenter");
      _HallSteelTSRShieldCenter = _config.getHep3Vector("fluxcrv.HallSteelTSRShieldCenter");
      _HallSteelTSLShieldCenter = _config.getHep3Vector("fluxcrv.HallSteelTSLShieldCenter");
      _HallSteelTSTShieldCenter = _config.getHep3Vector("fluxcrv.HallSteelTSTShieldCenter");
    }

    _scintillatorLayersPerModule  = _config.getInt("crs.scintillatorLayersPerModule");
    _scintillatorBarsPerFullLayer = _config.getInt("crs.scintillatorBarsPerFullLayer");

    if ( _scintillatorBarsPerFullLayer%2!=0) {
      throw cet::exception("GEOM")
        << "crs.scintillatorBarsPerFullLayer number should be even\n";
    }

    _config.getVectorDouble("crs.scintillatorBarHalfLengths",_scintillatorBarHalfLengths,3);
    _scintillatorBarMaterialName  = _config.getString("crs.scintillatorBarMaterialName");
    _scintillatorLayerShift       = _config.getDouble("crs.scintillatorLayerShift");
    _scintillatorLayerGap         = _config.getDouble("crs.scintillatorLayerGap");

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
    _config.getVectorInt("crs.shieldT_NumberOfModules",_shieldT_NumberOfModules,2);
    _config.getVectorInt("crs.shieldTSR_NumberOfModules",_shieldTSR_NumberOfModules,2);
    _config.getVectorInt("crs.shieldTSL_NumberOfModules",_shieldTSL_NumberOfModules,2);
    _config.getVectorInt("crs.shieldTST_NumberOfModules",_shieldTST_NumberOfModules,2);

    _scintillatorShieldOffset = _config.getHep3Vector("crs.scintillatorShieldOffset");

    _config.getVectorDouble("crs.shieldR_Offset",_shieldR_Offset,3);
    _config.getVectorDouble("crs.shieldL_Offset",_shieldL_Offset,3);
    _config.getVectorDouble("crs.shieldD_Offset",_shieldD_Offset,3);
    _config.getVectorDouble("crs.shieldT_Offset",_shieldT_Offset,3);
    _config.getVectorDouble("crs.shieldTSR_Offset",_shieldTSR_Offset,3);
    _config.getVectorDouble("crs.shieldTSL_Offset",_shieldTSL_Offset,3);
    _config.getVectorDouble("crs.shieldTST_Offset",_shieldTST_Offset,3);

    _config.getVectorDouble("crs.moduleUnistrutHalfLengths",_moduleUnistrutHalfLengths,3);
    //    _wallUnistrutHalfThickness  = _config.getDouble("crs.wallUnistrutHalfThickness");

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

    // the "universal", mainly bar, parameters are kept in  "details"

    // calculate "global" parameters
    calculateCommonCRSScintillatorParameters();

    makeDetails();

    // shield counter
    int ishield=0;

    // we generate/assign the names here:
    // R Shield is the "Right" shield

    // we may want to have an enumeration/function translating from
    // R, L, T, D, TSR, TSL, TST to numbers,
    // as of now the "translation" is done by using ishield:
    // 0  1  2  3  4  5  6

    std::string name = "CRSScintillatorRShield";

    std::vector<int> numberOfModules = _shieldR_NumberOfModules;

    // the shields are placed in a "nominal" postion first

    CLHEP::Hep3Vector CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldR_Offset[0],_shieldR_Offset[1],_shieldR_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;
    }

    // we set the rotation angles of the shields so that they "start" from the downstream end
    // but we do not actually rotate the shields, modules or layers as of now, we do it for bars only

    if (_hasPassiveShield) {
        
      CLHEP::Hep3Vector diff = CRSScintillatorShieldOffset + _scintillatorShieldOffset - 
        _HallSteelRShieldCenter - _HallSteelOffset;

      if ( fabs(diff.x())<0. ){
        throw cet::exception("GEOM")
          << "The CRPassiveShield && hasCRActiveShield configurations are inconsistent "
          << name
          << "\n";
      }
    }

    std::vector<double> CRSScintillatorShieldRotationAngles; // x,y,z
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);

    // the constructors are "simple", most work is done in the maker

    _crs->_scintillatorShields[name] =
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset, // in Mu2e
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorLShield";

    numberOfModules = _shieldL_NumberOfModules;

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldL_Offset[0],_shieldL_Offset[1],_shieldL_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;

    }

    if (_hasPassiveShield) {
        
      CLHEP::Hep3Vector diff = CRSScintillatorShieldOffset + _scintillatorShieldOffset - 
        _HallSteelRShieldCenter - _HallSteelOffset;

      if ( fabs(diff.x())<0. ){
        throw cet::exception("GEOM")
          << "The CRPassiveShield && hasCRActiveShield configurations are inconsistent "
          << name
          << "\n";
      }
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
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorTShield";

    numberOfModules = _shieldT_NumberOfModules;

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldT_Offset[0],_shieldT_Offset[1],_shieldT_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;
    }

    if (_hasPassiveShield) {
        
      CLHEP::Hep3Vector diff = CRSScintillatorShieldOffset + _scintillatorShieldOffset - 
        _HallSteelRShieldCenter - _HallSteelOffset;

      if ( fabs(diff.x())<0. ){
        throw cet::exception("GEOM")
          << "The CRPassiveShield && hasCRActiveShield configurations are inconsistent "
          << name
          << "\n";
      }
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
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorDShield";

    numberOfModules = _shieldD_NumberOfModules;

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldD_Offset[0],_shieldD_Offset[1],_shieldD_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;
    }

    // no good check against overlaps yet

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::halfpi);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::halfpi);

    _crs->_scintillatorShields[name] =
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset,
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    name = "CRSScintillatorTSRShield";

    numberOfModules = _shieldTSR_NumberOfModules;

    // the shields are placed in a "nominal" postion first

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldTSR_Offset[0],_shieldTSR_Offset[1],_shieldTSR_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;
    }

    // no good check against overlaps yet

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::halfpi);
    CRSScintillatorShieldRotationAngles.push_back(0.);

    // the constructors are "simple", most work is done in the maker

    _crs->_scintillatorShields[name] =
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset, // in Mu2e
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorTSLShield";

    numberOfModules = _shieldTSL_NumberOfModules;

    // the shields are placed in a "nominal" postion first

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldTSL_Offset[0],_shieldTSL_Offset[1],_shieldTSL_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;
    }

    // no good check against overlaps yet

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::pi);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::halfpi);
    CRSScintillatorShieldRotationAngles.push_back(0);

    // the constructors are "simple", most work is done in the maker

    _crs->_scintillatorShields[name] =
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset, // in Mu2e
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //

    name = "CRSScintillatorTSTShield";

    numberOfModules = _shieldTST_NumberOfModules;

    // the shields are placed in a "nominal" postion first

    CRSScintillatorShieldOffset =
      CLHEP::Hep3Vector(_shieldTST_Offset[0],_shieldTST_Offset[1],_shieldTST_Offset[2]);

    if ( _diagLevel > 0) {
      cout << __func__ << " CRSScintillatorShieldOffset : " << name << " : " <<
        CRSScintillatorShieldOffset << endl;
    }

    // no good check against overlaps yet

    // around x, y, z
    CRSScintillatorShieldRotationAngles.clear();
    CRSScintillatorShieldRotationAngles.reserve(3);
    CRSScintillatorShieldRotationAngles.push_back(0.);
    CRSScintillatorShieldRotationAngles.push_back(-CLHEP::halfpi);
    CRSScintillatorShieldRotationAngles.push_back(CLHEP::halfpi);

    // the constructors are "simple", most work is done in the maker

    _crs->_scintillatorShields[name] =
      CRSScintillatorShield(ishield,
                            name,
                            CRSScintillatorShieldRotationAngles,
                            CRSScintillatorShieldOffset + _scintillatorShieldOffset, // in Mu2e
                            _scintillatorShieldHalfThickness,
                            numberOfModules);

    ++ishield;

    //


    for (std::map<std::string,CRSScintillatorShield>::iterator itshield=_crs->_scintillatorShields.begin();
         itshield!=_crs->_scintillatorShields.end(); ++itshield) {

      if ( _diagLevel > 0) {
        cout << __func__ << " shield._name        : " << (itshield->second)._name << endl;
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

      // calculate shields position in the mu2e frame

      CLHEP::HepRotationX RX(shield._globalRotationAngles[0]);
      CLHEP::HepRotationY RY(shield._globalRotationAngles[1]);
      CLHEP::HepRotationZ RZ(shield._globalRotationAngles[2]);

      CLHEP::HepRotation shieldRotation(RX*RY*RZ);

      CLHEP::Hep3Vector shieldGlobalCRSOffsetRotated =
        (shieldRotation * (shield._globalOffset - _scintillatorShieldOffset)) + _scintillatorShieldOffset;

      cout << __func__ << " shieldGlobalCRSOffsetRotated : " << shield._name << " : "
           << shieldGlobalCRSOffsetRotated << endl;

      if (shield._name.compare("CRSScintillatorDShield")==0 ||
          shield._name.compare("CRSScintillatorUShield")==0) {
        cout << __func__ << " " << shield._name << " extent    : " <<
          shieldGlobalCRSOffsetRotated[CLHEP::Hep3Vector::Z] - shield._halfThickness << ", " <<
          shieldGlobalCRSOffsetRotated[CLHEP::Hep3Vector::Z] + shield._halfThickness << endl;
      } else {
        cout << __func__ << " " << shield._name << " extent    : " <<
          shieldGlobalCRSOffsetRotated[CLHEP::Hep3Vector::Z] - shieldHalfWidth << ", " <<
          shieldGlobalCRSOffsetRotated[CLHEP::Hep3Vector::Z] + shieldHalfWidth << endl;
      }
    }

    // module counter for the given shield
    int imodule = 0;
    // we place/make the full modules first

    double firstFullModuleLocalOffsetL = -shieldHalfWidth + _scintillatorFullModuleHalfWidth;

    for (int ii = 0; ii<numberOfFullModules; ++ii) {

      if ( _diagLevel > 1) {
        cout << __func__ << " creating module          : " << imodule << endl;
      }

      // module local offsets distance from the center of the shield for a given module

      double moduleLocalOffsetL = firstFullModuleLocalOffsetL +
        ii*(2.*_scintillatorFullModuleHalfWidth-_scintillatorModuleOverlap);

      if ( _diagLevel > 1) {
        cout << __func__ << " moduleLocalOffsetL               : " << moduleLocalOffsetL << endl;
      }

      // creating an empty module (need to know if full/half)
      // than modifying it

      shield._modules.push_back(CRSScintillatorModule( CRSScintillatorModuleId(shield._id,imodule),
                                                       _scintillatorBarsPerFullLayer));

      CRSScintillatorModule& module = shield._modules.back();

      module._layers.reserve(_scintillatorLayersPerModule);

/*
      // we need to shift the modules transversly as well
      // we shall assume that the first (0th) one will be touching the steel (%2 below)
      module._globalOffset = CLHEP::Hep3Vector( imodule%2 != 0 ?
                                               _scintillatorModuleCoreHalfThickness :
                                               -_scintillatorModuleCoreHalfThickness,
                                               0., moduleLocalOffsetL)
	+ shield._globalOffset;
*/
      module._globalOffset = CLHEP::Hep3Vector(0., 0., moduleLocalOffsetL)
	                   + shield._globalOffset;

      module._globalRotationAngles = shield._globalRotationAngles;

      ++imodule;

    }

    // same for half modules (may create a function for the two)
/*
    double firstHalfModuleLocalOffsetL = -shieldHalfWidth +
      numberOfFullModules*(2.*_scintillatorFullModuleHalfWidth - _scintillatorModuleOverlap)
      + _scintillatorHalfModuleHalfWidth;

    if ( _diagLevel > 1) {
      cout << __func__ << " _scintillatorModuleOverlap       : " << _scintillatorModuleOverlap << endl;
    }

    for (int ii = 0; ii<numberOfHalfModules; ++ii) {

      if ( _diagLevel > 1) {
        cout << __func__ << " creating module          : " << imodule << endl;
      }

      double moduleLocalOffsetL = firstHalfModuleLocalOffsetL +
        ii*(2.0*_scintillatorHalfModuleHalfWidth-_scintillatorModuleOverlap);

      if ( _diagLevel > 1) {
        cout << __func__ << " moduleLocalOffsetL               : " << moduleLocalOffsetL << endl;
      }

      // creating an empty module than modifying it

      shield._modules.push_back(CRSScintillatorModule( CRSScintillatorModuleId(shield._id,imodule),
                                                       _scintillatorBarsPerFullLayer/2 ));

      CRSScintillatorModule& module = shield._modules.back();

      module._layers.reserve(_scintillatorLayersPerModule);

      module._globalOffset = CLHEP::Hep3Vector(0., 0., moduleLocalOffsetL) 
	                   + shield._globalOffset;

      module._globalOffset = CLHEP::Hep3Vector( imodule%2 != 0 ?
                                               _scintillatorModuleCoreHalfThickness :
                                               -_scintillatorModuleCoreHalfThickness,
                                               0., moduleLocalOffsetL) 
	+ shield._globalOffset;


      module._globalRotationAngles = shield._globalRotationAngles;

      ++imodule;

    }
*/
    // set global offsets, make layers

    for (int ii = 0; ii<numberOfModules; ++ii) {

      if ( _diagLevel > 1) {
        cout << __func__ << " making layers for module : " << ii << endl;
      }

      CRSScintillatorModule& module = shield._modules.at(ii);
      if ( _diagLevel > 1) {
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

      if ( _diagLevel > 1) {
        cout << __func__ << " making layer : " << ii << endl;
      }

      module._layers.push_back(CRSScintillatorLayer( CRSScintillatorLayerId(module._id,ii),
                                                     numberOfBarsPerLayer));
      CRSScintillatorLayer& layer = module._layers.back();

      // fill in the layer offsets etc...

      // we will place innermost layers shifted most to the left
      layer._globalOffset = CLHEP::Hep3Vector(firtstLayerLocalOffsetT + ii*scintillatorLayerShiftT,
					      0.,
					      firtstLayerLocalOffsetL + ii*_scintillatorLayerShift)
	+ module._globalOffset;

      layer._globalRotationAngles = module._globalRotationAngles;

      layer._nBars = module._nBarsPerLayer;

      if ( _diagLevel > 2) {
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

      if ( _diagLevel > 2) {
        cout << __func__ << " making bar   : " << ii << endl;
      }

      // calculate the bar offsets, ids, indeces and enter it into the container

      CRSScintillatorBarIndex index(_crs->_allCRSScintillatorBars.size());

      // local longitudinal offset wrt the center of the layer

      _crs->_allCRSScintillatorBars.push_back(CRSScintillatorBar(CRSScintillatorBarId(layer._id,ii),
                                                                 index));

      CRSScintillatorBar& bar = _crs->_allCRSScintillatorBars.back();

      CLHEP::Hep3Vector barLocalOffset(0.,0.,
				       firtstBarLocalOffsetL + ii*(2.*_scintillatorBarHalfLengths[2]+
								   _scintillatorLayerGap));

      bar._globalRotationAngles = layer._globalRotationAngles;

      // untill now the object positions were calculated locally
      // wrt to the the module etc... in a "nominal" position

      // we will now calculate position of the bar wrt the CRS origin and rotate it wrt that origin

      CLHEP::Hep3Vector barCRSOffset = layer._globalOffset + barLocalOffset - _scintillatorShieldOffset;

      CLHEP::HepRotationX RX(bar._globalRotationAngles[0]);
      CLHEP::HepRotationY RY(bar._globalRotationAngles[1]);
      CLHEP::HepRotationZ RZ(bar._globalRotationAngles[2]);

      CLHEP::HepRotation barRotation(RX*RY*RZ);

      CLHEP::Hep3Vector barCRSOffsetRotated = barRotation * barCRSOffset;

      //  bar._globalOffset = layer._globalOffset + bar._localOffset;

      //  note that difference, see the comments above, rotation done for the bars only so far

      bar._globalOffset = barCRSOffsetRotated + _scintillatorShieldOffset;

      layer._bars.push_back(&bar);
      bar._detail = &_crs->_barDetails;
      layer._indices.push_back(index);

      if ( _diagLevel > 3) {
        cout << __func__ << " barCRSOffset        : " << barCRSOffset        << endl;
        cout << __func__ << " barRotation         : " << barRotation         << endl;
        cout << __func__ << " barCRSOffsetRotated : " << barCRSOffsetRotated << endl;
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
      throw cet::exception("GEOM")
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

    // we need to allocate space for all the bars at once to avoid pointer invalidatation

    _totalNumberOfBars =
      (_shieldR_NumberOfModules[0] +
       _shieldL_NumberOfModules[0] +
       _shieldD_NumberOfModules[0] +
       _shieldT_NumberOfModules[0] +
       _shieldTSR_NumberOfModules[0] +
       _shieldTSL_NumberOfModules[0] +
       _shieldTST_NumberOfModules[0])*_scintillatorLayersPerModule*_scintillatorBarsPerFullLayer +
      (_shieldR_NumberOfModules[1] +
       _shieldL_NumberOfModules[1] +
       _shieldD_NumberOfModules[1] +
       _shieldT_NumberOfModules[1] +
       _shieldTSR_NumberOfModules[1] +
       _shieldTSL_NumberOfModules[1] +
       _shieldTST_NumberOfModules[1])*_scintillatorLayersPerModule*(_scintillatorBarsPerFullLayer/2);

    if ( _diagLevel > 0) {
      cout << __func__ << " _totalNumberOfBars : " << _totalNumberOfBars << endl;
    }
    _crs->_allCRSScintillatorBars.reserve(_totalNumberOfBars);

  }

  void CosmicRayShieldMaker::makeCRSSteelShield(SimpleConfig const & _config) {
    // first make the steel (fluxreturn)

    double HallSteelDSShieldHalfWidth=(_HallSteelLShieldCenter.x()-_HallSteelRShieldCenter.x())*0.5+_HallSteelHalfThick;
    double HallSteelTSTShieldHalfWidth=(_HallSteelTSLShieldCenter.z()-_HallSteelTSRShieldCenter.z())*0.5+_HallSteelHalfThick;

    double CRSSteelRShieldDims[3] ={_HallSteelHalfThick,_HallSteelHalfSideShieldHeight,_HallSteelHalfRShieldLength};
    double CRSSteelLShieldDims[3] ={_HallSteelHalfThick,_HallSteelHalfSideShieldHeight,_HallSteelHalfLShieldLength};
    double CRSSteelTShieldDims[3] ={HallSteelDSShieldHalfWidth,_HallSteelHalfThick,_HallSteelHalfTShieldLength};
    double CRSSteelDShieldDims[3] ={HallSteelDSShieldHalfWidth,_HallSteelHalfSideShieldHeight,_HallSteelHalfThick};
    double CRSSteelTSRShieldDims[3] ={_HallSteelHalfTSRShieldLength,_HallSteelHalfSideShieldHeight,_HallSteelHalfThick};
    double CRSSteelTSLShieldDims[3] ={_HallSteelHalfTSLShieldLength,_HallSteelHalfSideShieldHeight,_HallSteelHalfThick};
    double CRSSteelTSTShieldDims[3] ={_HallSteelHalfTSTShieldLength,_HallSteelHalfThick,HallSteelTSTShieldHalfWidth};

    // finaly create the steel shield objects (we invent/assign their names here...)

    std::string name = "CRSSteelRShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelRShieldCenter + _HallSteelOffset,CRSSteelRShieldDims);
    //                                         HepRotation,    global offset in Mu2e

    name = "CRSSteelLShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelLShieldCenter + _HallSteelOffset,CRSSteelLShieldDims);

    name = "CRSSteelTShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelTShieldCenter + _HallSteelOffset,CRSSteelTShieldDims);

    double downStreamHoleRadius = _config.getBool("hasMBS",false) ?
      GeomHandle<MBS>()->getEnvelopeRmax() : 0. ;

    name = "CRSSteelDShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelDShieldCenter + _HallSteelOffset,CRSSteelDShieldDims,
                                               downStreamHoleRadius, CLHEP::Hep3Vector(0.,-_HallSteelDShieldCenter.y(),0.));

    name = "CRSSteelTSRShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelTSRShieldCenter + _HallSteelOffset,CRSSteelTSRShieldDims);

    name = "CRSSteelTSLShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelTSLShieldCenter + _HallSteelOffset,CRSSteelTSLShieldDims);

    name = "CRSSteelTSTShield";
    _crs->_steelShields[name] = CRSSteelShield(name,0,_HallSteelTSTShieldCenter + _HallSteelOffset,CRSSteelTSTShieldDims);
  }

} // namespace mu2e
