//
// Construct and return Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
//
// Notes
// See mu2e-doc-XXXX for naming conventions etc.

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "GeometryService/inc/STMMaker.hh"
#include "STMGeom/inc/STM.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  STMMaker::STMMaker(SimpleConfig const & _config,
                     double solenoidOffset)
  {
    // if( ! _config.getBool("hasSTM",false) ) return;

    // create an empty STM
    _stm = unique_ptr<STM>(new STM());

    // access its object through a reference

    STM & stm = *_stm.get();

    parseConfig(_config);

    // now create the specific components

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    const CLHEP::Hep3Vector &dsP( ds->position() );    

    GeomHandle<CosmicRayShield> CRS;
    std::vector<double> crvd_halfLengths = CRS->getSectorHalfLengths("D");
    CLHEP::Hep3Vector   crvd_position    = CRS->getSectorPosition("D");
    const double z_crv_max = crvd_position.z() + crvd_halfLengths[2];
    //const double z_crv_max = 20000.0;
    
    //Create a reference position (most things in the STM geometry will be defined w.r.t. this position)
    // Our reference z is the downstream edge of the CRV-D
    const CLHEP::Hep3Vector _STMMOffsetInMu2e(dsP.x(), 0.0, z_crv_max );
    const CLHEP::HepRotation _magnetRotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector _magnetOffsetInMu2e  = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);
    //if (_magnetBuild){
      stm._pSTMMagnetParams = std::unique_ptr<PermanentMagnet>
        (new PermanentMagnet(_magnetBuild,
                             _magnetHalfWidth,
                             _magnetHalfHeight,
                             _magnetHalfLength,
                             _magnetHoleHalfWidth,
                             _magnetHoleHalfHeight,
                             _magnetOffsetInMu2e,
                             _magnetRotation,
                             _magnetMaterial,
                             _magnetField,
                             _magnetFieldVisible
                            ));
    //}

    
    const CLHEP::HepRotation _FOVCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _FOVCollimatorOffsetInMu2e = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+2.0*_magnetHalfLength+2.0*_pipeDnStrHalfLength+_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength);    
    //if (_FOVCollimatorBuild){
      stm._pSTMFOVCollimatorParams = std::unique_ptr<STMCollimator>
        (new STMCollimator(_FOVCollimatorBuild,           
                           _FOVCollimatorHalfWidth,       
                           _FOVCollimatorHalfHeight,      
                           _FOVCollimatorHalfLength,      
                           _FOVCollimatorLinerBuild,  
                           _FOVCollimatorLinerHalfWidth,  
                           _FOVCollimatorLinerHalfHeight, 
                           _FOVCollimatorLinerHalfLength, 
                           _FOVCollimatorLinerCutOutHalfLength,
                           _FOVCollimatorHole1xOffset,    
                           _FOVCollimatorHole1RadiusUpStr,
                           _FOVCollimatorHole1RadiusDnStr,
                           _FOVCollimatorHole1LinerBuild,
                           _FOVCollimatorHole1LinerThickness,
                           _FOVCollimatorHole2Build,      
                           _FOVCollimatorHole2xOffset,    
                           _FOVCollimatorHole2RadiusUpStr,
                           _FOVCollimatorHole2RadiusDnStr,
                           _FOVCollimatorHole2LinerBuild,
                           _FOVCollimatorHole2LinerThickness,
                           _FOVCollimatorOffsetInMu2e,
                           _FOVCollimatorRotation,
                           _FOVCollimatorMaterial,        
                           _FOVCollimatorLinerMaterial,
                           _FOVCollimatorHoleLinerMaterial
                          ));       
    //}
    
    
    const CLHEP::HepRotation _pipeRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _pipeOffsetInMu2e = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);    
    //if (_pipeBuild){
      stm._pSTMTransportPipeParams = std::unique_ptr<TransportPipe>
        (new TransportPipe(_pipeBuild,                
                           _pipeRadiusIn,             
                           _pipeRadiusOut,            
                           _pipeMaterial,             
                           _pipeGasMaterial,          
                           _pipeUpStrSpace,           
                           _pipeDnStrHalfLength,      
                           _pipeUpStrWindowMaterial,  
                           _pipeUpStrWindowHalfLength,
                           _pipeDnStrWindowMaterial,  
                           _pipeDnStrWindowHalfLength,
                           _pipeFlangeHalfLength,     
                           _pipeFlangeOverhangR,
                           _pipeOffsetInMu2e,
                           _pipeRotation
                          ));       
    //}

    double _magnetTableTopHalfWidth = 0.0;
    if ( _magnetBuild && !_FOVCollimatorBuild) _magnetTableTopHalfWidth = _magnetHalfWidth;
    if (!_magnetBuild &&  _FOVCollimatorBuild) _magnetTableTopHalfWidth = _FOVCollimatorHalfWidth;
    if ( _magnetBuild &&  _FOVCollimatorBuild) _magnetTableTopHalfWidth = std::max(_magnetHalfWidth,_FOVCollimatorHalfWidth);
    _magnetTableTopHalfWidth += _magnetTableTopExtraWidth;
    
    double _magnetTableTopHalfLength = 0.0;
    if ( _magnetBuild )        _magnetTableTopHalfLength += _magnetHalfLength;
    if ( _pipeBuild )          _magnetTableTopHalfLength += _pipeDnStrHalfLength;
    if ( _FOVCollimatorBuild ) _magnetTableTopHalfLength += 0.5*_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength;
    if ( _shieldBuild )        _magnetTableTopHalfLength += 0.5*_shieldDnStrSpace+_shieldDnStrWallHalfLength;
    _magnetTableTopHalfLength += _magnetTableTopExtraLength;
    
    const CLHEP::HepRotation _magnetTableRotation     = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _magnetTableOffsetInMu2e = _STMMOffsetInMu2e - CLHEP::Hep3Vector(0.0,_magnetHalfHeight+_magnetTableTopHalfHeight,0.0);
    if ( _magnetBuild )        _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_magnetUpStrSpace+_magnetHalfLength);    
    if ( _pipeBuild )          _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,_pipeDnStrHalfLength);    
    if ( _FOVCollimatorBuild ) _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,0.5*_FOVCollimatorUpStrSpace+_FOVCollimatorHalfLength);    
    if ( _shieldBuild )        _magnetTableOffsetInMu2e += CLHEP::Hep3Vector(0.0,0.,-0.5*_shieldDnStrSpace-_shieldDnStrWallHalfLength);    
    
    //if (_magnetTableBuild && (_magnetBuild||_FOVCollimatorBuild) ){
      stm._pSTMMagnetSupportTableParams = std::unique_ptr<SupportTable>
        (new SupportTable( _magnetTableBuild,
                           _magnetTableTopHalfWidth,
                           _magnetTableTopHalfHeight,
                           _magnetTableTopHalfLength,
                           _magnetTableLegRadius,
                           _magnetTableOffsetInMu2e,
                           _magnetTableRotation,
                           _magnetTableMaterial
                          ));        
    //}
    
    
    //The STM geometry must fit inside the detector hall, so find the z of the East hall wall
    GeomHandle<Mu2eHall> hall;
    const double z_hall_inside_max = hall->getWallExtentz("dsArea",1)/CLHEP::mm;//the integer allows you to specify which side of which wall you want the z for: 1 = west side of east wall (i.e. the z of the inside surface of the east wall)
    const CLHEP::Hep3Vector BeamAxisAtEastWallInMu2e(dsP.x(), 0.0, z_hall_inside_max );
    
    const CLHEP::HepRotation _SSCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _SSCollimatorOffsetInMu2e = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0.0,0.,-_stmZAllowed+_SSCollimatorHalfLength);    
    //if (_SSCollimatorBuild){
      stm._pSTMSSCollimatorParams = std::unique_ptr<STMCollimator>
        (new STMCollimator(_SSCollimatorBuild,           
                           _SSCollimatorHalfWidth,       
                           _SSCollimatorHalfHeight,      
                           _SSCollimatorHalfLength,      
                           _SSCollimatorLinerBuild,  
                           _SSCollimatorLinerHalfWidth,  
                           _SSCollimatorLinerHalfHeight, 
                           _SSCollimatorLinerHalfLength,
                           _SSCollimatorLinerCutOutHalfLength, 
                           _SSCollimatorHole1xOffset,    
                           _SSCollimatorHole1RadiusUpStr,
                           _SSCollimatorHole1RadiusDnStr,
                           _SSCollimatorHole1LinerBuild,
                           _SSCollimatorHole1LinerThickness,
                           _SSCollimatorHole2Build,      
                           _SSCollimatorHole2xOffset,    
                           _SSCollimatorHole2RadiusUpStr,
                           _SSCollimatorHole2RadiusDnStr,
                           _SSCollimatorHole2LinerBuild,
                           _SSCollimatorHole2LinerThickness,
                           _SSCollimatorOffsetInMu2e,
                           _SSCollimatorRotation,
                           _SSCollimatorMaterial,        
                           _SSCollimatorLinerMaterial,
                           _SSCollimatorHoleLinerMaterial
                          ));       
    //}
     

    double _detectorTableTopHalfWidth = _SSCollimatorHalfWidth + _detectorTableTopExtraWidth;    
    double _detectorTableTopHalfLength = 0.5*_stmZAllowed - 1.0;
    const CLHEP::HepRotation _detectorTableRotation = CLHEP::HepRotation::IDENTITY;
    CLHEP::Hep3Vector  _detectorTableOffsetInMu2e = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0.0,-_SSCollimatorHalfHeight-_detectorTableTopHalfHeight, -_stmZAllowed+_detectorTableTopHalfLength);    
     
    //if (_detectorTableBuild && _SSCollimatorBuild ){
      stm._pSTMDetectorSupportTableParams = std::unique_ptr<SupportTable>
        (new SupportTable( _detectorTableBuild,
                           _detectorTableTopHalfWidth,
                           _detectorTableTopHalfHeight,
                           _detectorTableTopHalfLength,
                           _detectorTableLegRadius,
                           _detectorTableOffsetInMu2e,
                           _detectorTableRotation,
                           _detectorTableMaterial
                          ));        
    //}
    
                      
    const CLHEP::HepRotation _detector1Rotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _detector1OffsetInMu2e  = _SSCollimatorOffsetInMu2e + CLHEP::Hep3Vector(_detector1xOffset, 0.0, _SSCollimatorHalfLength+_detector1CanUpStrSpace+_detector1CanHalfLength);
    //if (_detector1Build){
      stm._pSTMDetector1Params = std::unique_ptr<GeDetector>
        (new GeDetector(_detector1Build,
                        _detector1CrystalMaterial,
                        _detector1CrystalRadiusIn,
                        _detector1CrystalRadiusOut,
                        _detector1CrystalHalfLength,
                        _detector1CanMaterial,
                        _detector1CanRadiusIn,
                        _detector1CanRadiusOut,
                        _detector1CanHalfLength,
                        _detector1CanUpStrWindowMaterial,
                        _detector1CanUpStrWindowHalfLength,
                        _detector1CanGasMaterial,
                        _detector1OffsetInMu2e,
                        _detector1Rotation
                       ));
    //}
    
    const CLHEP::HepRotation _detector2Rotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _detector2OffsetInMu2e  = _SSCollimatorOffsetInMu2e + CLHEP::Hep3Vector(_detector2xOffset, 0.0, _SSCollimatorHalfLength+_detector2CanUpStrSpace+_detector2CanHalfLength);
    //if (_detector2Build){
      stm._pSTMDetector2Params = std::unique_ptr<GeDetector>
        (new GeDetector(_detector2Build,
                        _detector2CrystalMaterial,
                        _detector2CrystalRadiusIn,
                        _detector2CrystalRadiusOut,
                        _detector2CrystalHalfLength,
                        _detector2CanMaterial,
                        _detector2CanRadiusIn,
                        _detector2CanRadiusOut,
                        _detector2CanHalfLength,
                        _detector2CanUpStrWindowMaterial,
                        _detector2CanUpStrWindowHalfLength,
                        _detector2CanGasMaterial,
                        _detector2OffsetInMu2e,
                        _detector2Rotation
                       ));
    //} 
    

    const CLHEP::HepRotation _shieldRotation = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _shieldOffsetInMu2e  = _FOVCollimatorOffsetInMu2e + CLHEP::Hep3Vector(0.0, 0.0, -_FOVCollimatorHalfLength);
    //if (_shieldBuild){
      stm._pSTMShieldPipeParams = std::unique_ptr<ShieldPipe>
        (new ShieldPipe(_shieldBuild,
                        _shieldRadiusIn,
                        _shieldLinerWidth,
                        _shieldRadiusOut,
                        _shieldPipeHalfLength,
                        _shieldMaterialLiner,
                        _shieldMaterial,
                        _shieldUpStrSpace,
                        _shieldDnStrSpace,
                        _shieldDnStrWallHalfLength,
                        _shieldOffsetInMu2e, //This is upstream edge of FOV collimator for now.
                        _shieldRotation
                       ));
    //}
    
  
  }

  void STMMaker::parseConfig( SimpleConfig const & _config ){

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "stmMagnetField", "stm.magnet.field");


    _verbosityLevel            = _config.getInt("stm.verbosityLevel",0);
    _stmZAllowed               = _config.getDouble("stm.z.allowed");
    
    _magnetBuild               = _config.getBool(  "stm.magnet.build",false);
    _magnetUpStrSpace          = _config.getDouble("stm.magnet.UpStrSpace");
    _magnetHalfLength          = _config.getDouble("stm.magnet.halfLength");
    _magnetHalfWidth           = _config.getDouble("stm.magnet.halfWidth");
    _magnetHalfHeight          = _config.getDouble("stm.magnet.halfHeight");
    _magnetHoleHalfWidth       = _config.getDouble("stm.magnet.holeHalfWidth");
    _magnetHoleHalfHeight      = _config.getDouble("stm.magnet.holeHalfHeight");    
    _magnetMaterial            = _config.getString("stm.magnet.material");
    _magnetField               = _config.getDouble("stm.magnet.field");
    //_magnetFieldVisible        = _config.getBool(  "stm.magnet.fieldVisible",false);
    _magnetFieldVisible        = geomOptions->isVisible("stmMagnetField");
    
    _FOVCollimatorBuild            = _config.getBool(  "stm.FOVcollimator.build");
    _FOVCollimatorMaterial         = _config.getString("stm.FOVcollimator.material");       
    _FOVCollimatorUpStrSpace       = _config.getDouble("stm.FOVcollimator.UpStrSpace");
    _FOVCollimatorHalfWidth        = _config.getDouble("stm.FOVcollimator.halfWidth");
    _FOVCollimatorHalfHeight       = _config.getDouble("stm.FOVcollimator.halfHeight");
    _FOVCollimatorHalfLength       = _config.getDouble("stm.FOVcollimator.halfLength");
    _FOVCollimatorLinerBuild       = _config.getBool(  "stm.FOVcollimator.liner.build");
    _FOVCollimatorLinerMaterial    = _config.getString("stm.FOVcollimator.liner.material");       
    _FOVCollimatorLinerHalfWidth   = _config.getDouble("stm.FOVcollimator.liner.halfWidth");
    _FOVCollimatorLinerHalfHeight  = _config.getDouble("stm.FOVcollimator.liner.halfHeight");
    _FOVCollimatorLinerHalfLength  = _config.getDouble("stm.FOVcollimator.liner.halfLength");
    _FOVCollimatorLinerCutOutHalfLength  = _config.getDouble("stm.FOVcollimator.liner.cutOutHalfLength");
    _FOVCollimatorHole1xOffset     = _config.getDouble("stm.FOVcollimator.hole1.xoffset");
    _FOVCollimatorHole1RadiusUpStr = _config.getDouble("stm.FOVcollimator.hole1.radiusUpStr");
    _FOVCollimatorHole1RadiusDnStr = _config.getDouble("stm.FOVcollimator.hole1.radiusDnStr");
    _FOVCollimatorHole1LinerBuild     = _config.getBool(  "stm.FOVcollimator.hole1.liner.build");
    _FOVCollimatorHole1LinerThickness = _config.getDouble("stm.FOVcollimator.hole1.liner.thickness");
    _FOVCollimatorHole2Build       = _config.getBool(  "stm.FOVcollimator.hole2.build");
    _FOVCollimatorHole2xOffset     = _config.getDouble("stm.FOVcollimator.hole2.xoffset");
    _FOVCollimatorHole2RadiusUpStr = _config.getDouble("stm.FOVcollimator.hole2.radiusUpStr");
    _FOVCollimatorHole2RadiusDnStr = _config.getDouble("stm.FOVcollimator.hole2.radiusDnStr");
    _FOVCollimatorHole2LinerBuild     = _config.getBool(  "stm.FOVcollimator.hole2.liner.build");
    _FOVCollimatorHole2LinerThickness = _config.getDouble("stm.FOVcollimator.hole2.liner.thickness");
    _FOVCollimatorHoleLinerMaterial= _config.getString("stm.FOVcollimator.hole.liner.material");
    
    _pipeBuild                 = _config.getBool(  "stm.pipe.build");
    _pipeRadiusIn              = _config.getDouble("stm.pipe.rIn");
    _pipeRadiusOut             = _config.getDouble("stm.pipe.rOut");
    _pipeMaterial              = _config.getString("stm.pipe.material");
    _pipeGasMaterial           = _config.getString("stm.pipe.gas.material");
    _pipeUpStrSpace            = _config.getDouble("stm.pipe.UpStrSpace");
    _pipeDnStrHalfLength       = _config.getDouble("stm.pipe.DnStrHalfLength");
    _pipeUpStrWindowMaterial   = _config.getString("stm.pipe.UpStrWindow.material");
    _pipeUpStrWindowHalfLength = _config.getDouble("stm.pipe.UpStrWindow.halfLength");
    _pipeDnStrWindowMaterial   = _config.getString("stm.pipe.DnStrWindow.material");
    _pipeDnStrWindowHalfLength = _config.getDouble("stm.pipe.DnStrWindow.halfLength");
    _pipeFlangeHalfLength      = _config.getDouble("stm.pipe.flange.halfLength");
    _pipeFlangeOverhangR       = _config.getDouble("stm.pipe.flange.overhangR");
    
    _magnetTableBuild          = _config.getBool(  "stm.magnet.stand.build",false);
    _magnetTableMaterial       = _config.getString("stm.magnet.stand.material");       
    _magnetTableTopExtraWidth  = _config.getDouble("stm.magnet.stand.topExtraWidth");
    _magnetTableTopExtraLength = _config.getDouble("stm.magnet.stand.topExtraLength");
    _magnetTableTopHalfHeight  = _config.getDouble("stm.magnet.stand.topHalfHeight");
    _magnetTableLegRadius      = _config.getDouble("stm.magnet.stand.legRadius");
    
    _SSCollimatorBuild            = _config.getBool(  "stm.SScollimator.build");
    _SSCollimatorMaterial         = _config.getString("stm.SScollimator.material");       
    _SSCollimatorUpStrSpace       = _config.getDouble("stm.SScollimator.UpStrSpace");
    _SSCollimatorHalfWidth        = _config.getDouble("stm.SScollimator.halfWidth");
    _SSCollimatorHalfHeight       = _config.getDouble("stm.SScollimator.halfHeight");
    _SSCollimatorHalfLength       = _config.getDouble("stm.SScollimator.halfLength");
    _SSCollimatorLinerBuild       = _config.getBool(  "stm.SScollimator.liner.build");
    _SSCollimatorLinerMaterial    = _config.getString("stm.SScollimator.liner.material");       
    _SSCollimatorLinerHalfWidth   = _config.getDouble("stm.SScollimator.liner.halfWidth");
    _SSCollimatorLinerHalfHeight  = _config.getDouble("stm.SScollimator.liner.halfHeight");
    _SSCollimatorLinerHalfLength  = _config.getDouble("stm.SScollimator.liner.halfLength");
    _SSCollimatorLinerCutOutHalfLength  = _config.getDouble("stm.SScollimator.liner.cutOutHalfLength");
    _SSCollimatorHole1xOffset     = _config.getDouble("stm.SScollimator.hole1.xoffset");
    _SSCollimatorHole1RadiusUpStr = _config.getDouble("stm.SScollimator.hole1.radiusUpStr");
    _SSCollimatorHole1RadiusDnStr = _config.getDouble("stm.SScollimator.hole1.radiusDnStr");
    _SSCollimatorHole1LinerBuild     = _config.getBool(  "stm.SScollimator.hole1.liner.build");
    _SSCollimatorHole1LinerThickness = _config.getDouble("stm.SScollimator.hole1.liner.thickness");
    _SSCollimatorHole2Build       = _config.getBool(  "stm.SScollimator.hole2.build");
    _SSCollimatorHole2xOffset     = _config.getDouble("stm.SScollimator.hole2.xoffset");
    _SSCollimatorHole2RadiusUpStr = _config.getDouble("stm.SScollimator.hole2.radiusUpStr");
    _SSCollimatorHole2RadiusDnStr = _config.getDouble("stm.SScollimator.hole2.radiusDnStr");    
    _SSCollimatorHole2LinerBuild     = _config.getBool(  "stm.SScollimator.hole2.liner.build");
    _SSCollimatorHole2LinerThickness = _config.getDouble("stm.SScollimator.hole2.liner.thickness");
    _SSCollimatorHoleLinerMaterial = _config.getString("stm.SScollimator.hole.liner.material");       
    
    _detectorTableBuild          = _config.getBool(  "stm.detector.stand.build",false);
    _detectorTableMaterial       = _config.getString("stm.detector.stand.material");       
    _detectorTableTopExtraWidth  = _config.getDouble("stm.detector.stand.topExtraWidth");
    _detectorTableTopExtraLength = _config.getDouble("stm.detector.stand.topExtraLength");
    _detectorTableTopHalfHeight  = _config.getDouble("stm.detector.stand.topHalfHeight");
    _detectorTableLegRadius      = _config.getDouble("stm.detector.stand.legRadius");
    
    _detector1Build                    = _config.getBool(  "stm.det1.build",false);
    _detector1CrystalMaterial          = _config.getString("stm.det1.material");   
    _detector1CrystalRadiusIn          = _config.getDouble("stm.det1.rIn");
    _detector1CrystalRadiusOut         = _config.getDouble("stm.det1.rOut");
    _detector1CrystalHalfLength        = _config.getDouble("stm.det1.halfLength");
    _detector1xOffset                  = _config.getDouble("stm.det1.xoffset");
    _detector1CanMaterial              = _config.getString("stm.det1.can.material"); 
    _detector1CanRadiusIn              = _config.getDouble("stm.det1.can.rIn");
    _detector1CanRadiusOut             = _config.getDouble("stm.det1.can.rOut");
    _detector1CanHalfLength            = _config.getDouble("stm.det1.can.halfLength");
    _detector1CanUpStrSpace            = _config.getDouble("stm.det1.can.UpStrSpace");
    _detector1CanUpStrWindowMaterial   = _config.getString("stm.det1.can.UpStrWindowMaterial"); 
    _detector1CanUpStrWindowHalfLength = _config.getDouble("stm.det1.can.UpStrWindowHalfLength");
    _detector1CanGasMaterial           = _config.getString("stm.det1.can.gas"); 

    _detector2Build                    = _config.getBool(  "stm.det2.build",false);
    _detector2CrystalMaterial          = _config.getString("stm.det2.material");   
    _detector2CrystalRadiusIn          = _config.getDouble("stm.det2.rIn");
    _detector2CrystalRadiusOut         = _config.getDouble("stm.det2.rOut");
    _detector2CrystalHalfLength        = _config.getDouble("stm.det2.halfLength");
    _detector2xOffset                  = _config.getDouble("stm.det2.xoffset");
    _detector2CanMaterial              = _config.getString("stm.det2.can.material"); 
    _detector2CanRadiusIn              = _config.getDouble("stm.det2.can.rIn");
    _detector2CanRadiusOut             = _config.getDouble("stm.det2.can.rOut");
    _detector2CanHalfLength            = _config.getDouble("stm.det2.can.halfLength");
    _detector2CanUpStrSpace            = _config.getDouble("stm.det2.can.UpStrSpace");
    _detector2CanUpStrWindowMaterial   = _config.getString("stm.det2.can.UpStrWindowMaterial"); 
    _detector2CanUpStrWindowHalfLength = _config.getDouble("stm.det2.can.UpStrWindowHalfLength");
    _detector2CanGasMaterial           = _config.getString("stm.det2.can.gas"); 
    
    _shieldBuild                = _config.getBool(  "stm.shield.build",false);
    _shieldRadiusIn             = _config.getDouble("stm.shield.rIn");
    _shieldLinerWidth           = _config.getDouble("stm.shield.widthLiner");
    _shieldRadiusOut            = _config.getDouble("stm.shield.rOut");
    _shieldPipeHalfLength       = _config.getDouble("stm.shield.pipe.halfLength");
    _shieldMaterialLiner        = _config.getString("stm.shield.materialLiner");
    _shieldMaterial             = _config.getString("stm.shield.material"); 
    _shieldUpStrSpace           = _config.getDouble("stm.shield.UpStrSpace");
    _shieldDnStrSpace           = _config.getDouble("stm.shield.DnStrSpace");
    _shieldDnStrWallHalfLength  = _config.getDouble("stm.shield.DnStrWall.halfLength");

    
  }

} // namespace mu2e
