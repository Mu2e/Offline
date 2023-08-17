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
#include "Offline/GeometryService/inc/STMMaker.hh"
#include "Offline/STMGeom/inc/STM.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"

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
    double magnetZOffset = _magnetUpStrSpace+_magnetHalfLength;
    //calculate the magnet position assuming the shield pipe is flush to the wall
    if(_config.getBool("stm.magnet.usePipeAsOrigin", false)) {
      magnetZOffset = _magnetHalfLength + _shieldDnStrWallGap + 2.*_shieldPipeHalfLength + _shieldUpStrWallGap;
      if(!_shieldMatchPipeBlock) magnetZOffset += 2.*_shieldDnStrWallHalfLength;
    }
    const CLHEP::Hep3Vector _magnetOffsetInMu2e  = _STMMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,magnetZOffset);
    const CLHEP::Hep3Vector _magnetHoleOffset  = CLHEP::Hep3Vector(_magnetHoleXOffset,_magnetHoleYOffset, 0.);
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
                             _magnetHoleOffset,
                             _magnetMaterial,
                             _magnetHasLiner,
                             _magnetField,
                             _magnetFieldVisible
                            ));
    //}


    const CLHEP::HepRotation _FOVCollimatorRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _FOVCollimatorOffsetInMu2e = _magnetOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_magnetHalfLength + _FOVCollimatorUpStrSpace + _FOVCollimatorHalfLength);
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


      ////////////////////////////////////
      // STM Downstream Area
      //
    //The STM geometry must fit inside the detector hall, so find the z of the East hall wall
    GeomHandle<Mu2eHall> hall;
    const double z_hall_inside_max = hall->getWallExtentz("dsArea",1)/CLHEP::mm;//the integer allows you to specify which side of which wall you want the z for: 1 = west side of east wall (i.e. the z of the inside surface of the east wall)
    const CLHEP::Hep3Vector BeamAxisAtEastWallInMu2e(dsP.x(), 0.0, z_hall_inside_max );
    const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
    const CLHEP::Hep3Vector FloorAtEastWallInMu2e = BeamAxisAtEastWallInMu2e - CLHEP::Hep3Vector(0.0, yExtentLow, 0.0);

    // Define the envelope w.r.t the floor at the east wall
    const CLHEP::HepRotation _stmDnStrEnvRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector _stmDnStrEnvPositionInMu2e = FloorAtEastWallInMu2e + CLHEP::Hep3Vector(0.0, +_stmDnStrEnvHalfHeight, -_stmDnStrEnvHalfLength);
    stm._pSTMDnStrEnvParams = std::unique_ptr<STMDownstreamEnvelope>
      (new STMDownstreamEnvelope(_stmDnStrEnvBuild,
                                 _stmDnStrEnvHalfWidth,
                                 _stmDnStrEnvHalfHeight,
                                 _stmDnStrEnvHalfLength,
                                 _stmDnStrEnvPositionInMu2e,
                                 _stmDnStrEnvRotation,
                                 _stmDnStrEnvMaterial
                                 ));

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
                        _shieldHasLiner,
                        _shieldLinerWidth,
                        _shieldRadiusOut,
                        _shieldPipeHalfLength,
                        _shieldMaterialLiner,
                        _shieldMaterial,
                        _shieldMatchPipeBlock,
                        _shieldUpStrSpace,
                        _shieldDnStrSpace,
                        _shieldDnStrWallHalfLength,
                        _shieldDnStrWallHoleRadius,
                        _shieldDnStrWallHalfHeight,
                        _shieldDnStrWallHalfWidth,
                        _shieldDnStrWallGap,
                        _shieldDnStrWallMaterial,
                        _shieldOffsetInMu2e, //This is upstream edge of FOV collimator for now.
                        _shieldRotation
                       ));
    //}


    double fW_x = -0.15*25.4 - _STM_SSCdelta_WlR/2 + _STM_SSCdelta_WlL/2;
    double fW_z_1 =  -stm._copperdepth - stm._leaddepth2 - stm._BPdepth*2 - stm._leaddepth -(_STM_SSCWdepth_b+_STM_SSCWdepth_f)/2 - _STM_SSCleak;

   
    const CLHEP::Hep3Vector  _STMShieldingRef = BeamAxisAtEastWallInMu2e + CLHEP::Hep3Vector(0., 0., -_stmZAllowed+(_STM_SSCWdepth_b+_STM_SSCWdepth_f)/2);
    const CLHEP::HepRotation _STM_SSCRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _STM_SSCOffsetInMu2e = _STMShieldingRef  + CLHEP::Hep3Vector(fW_x, 0, fW_z_1 - (_STM_SSCWdepth_b+_STM_SSCWdepth_f)/2);
 
          stm._pSTM_SSCParams = std::unique_ptr<STM_SSC>
                  (new STM_SSC(_STM_SSCBuild,
                               _STM_SSCdelta_WlR,
                               _STM_SSCdelta_WlL,
                               _STM_SSCW_length,
                               _STM_SSCW_height,
                               _STM_SSCleak,
                               _STM_SSCWdepth_f,
                               _STM_SSCWdepth_b,
                               _STM_SSCAperture_HPGe1,
                               _STM_SSCAperture_HPGe2,
                               _STM_SSCAperture_LaBr1,
                               _STM_SSCAperture_LaBr2,
                               _STM_SSCoffset_Spot,
                               _STM_SSCOffsetInMu2e,
                               _STM_SSCRotation,
                               _STM_SSCMaterial));
   
    const double offset_HPGe=_HPGeoffset_HPGe;

    CLHEP::HepRotation rotHPGe = CLHEP::HepRotation::IDENTITY;
    rotHPGe.rotateY(45*CLHEP::degree);

    const CLHEP::HepRotation _HPGeRotation     = rotHPGe;
    const CLHEP::Hep3Vector  _HPGeOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(-_STM_SSCoffset_Spot + offset_HPGe, 0., _HPGeZ_HPGe);
          stm._pSTMHPGeDetectorParams = std::unique_ptr<HPGeDetector>
          (new HPGeDetector(_HPGeBuild,
                            _HPGecrystalMaterial,
                            _HPGeholeMaterial,
                            _HPGewindowMaterial,
                            _HPGewallMaterial,
                            _HPGecapsuleMaterial, 
                            _HPGeEndcupR,
                            _HPGeEndcupL,
                            _HPGeCrystalR,
                            _HPGeCrystalL,
                            _HPGeZ_HPGe,
                            _HPGeHoleR,
                            _HPGeHoleL,
                            _HPGeCapsule_Wallthick,
                            _HPGeCapsule_Windowthick,
                            _HPGeCapsule_Endthick,
                            _HPGeCapsule_Walllength,
                            _HPGeWindowD,
                            _HPGeEndcupD,
                            _HPGeAirD,
                            _HPGeoffset_HPGe,
                            _HPGeOffsetInMu2e,
                            _HPGeRotation));


    const double offset_LaBr=_LaBroffset_LaBr;
    const CLHEP::HepRotation _LaBrRotation     = CLHEP::HepRotation::IDENTITY;
    const CLHEP::Hep3Vector  _LaBrOffsetInMu2e = _STMShieldingRef + CLHEP::Hep3Vector(_STM_SSCoffset_Spot + offset_LaBr, 0., _LaBrZ_LaBr);
          stm._pSTMLaBrDetectorParams = std::unique_ptr<LaBrDetector>
          (new LaBrDetector(_LaBrBuild,
                            _LaBrcrystalMaterial,
                            _LaBrwindowMaterial,
                            _LaBrwallMaterial,
                            _LaBrEndcupR,
                            _LaBrEndcupL,
                            _LaBrCrystalR,
                            _LaBrCrystalL,
                            _LaBrZ_LaBr,
                            _LaBrWindowD,
                            _LaBrEndcupD,
                            _LaBrAirD,
                            _LaBroffset_LaBr,
                            _LaBrOffsetInMu2e,
                            _LaBrRotation));


          stm._pSTMFrontShieldingParams = std::unique_ptr<FrontShielding>
          (new FrontShielding(_FrontShieldingBuild));

          stm._pSTMLeftShieldingParams = std::unique_ptr<LeftShielding>
          (new LeftShielding(_LeftShieldingBuild));

          stm._pSTMRightShieldingParams = std::unique_ptr<RightShielding>
          (new RightShielding(_RightShieldingBuild));

          stm._pSTMTopShieldingParams = std::unique_ptr<TopShielding>
          (new TopShielding(_TopShieldingBuild));

          stm._pSTMBottomShieldingParams = std::unique_ptr<BottomShielding>
          (new BottomShielding(_BottomShieldingBuild));

          stm._pSTMInnerShieldingParams = std::unique_ptr<InnerShielding>
          (new InnerShielding(_InnerShieldingBuild));

          stm._pSTMBackShieldingParams = std::unique_ptr<BackShielding>
          (new BackShielding(_BackShieldingBuild));

          stm._pSTMElectronicShieldingParams = std::unique_ptr<ElectronicShielding>
          (new ElectronicShielding(_ElectronicShieldingBuild,   
                                   _ElectronicShieldingSiGridX,                            
                                   _ElectronicShieldingSiGridY,
                                   _ElectronicShieldingSiGridZ,
                                   _ElectronicShieldingSiXcenter,
                                   _ElectronicShieldingSiYcenter,
                                   _ElectronicShieldingSiZcenter,
                                   _ElectronicShieldingConcreteT));

          stm._pSTMSTM_AbsorberParams = std::unique_ptr<STM_Absorber>
          (new STM_Absorber(_STM_AbsorberBuild,
                            _STM_Absorber_hW,
                            _STM_Absorber_hH,
                            _STM_Absorber_hT,
                            _STM_Absorber_GaptoSSC));

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
    _magnetHoleXOffset         = _config.getDouble("stm.magnet.holeXOffset", 0.);
    _magnetHoleYOffset         = _config.getDouble("stm.magnet.holeYOffset", 0.);
    _magnetMaterial            = _config.getString("stm.magnet.material");
    _magnetHasLiner            = _config.getBool("stm.magnet.hasLiner", true);
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
    _shieldHasLiner             = _config.getBool(  "stm.shield.hasLiner", true /*true default for backwards compatibility*/);
    _shieldLinerWidth           = _config.getDouble("stm.shield.widthLiner");
    _shieldRadiusOut            = _config.getDouble("stm.shield.rOut");
    _shieldPipeHalfLength       = _config.getDouble("stm.shield.pipe.halfLength");
    _shieldMaterialLiner        = _config.getString("stm.shield.materialLiner");
    _shieldMaterial             = _config.getString("stm.shield.material");
    _shieldMatchPipeBlock       = _config.getBool  ("stm.shield.matchPipeBlock", false);
    _shieldUpStrSpace           = _config.getDouble("stm.shield.UpStrSpace");
    _shieldDnStrSpace           = _config.getDouble("stm.shield.DnStrSpace");
    _shieldDnStrWallHalfLength  = _config.getDouble("stm.shield.DnStrWall.halfLength");
    _shieldDnStrWallHoleRadius  = _config.getDouble("stm.shield.DnStrWall.holeRadius", -1.);
    _shieldDnStrWallHalfHeight  = _config.getDouble("stm.shield.DnStrWall.halfHeight", -1.);
    _shieldDnStrWallHalfWidth   = _config.getDouble("stm.shield.DnStrWall.halfWidth", -1.);
    _shieldDnStrWallGap         = _config.getDouble("stm.shield.DnStrWall.gap", 0.);
    _shieldUpStrWallGap         = _config.getDouble("stm.shield.UpStrWall.gap", 0.); //only if using pipe as origin
    _shieldDnStrWallMaterial    = _config.getString("stm.shield.DnStrWall.material", _shieldMaterial);

    _stmDnStrEnvBuild       = _config.getBool("stm.downstream.build");
    _stmDnStrEnvHalfLength  = _config.getDouble("stm.downstream.halfLength");
    _stmDnStrEnvHalfWidth   = _config.getDouble("stm.downstream.halfWidth");
    _stmDnStrEnvHalfHeight  = _config.getDouble("stm.downstream.halfHeight");
    _stmDnStrEnvMaterial    = _config.getString("stm.downstream.material");

    _STM_SSCBuild          = _config.getBool(  "stm.STM_SSC.build");
    _STM_SSCdelta_WlR      = _config.getDouble("stm.STM_SSC.delta_WlR");
    _STM_SSCdelta_WlL      = _config.getDouble("stm.STM_SSC.delta_WlL");
    _STM_SSCW_length       = _config.getDouble("stm.STM_SSC.W_length");
    _STM_SSCW_height       = _config.getDouble("stm.STM_SSC.W_height");
    _STM_SSCleak           = _config.getDouble("stm.STM_SSC.leak");
    _STM_SSCWdepth_f       = _config.getDouble("stm.STM_SSC.Wdepth_f");
    _STM_SSCWdepth_b       = _config.getDouble("stm.STM_SSC.Wdepth_b");
    _STM_SSCAperture_HPGe1 = _config.getDouble("stm.STM_SSC.Aperture_HPGe1");
    _STM_SSCAperture_HPGe2 = _config.getDouble("stm.STM_SSC.Aperture_HPGe2");
    _STM_SSCAperture_LaBr1 = _config.getDouble("stm.STM_SSC.Aperture_LaBr1");
    _STM_SSCAperture_LaBr2 = _config.getDouble("stm.STM_SSC.Aperture_LaBr2");
    _STM_SSCoffset_Spot    = _config.getDouble("stm.STM_SSC.offset_Spot");
    _STM_SSCMaterial       = _config.getString("stm.STM_SSC.material");

    _HPGeBuild                  = _config.getBool("stm.HPGe.build");
    _HPGecrystalMaterial        = _config.getString("stm.HPGe.crystalMaterial");
    _HPGeholeMaterial           = _config.getString("stm.HPGe.holeMaterial");
    _HPGewindowMaterial         = _config.getString("stm.HPGe.windowMaterial");
    _HPGewallMaterial           = _config.getString("stm.HPGe.wallMaterial");
    _HPGecapsuleMaterial        = _config.getString("stm.HPGe.capsuleMaterial");
    _HPGeEndcupR                = _config.getDouble("stm.HPGe.EndcupR");
    _HPGeEndcupL                = _config.getDouble("stm.HPGe.EndcupL");
    _HPGeCrystalR               = _config.getDouble("stm.HPGe.CrystalR");
    _HPGeCrystalL               = _config.getDouble("stm.HPGe.CrystalL");
    _HPGeZ_HPGe                 = _config.getDouble("stm.HPGe.Z_HPGe");
    _HPGeHoleR                  = _config.getDouble("stm.HPGe.HoleR");
    _HPGeHoleL                  = _config.getDouble("stm.HPGe.HoleL");
    _HPGeCapsule_Wallthick      = _config.getDouble("stm.HPGe.Capsule_Wallthick");
    _HPGeCapsule_Windowthick    = _config.getDouble("stm.HPGe.Capsule_Windowthick");
    _HPGeCapsule_Endthick       = _config.getDouble("stm.HPGe.Capsule_Endthick");
    _HPGeCapsule_Walllength     = _config.getDouble("stm.HPGe.Capsule_Walllength");
    _HPGeWindowD                = _config.getDouble("stm.HPGe.WindowD");
    _HPGeEndcupD                = _config.getDouble("stm.HPGe.EndcupD");
    _HPGeAirD                   = _config.getDouble("stm.HPGe.AirD");
    _HPGeoffset_HPGe            = _config.getDouble("stm.HPGe.offset_HPGe");

    _LaBrBuild                  = _config.getBool("stm.LaBr.build");
    _LaBrcrystalMaterial        = _config.getString("stm.LaBr.crystalMaterial");
    _LaBrwindowMaterial         = _config.getString("stm.LaBr.windowMaterial");
    _LaBrwallMaterial           = _config.getString("stm.LaBr.wallMaterial");
    _LaBrEndcupR                = _config.getDouble("stm.LaBr.EndcupR");
    _LaBrEndcupL                = _config.getDouble("stm.LaBr.EndcupL");
    _LaBrCrystalR               = _config.getDouble("stm.LaBr.CrystalR");
    _LaBrCrystalL               = _config.getDouble("stm.LaBr.CrystalL");
    _LaBrZ_LaBr                 = _config.getDouble("stm.LaBr.Z_LaBr");
    _LaBrWindowD                = _config.getDouble("stm.LaBr.WindowD");
    _LaBrEndcupD                = _config.getDouble("stm.LaBr.EndcupD");
    _LaBrAirD                   = _config.getDouble("stm.LaBr.AirD");
    _LaBroffset_LaBr            = _config.getDouble("stm.LaBr.offset_LaBr");

    _FrontShieldingBuild        = _config.getBool("stm.FrontShielding.build");
    _LeftShieldingBuild         = _config.getBool("stm.LeftShielding.build");
    _RightShieldingBuild        = _config.getBool("stm.RightShielding.build");
    _TopShieldingBuild          = _config.getBool("stm.TopShielding.build");
    _BottomShieldingBuild       = _config.getBool("stm.BottomShielding.build");
    _InnerShieldingBuild        = _config.getBool("stm.InnerShielding.build");
    _BackShieldingBuild         = _config.getBool("stm.BackShielding.build");

    _ElectronicShieldingBuild      = _config.getBool("stm.ElectronicShielding.build");
    _ElectronicShieldingSiGridX    = _config.getDouble("stm.ElectronicShielding.SiGridX");
    _ElectronicShieldingSiGridY    = _config.getDouble("stm.ElectronicShielding.SiGridY");
    _ElectronicShieldingSiGridZ    = _config.getDouble("stm.ElectronicShielding.SiGridZ");
    _ElectronicShieldingSiXcenter  = _config.getDouble("stm.ElectronicShielding.SiXcenter");
    _ElectronicShieldingSiYcenter  = _config.getDouble("stm.ElectronicShielding.SiYcenter");
    _ElectronicShieldingSiZcenter  = _config.getDouble("stm.ElectronicShielding.SiZcenter");
    _ElectronicShieldingConcreteT  = _config.getDouble("stm.ElectronicShielding.ConcreteT");


    _STM_AbsorberBuild          = _config.getBool("stm.STM_Absorber.build");
    _STM_Absorber_hW            = _config.getDouble("stm.STM_Absorber.hW");
    _STM_Absorber_hH            = _config.getDouble("stm.STM_Absorber.hH");
    _STM_Absorber_hT            = _config.getDouble("stm.STM_Absorber.hT");
    _STM_Absorber_GaptoSSC      = _config.getDouble("stm.STM_Absorber.GaptoSSC");

  }

} // namespace mu2e
