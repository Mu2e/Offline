#ifndef STMGeom_STMMaker_hh
#define STMGeom_STMMaker_hh
//
// Class to construct and return STM
//
// Author: Anthony Palladino
//

#include <memory>
#include <string>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class STM;
  class SimpleConfig;
  class Tube;

  class STMMaker {

  public:

    STMMaker( SimpleConfig const & config,
              double solenoidOffset);

    void parseConfig( SimpleConfig const & _config );

    // This is deprecated and will go away soon.
    // Still needed for root graphics version.
    const STM& getSTM() const { return *_stm;}

    // This is the accessor that will remain.
    std::unique_ptr<STM> getSTMPtr() { return std::move(_stm); }

  private:

    // hide automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    STMMaker( STMMaker const & );
    STMMaker const & operator= ( STMMaker const & );

    std::unique_ptr<STM> _stm;

    int          _verbosityLevel;
    double       _stmZAllowed;
    
    bool         _magnetBuild;
    double       _magnetUpStrSpace;
    double       _magnetHalfLength;
    double       _magnetHalfWidth;
    double       _magnetHalfHeight;
    double       _magnetHoleHalfLength;
    double       _magnetHoleHalfWidth;
    double       _magnetHoleHalfHeight;
    std::string  _magnetMaterial;
    double       _magnetField;
    bool         _magnetFieldVisible;
    
    bool        _FOVCollimatorBuild;
    std::string _FOVCollimatorMaterial;
    double      _FOVCollimatorUpStrSpace;
    double      _FOVCollimatorHalfWidth;
    double      _FOVCollimatorHalfHeight;
    double      _FOVCollimatorHalfLength;
    bool        _FOVCollimatorLinerBuild;
    std::string _FOVCollimatorLinerMaterial;
    double      _FOVCollimatorLinerHalfWidth;
    double      _FOVCollimatorLinerHalfHeight;
    double      _FOVCollimatorLinerHalfLength;
    double      _FOVCollimatorLinerCutOutHalfLength;
    double      _FOVCollimatorHole1xOffset;
    double      _FOVCollimatorHole1RadiusUpStr;
    double      _FOVCollimatorHole1RadiusDnStr;
    bool        _FOVCollimatorHole1LinerBuild;
    double      _FOVCollimatorHole1LinerThickness;
    bool        _FOVCollimatorHole2Build;
    double      _FOVCollimatorHole2xOffset;
    double      _FOVCollimatorHole2RadiusUpStr;
    double      _FOVCollimatorHole2RadiusDnStr;
    bool        _FOVCollimatorHole2LinerBuild;
    double      _FOVCollimatorHole2LinerThickness;
    std::string _FOVCollimatorHoleLinerMaterial;
    
    bool         _magnetTableBuild;
    std::string  _magnetTableMaterial;
    double       _magnetTableTopExtraWidth;
    double       _magnetTableTopExtraLength;
    double       _magnetTableTopHalfHeight;
    double       _magnetTableLegRadius;
    
    bool         _pipeBuild;
    double       _pipeRadiusIn;
    double       _pipeRadiusOut;
    std::string  _pipeMaterial;
    std::string  _pipeGasMaterial;
    double       _pipeUpStrSpace;
    double       _pipeDnStrHalfLength;
    std::string  _pipeUpStrWindowMaterial;
    double       _pipeUpStrWindowHalfLength;
    std::string  _pipeDnStrWindowMaterial;
    double       _pipeDnStrWindowHalfLength;
    double       _pipeFlangeHalfLength;
    double       _pipeFlangeOverhangR;

    bool        _SSCollimatorBuild;
    std::string _SSCollimatorMaterial;
    double      _SSCollimatorUpStrSpace;
    double      _SSCollimatorHalfWidth;
    double      _SSCollimatorHalfHeight;
    double      _SSCollimatorHalfLength;
    bool        _SSCollimatorLinerBuild;
    std::string _SSCollimatorLinerMaterial;
    double      _SSCollimatorLinerHalfWidth;
    double      _SSCollimatorLinerHalfHeight;
    double      _SSCollimatorLinerHalfLength;
    double      _SSCollimatorLinerCutOutHalfLength;
    double      _SSCollimatorHole1xOffset;
    double      _SSCollimatorHole1RadiusUpStr;
    double      _SSCollimatorHole1RadiusDnStr;
    bool        _SSCollimatorHole1LinerBuild;
    double      _SSCollimatorHole1LinerThickness;
    bool        _SSCollimatorHole2Build;
    double      _SSCollimatorHole2xOffset;
    double      _SSCollimatorHole2RadiusUpStr;
    double      _SSCollimatorHole2RadiusDnStr;    
    bool        _SSCollimatorHole2LinerBuild;
    double      _SSCollimatorHole2LinerThickness;
    std::string _SSCollimatorHoleLinerMaterial;
    
    bool         _detectorTableBuild;
    std::string  _detectorTableMaterial;
    double       _detectorTableTopExtraWidth;
    double       _detectorTableTopExtraLength;
    double       _detectorTableTopHalfHeight;
    double       _detectorTableLegRadius;
    
    bool         _detector1Build;
    std::string  _detector1CrystalMaterial;
    double       _detector1CrystalRadiusIn;
    double       _detector1CrystalRadiusOut;
    double       _detector1CrystalHalfLength;
    double       _detector1xOffset;
    std::string  _detector1CanMaterial;
    double       _detector1CanRadiusIn;
    double       _detector1CanRadiusOut;
    double       _detector1CanHalfLength;
    double       _detector1CanUpStrSpace;
    std::string  _detector1CanUpStrWindowMaterial;
    double       _detector1CanUpStrWindowHalfLength;
    std::string  _detector1CanGasMaterial;
    
    bool         _detector2Build;
    std::string  _detector2CrystalMaterial;
    double       _detector2CrystalRadiusIn;
    double       _detector2CrystalRadiusOut;
    double       _detector2CrystalHalfLength;
    double       _detector2xOffset;
    std::string  _detector2CanMaterial;
    double       _detector2CanRadiusIn;
    double       _detector2CanRadiusOut;
    double       _detector2CanHalfLength;
    double       _detector2CanUpStrSpace;
    std::string  _detector2CanUpStrWindowMaterial;
    double       _detector2CanUpStrWindowHalfLength;
    std::string  _detector2CanGasMaterial;
    
    bool         _shieldBuild;
    double       _shieldRadiusIn;
    double       _shieldLinerWidth;
    double       _shieldRadiusOut;
    double       _shieldPipeHalfLength;
    std::string  _shieldMaterialLiner;
    std::string  _shieldMaterial;
    double       _shieldUpStrSpace;
    double       _shieldDnStrSpace;
    double       _shieldDnStrWallHalfLength;
    
  };

}  //namespace mu2e

#endif /* STMGeom_STMMaker_hh */
