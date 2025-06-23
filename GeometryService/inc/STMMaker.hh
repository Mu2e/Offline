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
    ~STMMaker() = default;

    // delete automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    STMMaker( STMMaker const & ) = delete;
    STMMaker( STMMaker&&       ) = delete;
    STMMaker& operator=( STMMaker const & ) = delete;
    STMMaker& operator=( STMMaker&&       ) = delete;

    void parseConfig( SimpleConfig const & _config );

    // This is deprecated and will go away soon.
    // Still needed for root graphics version.
    const STM& getSTM() const { return *_stm;}

    // This is the accessor that will remain.
    std::unique_ptr<STM> getSTMPtr() { return std::move(_stm); }

  private:

    std::unique_ptr<STM> _stm;

    int          _verbosityLevel;
    double       _stmZAllowed;

    double       _stmReferenceZ;

    bool         _stmDnStrEnvBuild;
    double       _stmDnStrEnvHalfLength;
    double       _stmDnStrEnvHalfWidth;
    double       _stmDnStrEnvHalfHeight;
    std::string  _stmDnStrEnvMaterial;

    bool         _magnetBuild;
    double       _magnetUpStrSpace;
    double       _magnetHalfLength;
    double       _magnetHalfWidth;
    double       _magnetHalfHeight;
    double       _magnetHoleHalfLength;
    double       _magnetHoleHalfWidth;
    double       _magnetHoleHalfHeight;
    double       _magnetHoleXOffset;
    double       _magnetHoleYOffset;
    std::string  _magnetMaterial;
    bool         _magnetHasLiner;
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
    bool         _shieldHasLiner;
    double       _shieldLinerWidth;
    double       _shieldRadiusOut;
    double       _shieldPipeHalfLength;
    std::string  _shieldMaterialLiner;
    std::string  _shieldMaterial;
    bool         _shieldMatchPipeBlock;
    double       _shieldUpStrSpace;
    double       _shieldDnStrSpace;
    double       _shieldDnStrWallHalfLength;
    double       _shieldDnStrWallHoleRadius;
    double       _shieldDnStrWallHalfHeight;
    double       _shieldDnStrWallHalfWidth;
    double       _shieldDnStrWallGap;
    double       _shieldUpStrWallGap;
    double       _shieldPipeUpStrAirGap;
    std::string  _shieldDnStrWallMaterial;
    bool         _shieldBuildMatingBlock;

    bool        _STM_SSCBuild;
    bool        _STM_SSCVDBuild;
    double      _STM_SSCdelta_WlR;
    double      _STM_SSCdelta_WlL;
    double      _STM_SSCW_middle;
    double      _STM_SSCW_height;
    double      _STM_SSCWdepth_f;
    double      _STM_SSCWdepth_b;
    double      _STM_SSCAperture_HPGe1;
    double      _STM_SSCAperture_HPGe2;
    double      _STM_SSCAperture_LaBr1;
    double      _STM_SSCAperture_LaBr2;
    double      _STM_SSCoffset_Spot;
    double      _STM_SSCleak;
    double      _STM_SSCFrontToWall;
    double      _STM_SSCZGap;
    double      _STM_SSCZGapBack;
    std::string _STM_SSCMaterial;


    bool        _SSCSupportBuild;
    double      _SSCSupporttable_L;
    double      _SSCSupporttable_H;
    double      _SSCSupporttable_T;
    double      _SSCSupportleg_L;
    double      _SSCSupportleg_H;
    double      _SSCSupportleg_T;
    double      _SSCSupportbase_L;
    double      _SSCSupportbase_H;
    double      _SSCSupportbase_T;
    double      _SSCSupportwall_L;
    double      _SSCSupportwall_H;
    double      _SSCSupportwall_T;
    double      _SSCSupporthole_H;
    double      _SSCSupporthole_T;
    double      _SSCSupportFLeadStand_L;
    double      _SSCSupportFLeadStand_H;
    double      _SSCSupportFLeadStand_T;
    double      _SSCSupportFLeadShim_H;
    double      _SSCSupportFLeadShim_T;
    double      _SSCSupportFAluminumShim_T;
    double      _SSCSupportFAluminumExtra_L;
    double      _SSCSupportFAluminumExtra_H;


    bool    _FrontShieldingBuild;
    double  _FrontSHeightofRoom;
    double  _FrontStungstenlength;
    double  _FrontStungstendepth;
    double  _FrontSleaddepth1;
    double  _FrontSleaddepth2;
    double  _FrontSaluminumdepth;
    double  _FrontScopperdepth;
    double  _FrontSBPdepth;
    double  _FrontSfPb_lengthL;
    double  _FrontSfPb_lengthR;
    double  _FrontSGapForTop;
    double  _FrontSLeakForSSC;
    double  _FrontSCopperL;
    double  _FrontS_H;
    double  _FrontS_Thickness;
    double  _FrontS_Length;


    bool   _HPGeBuild;
    std::string _HPGecrystalMaterial;
    std::string _HPGeholeMaterial;
    std::string _HPGewindowMaterial;
    std::string _HPGewallMaterial;
    std::string _HPGecapsuleMaterial;
    double _HPGeEndcapR;
    double _HPGeEndcapL;
    double _HPGeCrystalR;
    double _HPGeCrystalL;
    double _HPGeZ_HPGe;
    double _HPGeHoleR;
    double _HPGeHoleL;
    double _HPGeCapsule_Wallthick;
    double _HPGeCapsule_Windowthick;
    double _HPGeCapsule_Endthick;
    double _HPGeCapsule_Walllength;
    double _HPGeWindowD;
    double _HPGeEndcapD;
    double _HPGeAirD;
    double _HPGeoffset_HPGe;

    bool   _LaBrBuild;
    std::string _LaBrcrystalMaterial;
    std::string _LaBrwindowMaterial;
    std::string _LaBrwallMaterial;
    double _LaBrEndcapR;
    double _LaBrEndcapL;
    double _LaBrCrystalR;
    double _LaBrCrystalL;
    double _LaBrZ_LaBr;
    double _LaBrWindowD;
    double _LaBrEndcapD;
    double _LaBrAirD;
    double _LaBroffset_LaBr;

    bool    _BottomShieldingBuild;
    double  _BottomSFront_LB;
    double  _BottomSFront_LB_inner;
    double  _BottomSfloor_Zlength;
    double  _BottomSleaddepth;
    double  _BottomScopperdepth;
    double  _BottomSBPdepth;


    bool    _LeftShieldingBuild;
    double  _LeftS_Length;
    double  _LeftSleaddepth;
    double  _LeftScopperdepth;
    double  _LeftSBPdepth;
    double  _LeftSXmin;

    bool    _RightShieldingBuild;
    double  _RightS_Length;
    double  _RightSleaddepth;
    double  _RightScopperdepth;
    double  _RightSBPdepth;
    double  _RightSXmax;


    bool    _TopShieldingBuild;
    bool    _TopShieldingSkirtBuild;
    double  _TopLiftBeam_L;
    double  _TopLiftBeam_H;
    double  _TopLiftBeam_T;
    double  _TopLiftBeam_Xmove;
    double  _TopSZlength;
    double  _TopSXlength;
    double  _TopSFront_LT;
    double  _TopTFZlength;
    double  _TopTFXlength;
    double  _TopTBZlength;
    double  _TopScontainerdepth;
    double  _TopSleaddepth;
    double  _TopScopperdepth;
    double  _TopSBPdepth;
    double  _TopSZHole;
    double  _TopSBarLeft;
    double  _TopSBarRight;
    double  _TopSGapLeft;
    double  _TopSGapRight;
    double  _TopSLeak;

    bool    _BackShieldingBuild;
    double  _BackSBPThick;
    double  _BackSBPLength;
    double  _BackSBPHeight;
    double  _BackS_dX;
    double  _BackS_dY;


    bool   _InnerShieldingBuild;

    bool     _ElectronicShieldingBuild;
    double   _ElectronicSSiGridX;
    double   _ElectronicSSiGridY;
    double   _ElectronicSSiGridZ;
    double   _ElectronicSSiXcenter;
    double   _ElectronicSSiYcenter;
    double   _ElectronicSSiZcenter;
    double   _ElectronicSConcreteT;
    double   _ElectronicSGapToSi;


    bool     _STM_AbsorberBuild;
    double   _STM_Absorber_hW;
    double   _STM_Absorber_hH;
    double   _STM_Absorber_hT;
    double   _STM_Absorber_GaptoSSC;

  };

}  //namespace mu2e

#endif /* STMGeom_STMMaker_hh */
