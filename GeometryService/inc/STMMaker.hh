#ifndef STMGeom_STMMaker_hh
#define STMGeom_STMMaker_hh
//
// Class to construct and return STM
//
// Author: Anthony Palladino
//

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

// Held by value below while the hand-stacked sections are assembled,
// so these need the definitions rather than a forward declaration.
#include "Offline/STMGeom/inc/BrickWall.hh"
#include "Offline/STMGeom/inc/FrontShieldingRight.hh"
#include "Offline/STMGeom/inc/FrontShieldingLeft.hh"
#include "Offline/STMGeom/inc/RightWall.hh"
#include "Offline/STMGeom/inc/LeftWall.hh"
#include "Offline/STMGeom/inc/TopWall.hh"
#include "Offline/STMGeom/inc/BottomWall.hh"
#include "Offline/STMGeom/inc/InnerShielding.hh"

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

    // updated, hand-stacked downstream shielding
    bool         _handstacked;

    // The standard lead bricks, shared by every structure that stacks them
    // Numbered by type, so a new size is a new entry rather than a new
    // member and a change to every caller.
    std::vector<std::string>       _leadBrickNames;
    std::vector<CLHEP::Hep3Vector> _leadBrickDims;
    double      _leadBrickWearY;
    double      _leadBrickWearXZ;
    std::string _leadBrickMaterial;

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
    // updated (hand-stacked) SSC: one block per slab, bores stepped in z
    double      _STM_SSCW_width;
    double      _STM_SSCr_LaBr_f;
    double      _STM_SSCr_HPGe_f;
    double      _STM_SSCr_LaBr_b;
    double      _STM_SSCr_HPGe_b;
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
    // updated (hand-stacked) cradle: five plates around the collimator
    double      _SSCSupportdepth;
    double      _SSCSupportside_T;
    double      _SSCSupportside_H;
    double      _SSCSupportplate_base_T;
    double      _SSCSupportbottom_T;
    double      _SSCSupporttop_T;
    double      _STM_SSCboreToBase;
    std::string _SSCSupportMaterial;

    // SSC front shield: the brick wall, the shelf and two poly blocks.
    // Centers are stored with every offset already applied.
    bool        _SSCFrontShieldBuild;
    std::vector<BrickWallBrick>    _SSCFrontShieldBricks;
    std::string _SSCFrontShieldShelfMaterial;
    CLHEP::Hep3Vector _SSCFrontShieldShelfDim;
    CLHEP::Hep3Vector _SSCFrontShieldShelfCenter;
    std::string _SSCFrontShieldPoly1Material;
    CLHEP::Hep3Vector _SSCFrontShieldPoly1Dim;
    CLHEP::Hep3Vector _SSCFrontShieldPoly1Center;
    double      _SSCFrontShieldPoly1BoreR;
    double      _SSCFrontShieldPoly1BoreDX;
    double      _SSCFrontShieldPoly1BoreDY;
    std::string _SSCFrontShieldPoly2Material;
    CLHEP::Hep3Vector _SSCFrontShieldPoly2Dim;
    CLHEP::Hep3Vector _SSCFrontShieldPoly2Center;

    // Front shielding, right. The lead courses and the poly sheets are
    // expanded from the config's type lists in parseConfig, so what is
    // held here is already the finished placement.
    bool        _FrontShieldingRightBuild;
    std::vector<BrickWall>                _FrontShieldingRightLeadLayers;
    std::vector<FrontShieldingRightSheet> _FrontShieldingRightSheets;
    // The four blocks are identical, so their material and size are
    // shared and only their centers vary.
    std::string _FrontShieldingRightBlockMaterial;
    CLHEP::Hep3Vector _FrontShieldingRightBlockHalfDim;
    std::vector<CLHEP::Hep3Vector>        _FrontShieldingRightBlockCenter;
    std::vector<FrontShieldingRightPipe>  _FrontShieldingRightPipes;
    FrontShieldingRightPlate              _FrontShieldingRightPlate;
    double      _FrontShieldingRightBackZ;

    // Landmarks of the front shielding that later sections butt
    // against. Resolved where the section is expanded and kept so a
    // wall downstream need not re-derive them from the layer sequence.
    double      _FrontShieldingRightBrickEndX;
    double      _FrontShieldingRightFrontZ;
    double      _FrontShieldingRightLastLeadBackZ;

    // The shield house's reference corner, where the walls meet the
    // front shielding. Resolved once before any wall is expanded, so
    // the right, left and top walls all key off the same point rather
    // than each re-deriving it -- and so none of them has to read
    // another wall's configuration keys.
    double      _ShieldHouseRefX;
    double      _ShieldHouseRefZ;

    // Front shielding, left. Nothing here is bored and nothing
    // interleaves, so it is three brick groups, four sheets and the
    // prism, all keyed off the right half.
    bool        _FrontShieldingLeftBuild;
    std::vector<FrontShieldingLeftSheet>  _FrontShieldingLeftSheets;
    std::vector<BrickWall>                _FrontShieldingLeftBrickGroups;
    FrontShieldingLeftPrism               _FrontShieldingLeftPrism;

    // Right wall. Five layers worked inward from the reference vertex,
    // expanded in parseConfig so what is held here is already placed.
    bool        _RightWallBuild;
    std::vector<BrickWall> _RightWallLeadLayers;
    std::vector<WallSheet> _RightWallSheets;
    RightWallEdgePrism     _RightWallEdgePrism;

    // Left wall. The same five kinds of layer, without the edge prism.
    bool        _LeftWallBuild;
    std::vector<BrickWall> _LeftWallLeadLayers;
    std::vector<WallSheet> _LeftWallSheets;

    // Top wall. One lead layer laid flat between two pairs of plates.
    bool        _TopWallBuild;
    std::vector<BrickWall> _TopWallLeadLayers;
    std::vector<WallSheet> _TopWallSheets;

    // Bottom wall. The baseplate plus two lead layers and four swept
    // prisms; each lead layer holds its courses and its strays.
    bool        _BottomWallBuild;
    BottomWallPlate               _BottomWallPlate;
    std::vector<BrickWall>        _BottomWallLeadLayers;
    std::vector<BottomWallPrism>  _BottomWallPrisms;

    // Inner shielding, updated. Not a stack of layers: every piece is
    // placed on its own, so it is three lists rather than a sequence.
    std::vector<InnerShieldingPrism>  _InnerShieldingPrisms;
    std::vector<InnerShieldingBox>    _InnerShieldingBoxes;
    std::vector<BrickWallBrick>       _InnerShieldingBricks;
    std::vector<BrickWallBore>        _InnerShieldingBores;

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
    double  _FrontSHole_r;
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
    double  _BackSPipeGap;
    std::string _BackSBPMaterial = "BP";  // read only in the hand-stacked geometry

    bool   _InnerShieldingBuild;

    bool     _ElectronicShieldingBuild;
    double   _ElectronicSSiGridX;
    double   _ElectronicSSiGridY;
    double   _ElectronicSSiGridZ;
    double   _ElectronicSSiXcenter;
    double   _ElectronicSSiYcenter;
    double   _ElectronicSSiZcenter;
    double   _ElectronicSConcreteT = 0.;
    double   _ElectronicSGapToSi = 0.;
    // hand-stacked only
    double   _ElectronicSXOffset = 0.;
    CLHEP::Hep3Vector _ElectronicSConcreteDims;
    int      _ElectronicSConcreteBlockN = 0;
    double   _ElectronicSConcreteBlockGap = 0.;
    double   _ElectronicSConcreteFromGround = 0.;
    double   _ElectronicSConcreteToWall = 0.;
    int      _ElectronicSSiNX = 0;
    int      _ElectronicSSiNY = 0;
    double   _ElectronicSSiFromGround = 0.;
    std::vector<double> _ElectronicSSiGapsToConcrete;


    bool     _STM_AbsorberBuild;
    double   _STM_Absorber_hW;
    double   _STM_Absorber_hH;
    double   _STM_Absorber_hT;
    double   _STM_Absorber_GaptoSSC;

  };

}  //namespace mu2e

#endif /* STMGeom_STMMaker_hh */
