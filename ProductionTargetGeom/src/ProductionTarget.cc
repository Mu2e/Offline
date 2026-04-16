#include "Offline/ProductionTargetGeom/inc/ProductionTarget.hh"
#include <numeric>

namespace mu2e {
  ProductionTarget::ProductionTarget(std::string tier1TargetType, int version, double rOut,
                                     double halfLength, double rotX,
                                     double rotY, const CLHEP::Hep3Vector& position,
                                     int    nFins = 0,
                                     double finHt = 0, double finThick = 0,
                                     double hubDisU = 0, double hubDisD = 0,
                                     double hubAngU = 0, double hubAngD = 0,
                                     double huboverU = 0, double huboverD = 0)
    : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    , _prodTargetPosition(position)
    , _tier1TargetType(tier1TargetType)
    , _version(version)
    , _rOut(rOut)
    , _halfLength(halfLength)
    , _envelHalfLength(halfLength)
    , _finHeight(finHt)
    , _finThickness(finThick)
    , _hubDistUS(hubDisU)
    , _hubDistDS(hubDisD)
    , _hubAngleUS(hubAngU)
    , _hubAngleDS(hubAngD)
    , _hubOverhangUS(huboverU)
    , _hubOverhangDS(huboverD)
  {
    _protonBeamRotation.rotateX(rotX).rotateY(rotY);
    _protonBeamInverseRotation = _protonBeamRotation.inverse();
  }

  ProductionTarget::ProductionTarget(
                                     std::string haymanTargetType, int version
                                     ,double productionTargetMotherOuterRadius
                                     ,double productionTargetMotherHalfLength
                                     ,double rOut
                                     ,double halfHaymanLength
                                     ,double rotHaymanX
                                     ,double rotHaymanY
                                     ,double rotHaymanZ
                                     ,const CLHEP::Hep3Vector& haymanProdTargetPosition
                                     ,std::string targetCoreMaterial
                                     ,std::string targetFinMaterial
                                     ,std::string targetVacuumMaterial
                                     ,std::string supportRingMaterial
                                     ,std::string spokeMaterial
                                     ,int numberOfTargetSections
                                     ,std::vector<double> startingSectionThickness
                                     ,std::vector<int> numberOfSegmentsPerSection
                                     ,std::vector<double> thicknessOfSegmentPerSection
                                     ,std::vector<double> heightOfRectangularGapPerSection
                                     ,std::vector<double> thicknessOfGapPerSection
                                     ,int nHaymanFins
                                     ,std::vector<double> finAngles
                                     ,double haymanFinThickness
                                     ,double finOuterRadius
                                     ,double supportRingLength
                                     ,double supportRingInnerRadius
                                     ,double supportRingOuterRadius
                                     ,double supportRingCutoutThickness
                                     ,double supportRingCutoutLength
                                     )
  : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    ,_haymanTargetType(haymanTargetType)
    ,_version(version)
    ,_productionTargetMotherOuterRadius(productionTargetMotherOuterRadius)
    ,_productionTargetMotherHalfLength(productionTargetMotherHalfLength)
    ,_rOut(rOut)
    ,_halfHaymanLength(halfHaymanLength)
    ,_rotHaymanX(rotHaymanX)
    ,_rotHaymanY(rotHaymanY)
    ,_rotHaymanZ(rotHaymanZ)
    ,_haymanProdTargetPosition(haymanProdTargetPosition)
    ,_targetCoreMaterial(targetCoreMaterial)
    ,_targetFinMaterial(targetFinMaterial)
    ,_targetVacuumMaterial(targetVacuumMaterial)
    ,_supportRingMaterial(supportRingMaterial)
    ,_spokeMaterial(spokeMaterial)
    ,_numberOfTargetSections(numberOfTargetSections)
    ,_startingSectionThickness(startingSectionThickness)
    ,_numberOfSegmentsPerSection(numberOfSegmentsPerSection)
    ,_thicknessOfSegmentPerSection(thicknessOfSegmentPerSection)
    ,_heightOfRectangularGapPerSection(heightOfRectangularGapPerSection)
    ,_thicknessOfGapPerSection(thicknessOfGapPerSection)
    ,_nHaymanFins(nHaymanFins)
    ,_finAngles(finAngles)
    ,_haymanFinThickness(haymanFinThickness)
    ,_finOuterRadius(finOuterRadius)
    ,_supportRingLength(supportRingLength)
    ,_supportRingInnerRadius(supportRingInnerRadius)
    ,_supportRingOuterRadius(supportRingOuterRadius)
    ,_supportRingCutoutThickness(supportRingCutoutThickness)
    ,_supportRingCutoutLength(supportRingCutoutLength)
  {
    //having this duplicated is inelegant but when it comes time to rip out Tier 1 I think it will be easier
    _protonBeamRotation.rotateX(rotHaymanX).rotateY(rotHaymanY).rotateZ(rotHaymanZ);
    _protonBeamInverseRotation = _protonBeamRotation.inverse();
    _halfLength = _productionTargetMotherHalfLength;
    _prodTargetPosition = _haymanProdTargetPosition;
  }

  ProductionTarget::ProductionTarget(
                                     std::string stickmanTargetType
                                     ,int version
                                     ,const StickmanEnvelopeParams& envelopeParams
                                     ,const StickmanPlateParams& plateParams
                                     ,const StickmanRodParams& rodParams
                                     ,const StickmanSpacerParams& spacerParams
                                     ,const StickmanSupportRingParams& supportRingParams
                                     )
  : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    ,_version(version)
    ,_productionTargetMotherOuterRadius(envelopeParams.productionTargetMotherOuterRadius)
    ,_productionTargetMotherHalfLength(envelopeParams.productionTargetMotherHalfLength)
    ,_targetVacuumMaterial(envelopeParams.targetVacuumMaterial)
    ,_stickmanTargetType(stickmanTargetType)
    ,_halfStickmanLength(envelopeParams.halfStickmanLength)
    ,_rotStickmanX(envelopeParams.rotStickmanX)
    ,_rotStickmanY(envelopeParams.rotStickmanY)
    ,_rotStickmanZ(envelopeParams.rotStickmanZ)
    ,_stickmanProdTargetPosition(envelopeParams.stickmanProdTargetPosition)
    ,_numberOfPlates(plateParams.numberOfPlates)
    ,_plateMaterial(plateParams.plateMaterial)
    ,_plateROut(plateParams.plateROut)
    ,_nStickmanFins(plateParams.nStickmanFins)
    ,_plateFinAngles(plateParams.plateFinAngles)
    ,_plateFinOuterRadius(plateParams.plateFinOuterRadius)
    ,_plateFinWidth(plateParams.plateFinWidth)
    ,_plateCenterToLugCenter(plateParams.plateCenterToLugCenter)
    ,_plateLugInnerRadius(plateParams.plateLugInnerRadius)
    ,_plateLugOuterRadius(plateParams.plateLugOuterRadius)
    ,_plateThickness(plateParams.plateThickness)
    ,_plateLugThickness(plateParams.plateLugThickness)
    ,_rodMaterial(rodParams.rodMaterial)
    ,_rodRadius(rodParams.rodRadius)
    ,_spacerMaterial(spacerParams.spacerMaterial)
    ,_spacerHalfLength(spacerParams.spacerHalfLength)
    ,_spacerOuterRadius(spacerParams.spacerOuterRadius)
    ,_spacerInnerRadius(spacerParams.spacerInnerRadius)
    ,_stickmanSupportRingMaterial(supportRingParams.stickmanSupportRingMaterial)
    ,_stickmanSupportRingLength(supportRingParams.stickmanSupportRingLength)
    ,_stickmanSupportRingInnerRadius(supportRingParams.stickmanSupportRingInnerRadius)
    ,_stickmanSupportRingOuterRadius(supportRingParams.stickmanSupportRingOuterRadius)
    ,_supportRingLugOuterRadius(supportRingParams.supportRingLugOuterRadius)
    ,_supportRingCutoutOffset(supportRingParams.supportRingCutoutOffset)
  {
    // rod half length, actual rod is longer than this since it inserts into the end rings, but for the geometry reconstruction, that part will be treated as part of the end ring.
    _rodHalfLength = std::accumulate(_plateLugThickness.begin(), _plateLugThickness.end(), 0.0) / 2.0 + 2 * _spacerHalfLength;

    // rotations
    _protonBeamRotation.rotateX(envelopeParams.rotStickmanX).rotateY(envelopeParams.rotStickmanY).rotateZ(envelopeParams.rotStickmanZ);
    _protonBeamInverseRotation = _protonBeamRotation.inverse(); // passive rotation matrix
    _halfLength = _productionTargetMotherHalfLength;
    _prodTargetPosition = _stickmanProdTargetPosition;
  }

  void ProductionTarget::configureStickman(const StickmanConfigParams& configParams) {
    // Configure plate fillet parameters
    _addFilletToPlateCore = configParams.addFilletToPlateCore;
    _addFilletToPlateLug = configParams.addFilletToPlateLug;
    _plateFilletRadius = configParams.plateFilletRadius;

    // Configure support ring fillet parameters
    _addFilletToSupportRingLug = configParams.addFilletToSupportRingLug;
    _supportRingLugFilletRadius = configParams.supportRingLugFilletRadius;

    // Configure support ring cutout parameters
    _addCutoutToSupportRing = configParams.addCutoutToSupportRing;
    _nSupportRingCutouts = configParams.nSupportRingCutouts;
    _supportRingCutoutAngles = configParams.supportRingCutoutAngles;
    _supportRingCutoutInnerRadius = configParams.supportRingCutoutInnerRadius;
    _supportRingCutoutTilt = configParams.supportRingCutoutTilt;

    // Configure support wheel parameters
    _supportsBuild = configParams.supportsBuild;
    _supportWheelRIn = configParams.supportWheelRIn;
    _supportWheelROut = configParams.supportWheelROut;
    _supportWheelHL = configParams.supportWheelHL;
    _supportWheelMaterial = configParams.supportWheelMaterial;
    _nSpokesPerSide = configParams.nSpokesPerSide;
    _supportWheelFeatureAngles = configParams.supportWheelFeatureAngles;
    _supportWheelFeatureArcs = configParams.supportWheelFeatureArcs;
    _supportWheelFeatureRIns = configParams.supportWheelFeatureRIns;
    _supportWheelRodHL = configParams.supportWheelRodHL;
    _supportWheelRodOffset = configParams.supportWheelRodOffset;
    _supportWheelRodPinOffset = configParams.supportWheelRodPinOffset;
    _supportWheelRodRadius = configParams.supportWheelRodRadius;
    _supportWheelRodRadialOffset = configParams.supportWheelRodRadialOffset;
    _supportWheelRodWireOffsetD = configParams.supportWheelRodWireOffsetD;
    _supportWheelRodWireOffsetU = configParams.supportWheelRodWireOffsetU;
    _supportWheelRodAngles = configParams.supportWheelRodAngles;
    _spokeTargetAnglesD = configParams.spokeTargetAnglesD;
    _spokeTargetAnglesU = configParams.spokeTargetAnglesU;
    _spokeRadius = configParams.spokeRadius;
    _spokeMaterial = configParams.spokeMaterial;
  }


}
