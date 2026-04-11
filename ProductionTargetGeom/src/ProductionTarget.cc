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
                                     ,double productionTargetMotherOuterRadius
                                     ,double productionTargetMotherHalfLength
                                     ,double rotStickmanX
                                     ,double rotStickmanY
                                     ,double rotStickmanZ
                                     ,double halfStickmanLength
                                     ,const CLHEP::Hep3Vector& stickmanProdTargetPosition
                                     ,std::string targetVacuumMaterial
                                     ,int numberOfPlates
                                     ,std::vector<std::string> plateMaterial
                                     ,std::vector<double> plateROut
                                     ,int nStickmanFins
                                     ,std::vector<double> plateFinAngles
                                     ,double plateFinOuterRadius
                                     ,double plateFinWidth
                                     ,double plateCenterToLugCenter
                                     ,double plateLugInnerRadius
                                     ,double plateLugOuterRadius
                                     ,std::vector<double> plateThickness
                                     ,std::vector<double> plateLugThickness
                                     ,std::string rodMaterial
                                     ,double rodRadius
                                     ,std::string spacerMaterial
                                     ,double spacerHalfLength
                                     ,double spacerOuterRadius
                                     ,double spacerInnerRadius
                                     ,std::string stickmanSupportRingMaterial
                                     ,double stickmanSupportRingLength
                                     ,double stickmanSupportRingInnerRadius
                                     ,double stickmanSupportRingOuterRadius
                                     ,double supportRingLugOuterRadius
                                     ,double supportRingCutoutOffset
                                     )
  : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    ,_version(version)
    ,_productionTargetMotherOuterRadius(productionTargetMotherOuterRadius)
    ,_productionTargetMotherHalfLength(productionTargetMotherHalfLength)
    ,_targetVacuumMaterial(targetVacuumMaterial)
    ,_stickmanTargetType(stickmanTargetType)
    ,_halfStickmanLength(halfStickmanLength)
    ,_rotStickmanX(rotStickmanX)
    ,_rotStickmanY(rotStickmanY)
    ,_rotStickmanZ(rotStickmanZ)
    ,_stickmanProdTargetPosition(stickmanProdTargetPosition)
    ,_numberOfPlates(numberOfPlates)
    ,_plateMaterial(plateMaterial)
    ,_plateROut(plateROut)
    ,_nStickmanFins(nStickmanFins)
    ,_plateFinAngles(plateFinAngles)
    ,_plateFinOuterRadius(plateFinOuterRadius)
    ,_plateFinWidth(plateFinWidth)
    ,_plateCenterToLugCenter(plateCenterToLugCenter)
    ,_plateLugInnerRadius(plateLugInnerRadius)
    ,_plateLugOuterRadius(plateLugOuterRadius)
    ,_plateThickness(plateThickness)
    ,_plateLugThickness(plateLugThickness)
    ,_rodMaterial(rodMaterial)
    ,_rodRadius(rodRadius)
    ,_spacerMaterial(spacerMaterial)
    ,_spacerHalfLength(spacerHalfLength)
    ,_spacerOuterRadius(spacerOuterRadius)
    ,_spacerInnerRadius(spacerInnerRadius)
    ,_stickmanSupportRingMaterial(stickmanSupportRingMaterial)
    ,_stickmanSupportRingLength(stickmanSupportRingLength)
    ,_stickmanSupportRingInnerRadius(stickmanSupportRingInnerRadius)
    ,_stickmanSupportRingOuterRadius(stickmanSupportRingOuterRadius)
    ,_supportRingLugOuterRadius(supportRingLugOuterRadius)
    ,_supportRingCutoutOffset(supportRingCutoutOffset)
  {
    // rod half length, actual rod is longer than this since it inserts into the end rings, but for the geometry reconstruction, that part will be treated as part of the end ring.
    _rodHalfLength = std::accumulate(_plateLugThickness.begin(), _plateLugThickness.end(), 0.0) / 2.0 + 2 * _spacerHalfLength;

    // rotations
    _protonBeamRotation.rotateX(rotStickmanX).rotateY(rotStickmanY).rotateZ(rotStickmanZ);
    _protonBeamInverseRotation = _protonBeamRotation.inverse(); // passive rotation matrix
    _halfLength = _productionTargetMotherHalfLength;
    _prodTargetPosition = _stickmanProdTargetPosition;
  }


}
