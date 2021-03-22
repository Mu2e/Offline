#include "ProductionTargetGeom/inc/ProductionTarget.hh"

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


}
