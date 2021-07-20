#ifndef TrackerGeom_StrawProperties_hh
#define TrackerGeom_StrawProperties_hh
namespace mu2e {
  struct StrawProperties {
    double _strawInnerRadius;
    double _strawOuterRadius;
    double _strawWallThickness;
    double _outerMetalThickness;
    double _innerMetal1Thickness;
    double _innerMetal2Thickness;
    double _wireRadius;
    double _wirePlateThickness;
    double strawInnerRadius() const{ return _strawInnerRadius; }
    double strawOuterRadius() const{ return _strawOuterRadius; }
    double strawWallThickness() const{ return _strawWallThickness; }
    double outerMetalThickness() const{ return _outerMetalThickness; }
    double innerMetal1Thickness() const{ return _innerMetal1Thickness; }
    double innerMetal2Thickness() const{ return _innerMetal2Thickness; }
    double wireRadius()           const { return _wireRadius; }
    double wirePlateThickness()   const { return _wirePlateThickness; }

  };
}
#endif
