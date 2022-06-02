#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/GeomPrimitives/inc/Box.hh"

namespace mu2e {

  // nameSuffix lets you give multiple PWC's unique names.
  // wireNumStart is the first number to be used when naming the sections of
  // gas corresponding to wires. This same wire numbering is used for the
  // copyNo argument when placing the wire gas in the geometry.
  PTMPWC::PTMPWC(std::string const& nameSuffix,
                     double frameHeight,
                     double frameWidth,
                     double frameThick,
                     double outerPlateThick,
                     std::string const& frameMaterialName,
                     double windowHeight,
                     double windowWidth,
                     double windowThick,
                     std::string const& windowMaterialName,
                     std::string const& gasMaterialName,
                     int numVertWires,
                     int numHorizWires,
                     CLHEP::Hep3Vector const& originInParent,
                     int wireNumStart,
                     double containerMargin,
                     int framesInDetector,
                     int outerPlatesInDetector,
                     double ground1Zframes,
                     double hv1Zframes,
                     double hv2Zframes,
                     double hv3Zframes,
                     double ground2Zframes) :
    _nameSuffix(nameSuffix),
    _originInParent(originInParent),
    _frameHeight(frameHeight),
    _frameWidth(frameWidth),
    _frameThick(frameThick),
    _outerPlateThick(outerPlateThick),
    _frameMaterialName(frameMaterialName),
    _windowMaterialName(windowMaterialName),
    _gasMaterialName(gasMaterialName),
    _wireNumStart(wireNumStart),
    _numVertWires(numVertWires),
    _numHorizWires(numHorizWires)
  {
    _detectorThick = (framesInDetector*frameThick) + (outerPlatesInDetector*outerPlateThick);
    _totalThick = _detectorThick + containerMargin;
    _totalHeight = frameHeight + containerMargin;
    _totalWidth = frameWidth + containerMargin;
    _ground1Z = ground1Zframes*frameThick;
    _hv1Z = hv1Zframes*frameThick;
    _hv2Z = hv2Zframes*frameThick;
    _hv3Z = hv3Zframes*frameThick;
    _ground2Z = ground2Zframes*frameThick;
    _gasInZ = 0.5*(_hv1Z + _ground1Z);
    _vertWireZ = 0.5*(_hv2Z + _hv1Z);
    _horizWireZ = 0.5*(_hv3Z + _hv2Z);
    _gasOutZ = 0.5*(_ground2Z + _hv3Z);

    _pwcWindow.reset(new Box(0.5*windowWidth, 0.5*windowHeight, 0.5*windowThick));
    double vertWireHalfHeight = 0.5*windowHeight / numVertWires;
    double vertWireHalfThick = 0.5*(_hv2Z - _hv1Z - windowThick);
    _vertWireGasSection.reset(new Box(0.5*windowWidth, vertWireHalfHeight, vertWireHalfThick));
    double horizWireHalfWidth = 0.5*windowWidth / numHorizWires;
    double horizWireHalfThick = 0.5*(_hv3Z - _hv2Z - windowThick);
    _horizWireGasSection.reset(new Box(horizWireHalfWidth, 0.5*windowHeight, horizWireHalfThick));
    double gasInHalfThick = 0.5*(_hv1Z - _ground1Z - windowThick);
    _gasSection1.reset(new Box(0.5*windowWidth, 0.5*windowHeight, gasInHalfThick));
    double gasOutHalfThick = 0.5*(_ground2Z - _hv3Z - windowThick);
    _gasSection4.reset(new Box(0.5*windowWidth, 0.5*windowHeight, gasOutHalfThick));

    for (int i=0; i<numVertWires; i++) {
      double gasY2 = (-0.5*windowHeight) + ((i+0.5)*2.*vertWireHalfHeight);
      _vertWireYpos.push_back(gasY2);
    }
    for (int i=0; i<numHorizWires; i++) {
      double gasX3 = (0.5*windowWidth) - ((i+0.5)*2.*horizWireHalfWidth);
      _horizWireXpos.push_back(gasX3);
    }
    _upstreamWindowSurfaceZ = _ground1Z - (0.5*windowThick);
  } // ctor

} // namespace mu2e
