#ifndef CRVConditions_CRVCalibPar_hh
#define CRVConditions_CRVCalibPar_hh

//
// a little container for CRV SiPM calibrations

namespace mu2e {

class CRVCalibPar {
 public:
  CRVCalibPar(float pedestal, float height, float area, float timeOffset) :
      _pedestal(pedestal), _height(height), _area(area),
      _timeOffset(timeOffset) {}
  float pedestal() const { return _pedestal; }
  float height() const { return _height; }
  float area() const { return _area; }
  float timeOffset() const { return _timeOffset; }

  float _pedestal;
  float _height;
  float _area;
  float _timeOffset;
};

}  // namespace mu2e

#endif
