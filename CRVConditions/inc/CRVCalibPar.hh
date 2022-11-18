#ifndef CRVConditions_CRVCalibPar_hh
#define CRVConditions_CRVCalibPar_hh

//
// a little container for CRV SiPM calibrations

namespace mu2e {

class CRVCalibPar {
 public:
  CRVCalibPar(float pedestal, float height, float area, float time) :
      _pedestal(pedestal), _height(height), _area(area), _time(time) {}
  float pedestal() const { return _pedestal; }
  float height() const { return _height; }
  float area() const { return _area; }
  float time() const { return _time; }

  float _pedestal;
  float _height;
  float _area;
  float _time;
};

}  // namespace mu2e

#endif
