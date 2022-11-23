#ifndef CRVConditions_CRVCalibPar_hh
#define CRVConditions_CRVCalibPar_hh

//
// a little container for CRV SiPM calibrations
//

namespace mu2e {

class CRVCalibPar {
 public:
  CRVCalibPar(float pedestal, float pulseHeight, float pulseArea,
              float timeOffset) :
      _pedestal(pedestal),
      _pulseHeight(pulseHeight), _pulseArea(pulseArea),
      _timeOffset(timeOffset) {}
  float pedestal() const { return _pedestal; }
  float pulseHeight() const { return _pulseHeight; }
  float pulseArea() const { return _pulseArea; }
  float timeOffset() const { return _timeOffset; }

  float _pedestal;
  float _pulseHeight;
  float _pulseArea;
  float _timeOffset;
};

}  // namespace mu2e

#endif
