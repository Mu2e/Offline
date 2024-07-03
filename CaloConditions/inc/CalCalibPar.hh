#ifndef CaloConditions_CalCalibPar_hh
#define CaloConditions_CalCalibPar_hh

//
// a container for Calo SiPM calibrations
// will provide both energy and time calibration for specific SiPM

namespace mu2e {

class CalCalibPar {
 public:
  CalCalibPar(float ADC2MeV, float timeOffset) :
      _ADC2MeV(ADC2MeV),
      _timeOffset(timeOffset) {}
  float ADC2MeV() const { return _ADC2MeV; }
  float timeOffset() const { return _timeOffset; }

  float _ADC2MeV;
  float _timeOffset;
};

} // namespace mu2e

#endif
