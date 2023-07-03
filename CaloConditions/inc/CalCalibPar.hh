#ifndef CaloConditions_CalCalibPar_hh
#define CaloConditions_CalCalibPar_hh

//
// a container for Calo SiPM calibrations
//

namespace mu2e {

class CalCalibPar {
 public:
  CalCalibPar(float ADC2MeV, int algID,
              float timeOffset) :
      _ADC2MeV(ADC2MeV),
      _algID(algID),
      _timeOffset(timeOffset) {}
  float ADC2MeV() const { return _ADC2MeV; }
  int algID() const { return _algID; }
  float timeOffset() const { return _timeOffset; }

  float _ADC2MeV;
  int  _algID;
  float _timeOffset;
};

} // namespace mu2e

#endif
