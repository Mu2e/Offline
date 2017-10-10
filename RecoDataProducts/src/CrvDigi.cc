#include "RecoDataProducts/inc/CrvDigi.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<CrvDigi::CrvSingleWaveform> &CrvDigi::GetSingleWaveforms(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _waveforms[SiPMNumber];
  }

  std::vector<CrvDigi::CrvSingleWaveform> &CrvDigi::GetSingleWaveforms(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _waveforms[SiPMNumber];
  }

  const std::vector<CrvDigi::CrvSingleWaveform> &CrvDigi::GetSingleWaveforms(int fiberNumber, int side) const
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _waveforms[SiPMNumber];
  }

  const std::vector<CrvDigi::CrvSingleWaveform> &CrvDigi::GetSingleWaveforms(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _waveforms[SiPMNumber];
  }

  double CrvDigi::GetDigitizationPrecision() const 
  {
    return _digitizationPrecision;
  }

  void CrvDigi::SetDigitizationPrecision(double digitizationPrecision) 
  {
    _digitizationPrecision=digitizationPrecision;
  }

  int CrvDigi::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvDigi::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }

}
