#include "MCDataProducts/inc/CrvWaveforms.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<double> &CrvWaveforms::GetWaveform(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _waveform[SiPMNumber];
  }

  std::vector<double> &CrvWaveforms::GetWaveform(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _waveform[SiPMNumber];
  }

  const std::vector<double> &CrvWaveforms::GetWaveform(int fiberNumber, int side) const
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _waveform[SiPMNumber];
  }

  const std::vector<double> &CrvWaveforms::GetWaveform(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _waveform[SiPMNumber];
  }

  double CrvWaveforms::GetStartTime(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _startTime[SiPMNumber];
  }

  double CrvWaveforms::GetStartTime(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _startTime[SiPMNumber];
  }

  void CrvWaveforms::SetStartTime(int fiberNumber, int side, double startTime) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    _startTime[SiPMNumber]=startTime;
  }

  void CrvWaveforms::SetStartTime(int SiPMNumber, double startTime) 
  {
    CheckSiPMNumber(SiPMNumber);
    _startTime[SiPMNumber]=startTime;
  }

  double CrvWaveforms::GetBinWidth() const 
  {
    return _binWidth;
  }

  void CrvWaveforms::SetBinWidth(double binWidth) 
  {
    _binWidth=binWidth;
  }

  int CrvWaveforms::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvWaveforms::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }

}
