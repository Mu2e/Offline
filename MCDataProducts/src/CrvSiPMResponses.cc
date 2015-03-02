#include "MCDataProducts/inc/CrvSiPMResponses.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<mu2e::CrvSiPMResponses::CrvSingleSiPMResponse> &CrvSiPMResponses::GetSiPMResponses(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvSiPMResponses[SiPMNumber];
  }

  std::vector<mu2e::CrvSiPMResponses::CrvSingleSiPMResponse> &CrvSiPMResponses::GetSiPMResponses(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvSiPMResponses[SiPMNumber];
  }

  const std::vector<mu2e::CrvSiPMResponses::CrvSingleSiPMResponse> &CrvSiPMResponses::GetSiPMResponses(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvSiPMResponses[SiPMNumber];
  }

  const std::vector<mu2e::CrvSiPMResponses::CrvSingleSiPMResponse> &CrvSiPMResponses::GetSiPMResponses(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvSiPMResponses[SiPMNumber];
  }

  unsigned int CrvSiPMResponses::GetNumberOfSiPMResponses(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvSiPMResponses[SiPMNumber].size();
  }

  unsigned int CrvSiPMResponses::GetNumberOfSiPMResponses(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvSiPMResponses[SiPMNumber].size();
  }

  double CrvSiPMResponses::GetFirstSiPMResponseTime() const
  {
    double firstTime = NAN;
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      if(_crvSiPMResponses[SiPM].size()==0) continue;
      double t = _crvSiPMResponses[SiPM].front()._time;
      if(isnan(firstTime) || t<firstTime) firstTime=t;
    }
    return firstTime;
  }

  double CrvSiPMResponses::GetLastSiPMResponseTime() const
  {
    double lastTime = NAN;
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      if(_crvSiPMResponses[SiPM].size()==0) continue;
      double t = _crvSiPMResponses[SiPM].back()._time;
      if(isnan(lastTime) || t<lastTime) lastTime=t;
    }
    return lastTime;
  }

  int CrvSiPMResponses::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvSiPMResponses::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }
}

