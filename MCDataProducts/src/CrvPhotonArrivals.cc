#include "MCDataProducts/inc/CrvPhotonArrivals.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<CrvPhotonArrivals::SinglePhoton> &CrvPhotonArrivals::GetPhotonArrivalTimes(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _times[SiPMNumber];
  }

  std::vector<CrvPhotonArrivals::SinglePhoton> &CrvPhotonArrivals::GetPhotonArrivalTimes(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _times[SiPMNumber];
  }

  const std::vector<CrvPhotonArrivals::SinglePhoton> &CrvPhotonArrivals::GetPhotonArrivalTimes(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _times[SiPMNumber];
  }

  const std::vector<CrvPhotonArrivals::SinglePhoton> &CrvPhotonArrivals::GetPhotonArrivalTimes(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _times[SiPMNumber];
  }

  unsigned int CrvPhotonArrivals::GetNumberOfPhotonArrivals(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _times[SiPMNumber].size();
  }

  unsigned int CrvPhotonArrivals::GetNumberOfPhotonArrivals(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _times[SiPMNumber].size();
  }

  double CrvPhotonArrivals::GetFirstPhotonArrivalTime() const
  {
    double firstTime = NAN;
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      for(size_t i=0; i<_times[SiPM].size(); i++)
      {
        double t = _times[SiPM][i]._time;
        if(isnan(firstTime) || t<firstTime) firstTime=t;
      }
    }
    return firstTime;
  }

  int CrvPhotonArrivals::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvPhotonArrivals::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }

}
