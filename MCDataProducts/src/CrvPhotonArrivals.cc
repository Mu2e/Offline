#include "MCDataProducts/inc/CrvPhotonArrivals.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<double> &CrvPhotonArrivals::GetPhotonArrivalTimes(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _times[SiPMNumber];
  }

  std::vector<double> &CrvPhotonArrivals::GetPhotonArrivalTimes(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _times[SiPMNumber];
  }

  const std::vector<double> &CrvPhotonArrivals::GetPhotonArrivalTimes(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _times[SiPMNumber];
  }

  const std::vector<double> &CrvPhotonArrivals::GetPhotonArrivalTimes(int SiPMNumber) const
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
      if(_times[SiPM].size()==0) continue;
      double t = *std::min_element(_times[SiPM].begin(),_times[SiPM].end());
      if(isnan(firstTime) || t<firstTime) firstTime=t;
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
