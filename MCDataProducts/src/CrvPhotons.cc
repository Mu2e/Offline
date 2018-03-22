#include "MCDataProducts/inc/CrvPhotons.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<CrvPhotons::SinglePhoton> &CrvPhotons::GetPhotons(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _photons[SiPMNumber];
  }

  std::vector<CrvPhotons::SinglePhoton> &CrvPhotons::GetPhotons(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _photons[SiPMNumber];
  }

  const std::vector<CrvPhotons::SinglePhoton> &CrvPhotons::GetPhotons(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _photons[SiPMNumber];
  }

  const std::vector<CrvPhotons::SinglePhoton> &CrvPhotons::GetPhotons(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _photons[SiPMNumber];
  }

  size_t CrvPhotons::GetNumberOfPhotons(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _photons[SiPMNumber].size();
  }

  size_t CrvPhotons::GetNumberOfPhotons(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _photons[SiPMNumber].size();
  }

  double CrvPhotons::GetFirstPhotonTime() const
  {
    double firstTime = NAN;
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      for(size_t i=0; i<_photons[SiPM].size(); i++)
      {
        double t = _photons[SiPM][i]._time;
        if(std::isnan(firstTime) || t<firstTime) firstTime=t;
      }
    }
    return firstTime;
  }

  int CrvPhotons::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvPhotons::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }

}
