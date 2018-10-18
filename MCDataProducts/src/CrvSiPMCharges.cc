#include "MCDataProducts/inc/CrvSiPMCharges.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace mu2e 
{

  std::vector<mu2e::CrvSiPMCharges::CrvSingleCharge> &CrvSiPMCharges::GetSiPMCharges(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvSiPMCharges[SiPMNumber];
  }

  std::vector<mu2e::CrvSiPMCharges::CrvSingleCharge> &CrvSiPMCharges::GetSiPMCharges(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvSiPMCharges[SiPMNumber];
  }

  const std::vector<mu2e::CrvSiPMCharges::CrvSingleCharge> &CrvSiPMCharges::GetSiPMCharges(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvSiPMCharges[SiPMNumber];
  }

  const std::vector<mu2e::CrvSiPMCharges::CrvSingleCharge> &CrvSiPMCharges::GetSiPMCharges(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvSiPMCharges[SiPMNumber];
  }

  size_t CrvSiPMCharges::GetNumberOfSiPMCharges(int fiberNumber, int side) const
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvSiPMCharges[SiPMNumber].size();
  }

  size_t CrvSiPMCharges::GetNumberOfSiPMCharges(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvSiPMCharges[SiPMNumber].size();
  }

  bool CrvSiPMCharges::IsEmpty() const
  {
    for(int SiPM=0; SiPM<4; SiPM++)
    {
      if(_crvSiPMCharges[SiPM].size()!=0) return false;
    }
    return true;
  }

  double CrvSiPMCharges::GetFirstSiPMChargeTime(int fiberNumber, int side) const
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    if(_crvSiPMCharges[SiPMNumber].size()==0) return NAN;
    return _crvSiPMCharges[SiPMNumber].front()._time;
  }

  double CrvSiPMCharges::GetFirstSiPMChargeTime(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    if(_crvSiPMCharges[SiPMNumber].size()==0) return NAN;
    return _crvSiPMCharges[SiPMNumber].front()._time;
  }

  int CrvSiPMCharges::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvSiPMCharges::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }
}

