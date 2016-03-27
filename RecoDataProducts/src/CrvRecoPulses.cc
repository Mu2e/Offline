#include "RecoDataProducts/inc/CrvRecoPulses.hh"
#include <stdexcept>
#include <algorithm>

namespace mu2e 
{

  std::vector<mu2e::CrvRecoPulses::CrvSingleRecoPulse> &CrvRecoPulses::GetRecoPulses(int fiberNumber, int side) 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvPulses[SiPMNumber];
  }

  std::vector<mu2e::CrvRecoPulses::CrvSingleRecoPulse> &CrvRecoPulses::GetRecoPulses(int SiPMNumber) 
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvPulses[SiPMNumber];
  }

  const std::vector<mu2e::CrvRecoPulses::CrvSingleRecoPulse> &CrvRecoPulses::GetRecoPulses(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvPulses[SiPMNumber];
  }

  const std::vector<mu2e::CrvRecoPulses::CrvSingleRecoPulse> &CrvRecoPulses::GetRecoPulses(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvPulses[SiPMNumber];
  }

  unsigned int CrvRecoPulses::GetNumberOfRecoPulses(int fiberNumber, int side) const 
  {
    int SiPMNumber = FindSiPMNumber(fiberNumber, side);
    return _crvPulses[SiPMNumber].size();
  }

  unsigned int CrvRecoPulses::GetNumberOfRecoPulses(int SiPMNumber) const
  {
    CheckSiPMNumber(SiPMNumber);
    return _crvPulses[SiPMNumber].size();
  }

  int CrvRecoPulses::FindSiPMNumber(int fiberNumber, int side)
  {
    if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
    if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
    int SiPM = 2*fiberNumber + side;
    return SiPM;
  }

  void CrvRecoPulses::CheckSiPMNumber(int SiPMNumber)
  {
    if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
  }

}

