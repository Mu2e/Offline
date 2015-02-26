#ifndef MCDataProducts_CrvPhotonArrivals_hh
#define MCDataProducts_CrvPhotonArrivals_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "StepPointMC.hh"

namespace mu2e 
{
  class CrvPhotonArrivals
  {
    public:

    CrvPhotonArrivals() {}

    std::vector<double> &GetPhotonArrivalTimes(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _times[SiPM];
    }

    std::vector<double> &GetPhotonArrivalTimes(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _times[SiPMNumber];
    }

    const std::vector<double> &GetPhotonArrivalTimes(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _times[SiPM];
    }

    const std::vector<double> &GetPhotonArrivalTimes(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _times[SiPMNumber];
    }

    unsigned int GetNumberOfPhotonArrivals(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _times[SiPM].size();
    }

    unsigned int GetNumberOfPhotonArrivals(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _times[SiPMNumber].size();
    }

    double GetFirstPhotonArrivalTime() const
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

    private:

    std::vector<double>  _times[4];
  };
}

#endif /* MCDataProducts_CrvPhotonArrivals_hh */
