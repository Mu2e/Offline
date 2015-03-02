#ifndef MCDataProducts_CrvPhotonArrivals_hh
#define MCDataProducts_CrvPhotonArrivals_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include <vector>

namespace mu2e 
{
  class CrvPhotonArrivals
  {
    public:

    CrvPhotonArrivals() {}

    std::vector<double> &GetPhotonArrivalTimes(int fiberNumber, int side); 
    std::vector<double> &GetPhotonArrivalTimes(int SiPMNumber);

    const std::vector<double> &GetPhotonArrivalTimes(int fiberNumber, int side) const;
    const std::vector<double> &GetPhotonArrivalTimes(int SiPMNumber) const;

    unsigned int GetNumberOfPhotonArrivals(int fiberNumber, int side) const; 
    unsigned int GetNumberOfPhotonArrivals(int SiPMNumber) const;

    double GetFirstPhotonArrivalTime() const;

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<double>  _times[4];
  };
}

#endif /* MCDataProducts_CrvPhotonArrivals_hh */
