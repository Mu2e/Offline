#ifndef MCDataProducts_CrvSiPMRespones_hh
#define MCDataProducts_CrvSiPMRespones_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include <vector>
#include <cmath>

namespace mu2e 
{
  class CrvSiPMResponses
  {
    public:

    CrvSiPMResponses() {}

    struct CrvSingleSiPMResponse
    {
      double _time;
      double _charge;
      CrvSingleSiPMResponse(double time, double charge) : _time(time), _charge(charge) {}
      CrvSingleSiPMResponse() : _time(NAN), _charge(NAN) {}  //to make ROOT happy
    };

    std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int fiberNumber, int side);
    std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int SiPMNumber);

    const std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int fiberNumber, int side) const;
    const std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int SiPMNumber) const;

    unsigned int GetNumberOfSiPMResponses(int fiberNumber, int side);
    unsigned int GetNumberOfSiPMResponses(int SiPMNumber);

    double GetFirstSiPMResponseTime() const;
    double GetLastSiPMResponseTime() const;

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<CrvSingleSiPMResponse> _crvSiPMResponses[4];
  };
}

#endif /* MCDataProducts_CrvSiPMRespones_hh */
