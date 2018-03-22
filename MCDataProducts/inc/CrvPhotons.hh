#ifndef MCDataProducts_CrvPhotons_hh
#define MCDataProducts_CrvPhotons_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include <vector>

namespace mu2e 
{
  class CrvPhotons
  {
    public:

    struct SinglePhoton
    {
      double _time;
      art::Ptr<StepPointMC> _step;
    };

    CrvPhotons() {}

    std::vector<SinglePhoton> &GetPhotons(int fiberNumber, int side); 
    std::vector<SinglePhoton> &GetPhotons(int SiPMNumber);

    const std::vector<SinglePhoton> &GetPhotons(int fiberNumber, int side) const;
    const std::vector<SinglePhoton> &GetPhotons(int SiPMNumber) const;

    size_t GetNumberOfPhotons(int fiberNumber, int side) const; 
    size_t GetNumberOfPhotons(int SiPMNumber) const;

    double GetFirstPhotonTime() const;

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<SinglePhoton> _photons[4];
  };
}

#endif /* MCDataProducts_CrvPhotons_hh */
