#ifndef MCDataProducts_CrvSiPMCharges_hh
#define MCDataProducts_CrvSiPMCharges_hh
//
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include <vector>
#include <cmath>

namespace mu2e 
{
  class CrvSiPMCharges
  {
    public:

    CrvSiPMCharges() {}

    struct CrvSingleCharge
    {
      double _time;
      double _charge;
      double _chargeInPEs;
      art::Ptr<StepPointMC> _step;
      CrvSingleCharge(double time, double charge, double chargeInPEs, art::Ptr<StepPointMC> step) : 
                                               _time(time), _charge(charge), _chargeInPEs(chargeInPEs), _step(step) {}
      CrvSingleCharge(double time, double charge, double chargeInPEs) :      //that's for dark noise (i.e. no StepPointMCs)
                                               _time(time), _charge(charge), _chargeInPEs(chargeInPEs) {}
      CrvSingleCharge() : _time(NAN), _charge(NAN), _chargeInPEs(NAN) {}  //to make ROOT happy
    };

    std::vector<CrvSingleCharge> &GetSiPMCharges(int fiberNumber, int side);
    std::vector<CrvSingleCharge> &GetSiPMCharges(int SiPMNumber);

    const std::vector<CrvSingleCharge> &GetSiPMCharges(int fiberNumber, int side) const;
    const std::vector<CrvSingleCharge> &GetSiPMCharges(int SiPMNumber) const;

    size_t GetNumberOfSiPMCharges(int fiberNumber, int side) const;
    size_t GetNumberOfSiPMCharges(int SiPMNumber) const;

    bool IsEmpty() const;

    double GetFirstSiPMChargeTime(int fiberNumber, int side) const;
    double GetFirstSiPMChargeTime(int SiPMNumber) const;

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<CrvSingleCharge> _crvSiPMCharges[4];
  };
}

#endif /* MCDataProducts_CrvSiPMCharges_hh */
