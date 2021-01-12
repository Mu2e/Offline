#ifndef MCDataProducts_CrvSiPMCharges_hh
#define MCDataProducts_CrvSiPMCharges_hh
//
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/CrvStep.hh"
#include <vector>

namespace mu2e 
{
  class CrvSiPMCharges
  {
    public:

    CrvSiPMCharges() : _SiPMNumber(-1) {} //automatically sets the barindex to invalid

    CrvSiPMCharges(mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
                                         _scintillatorBarIndex(scintillatorBarIndex),
                                         _SiPMNumber(SiPMNumber) {}

    struct SingleCharge
    {
      double _time;
      double _charge;
      double _chargeInPEs;
      art::Ptr<CrvStep> _step;
      SingleCharge(double time, double charge, double chargeInPEs, art::Ptr<CrvStep> step) : 
                                               _time(time), _charge(charge), _chargeInPEs(chargeInPEs), _step(step) {}
      SingleCharge(double time, double charge, double chargeInPEs) :   //that's for dark noise (i.e. no StepPointMCs)
                                               _time(time), _charge(charge), _chargeInPEs(chargeInPEs) {}
      SingleCharge() : _time(NAN), _charge(NAN), _chargeInPEs(NAN) {}  //to make ROOT happy
    };

    mu2e::CRSScintillatorBarIndex    GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                              GetSiPMNumber() const           {return _SiPMNumber;}
    std::vector<SingleCharge>       &GetCharges()                    {return _charges;}
    const std::vector<SingleCharge> &GetCharges() const              {return _charges;}

    private:

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
    std::vector<SingleCharge>      _charges;
  };

  typedef std::vector<mu2e::CrvSiPMCharges> CrvSiPMChargesCollection;
}

#endif /* MCDataProducts_CrvSiPMCharges_hh */
