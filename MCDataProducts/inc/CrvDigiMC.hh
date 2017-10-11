#ifndef MCDataProducts_CrvDigiMC_hh
#define MCDataProducts_CrvDigiMC_hh
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
  class CrvDigiMC
  {
    public:

    CrvDigiMC() {}

    struct CrvSingleWaveform
    {
      std::vector<double>                 _voltages;
      std::vector<art::Ptr<StepPointMC> > _steps;
      art::Ptr<SimParticle>               _simparticle;
      double                              _startTime;
    };

    std::vector<CrvSingleWaveform> &GetSingleWaveforms(int fiberNumber, int side);
    std::vector<CrvSingleWaveform> &GetSingleWaveforms(int SiPMNumber);

    const std::vector<CrvSingleWaveform> &GetSingleWaveforms(int fiberNumber, int side) const;
    const std::vector<CrvSingleWaveform> &GetSingleWaveforms(int SiPMNumber) const;

    double GetDigitizationPrecision() const; 
    void SetDigitizationPrecision(double precision);

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<CrvSingleWaveform> _waveforms[4];
    double _digitizationPrecision;
  };
}

#endif /* MCDataProducts_CrvDigiMC_hh */
