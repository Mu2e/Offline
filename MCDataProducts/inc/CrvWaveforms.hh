#ifndef MCDataProducts_CrvWaveforms_hh
#define MCDataProducts_CrvWaveforms_hh
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
  class CrvWaveforms
  {
    public:

    CrvWaveforms() {}

    struct CrvSingleWaveform
    {
      std::vector<double> _voltages;
      double _startTime;
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

#endif /* MCDataProducts_CrvWaveforms_hh */
