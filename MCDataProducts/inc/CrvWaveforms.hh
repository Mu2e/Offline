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

    std::vector<double> &GetWaveform(int fiberNumber, int side);
    std::vector<double> &GetWaveform(int SiPMNumber);

    const std::vector<double> &GetWaveform(int fiberNumber, int side) const;
    const std::vector<double> &GetWaveform(int SiPMNumber) const;

    double GetStartTime(int fiberNumber, int side) const;
    double GetStartTime(int SiPMNumber) const;

    void SetStartTime(int fiberNumber, int side, double startTime);
    void SetStartTime(int SiPMNumber, double startTime);

    double GetBinWidth() const; 

    void SetBinWidth(double binWidth);

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<double> _waveform[4];
    double _startTime[4];
    double _binWidth;
  };
}

#endif /* MCDataProducts_CrvWaveforms_hh */
