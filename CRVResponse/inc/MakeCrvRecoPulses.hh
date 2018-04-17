#ifndef MakeCrvRecoPulses_h
#define MakeCrvRecoPulses_h

#include <vector>

namespace mu2eCrv
{

class MakeCrvRecoPulses
{
  public:
  MakeCrvRecoPulses();
  void         SetWaveform(const std::vector<unsigned int> &waveform, unsigned int startTDC, 
                           double digitizationPeriod, double pedestal, double calibrationFactor, bool darkNoise);
  unsigned int GetNPulses();
  int          GetPEs(int pulse);
  double       GetPulseTime(int pulse);
  double       GetPulseHeight(int pulse);
  double       GetPulseWidth(int pulse);
  double       GetPulseFitChi2(int pulse);
  double       GetFitParam0(int pulse);
  double       GetFitParam1(int pulse);
  double       GetFitParam2(int pulse);
  double       GetT1(int pulse);
  double       GetT2(int pulse);
  double       GetLEtime(int pulse);
  double       GetLEfitChi2(int pulse);
  int          GetPeakBin(int pulse);

  private:
  std::vector<int>    _PEs;
  std::vector<double> _pulseTimes, _pulseHeights, _pulseWidths, _pulseFitChi2s;
  std::vector<double> _fitParams0, _fitParams1, _fitParams2, _t1s, _t2s;
  std::vector<double> _LEtimes, _LEfitChi2s;
  std::vector<int>    _peakBins;
};

}

#endif
