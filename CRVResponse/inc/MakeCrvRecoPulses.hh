#ifndef MakeCrvRecoPulses_h
#define MakeCrvRecoPulses_h

#include <vector>

namespace mu2eCrv
{

class MakeCrvRecoPulses
{
  private:
  MakeCrvRecoPulses();

  public:
  MakeCrvRecoPulses(double scale, bool useFittedPulseHeight, bool useFittedPulseTime);
  void         SetWaveform(const std::vector<double> &waveform, double startTime, double binWidth);
  unsigned int GetNPulses();
  int          GetPEs(int pulse);
  double       GetPulseTime(int pulse);
  double       GetPulseHeight(int pulse);
  double       GetPulseWidth(int pulse);
  double       GetFitParam0(int pulse);
  double       GetFitParam1(int pulse);
  double       GetFitParam2(int pulse);
  double       GetT1(int pulse);
  double       GetT2(int pulse);

  private:
  double _scale;
  bool   _useFittedPulseHeight, _useFittedPulseTime;

  std::vector<int>    _PEs;
  std::vector<double> _pulseTimes, _pulseHeights, _pulseWidths;
  std::vector<double> _fitParams0, _fitParams1, _fitParams2, _t1s, _t2s;
};

}

#endif
