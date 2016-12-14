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
  MakeCrvRecoPulses(double pulseThreshold, double leadingEdgeThreshold, double param0, double param1);
  void         SetWaveform(const std::vector<double> &waveform, double startTime, double binWidth);
  unsigned int GetNPulses();
  int          GetPEs(int pulse);
  double       GetLeadingEdge(int pulse);
  double       GetTimeOverThreshold(int pulse);
  double       GetPulseHeight(int pulse);
  double       GetPulseHeightLandau(int pulse);
  double       GetPeakTime(int pulse);
  double       GetPeakTimeLandau(int pulse);
  double       GetIntegral(int pulse);
  double       GetLandauParam0(int pulse);
  double       GetLandauParam1(int pulse);
  double       GetLandauParam2(int pulse);
  double       GetT1(int pulse);
  double       GetT2(int pulse);

  private:
  double _pulseThreshold;
  double _leadingEdgeThreshold;
  double _param0, _param1;

  std::vector<int>    _PEs;
  std::vector<double> _leadingEdges;
  std::vector<double> _pulseHeights, _pulseHeightsLandau;
  std::vector<double> _peakTimes, _peakTimesLandau;
  std::vector<double> _integrals, _landauParams0, _landauParams1, _landauParams2, _T1s, _T2s, _TOTs;
};

}

#endif
