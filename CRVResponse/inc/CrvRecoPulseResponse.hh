#include <vector>

class CrvRecoPulseResponse
{
  private:
  CrvRecoPulseResponse();

  public:
  CrvRecoPulseResponse(double pulseThreshold, double leadingEdgeThreshold, double integralFactor);
  void         SetWaveform(const std::vector<double> &waveform, double startTime, double binWidth);
  unsigned int GetNPulses();
  double       GetPEs(int pulse);
  double       GetLeadingEdge(int pulse);
  double       GetPulseHeight(int pulse);
  double       GetIntegral(int pulse);
  double       GetLandauParam0(int pulse);
  double       GetLandauParam1(int pulse);
  double       GetLandauParam2(int pulse);
  double       GetT1(int pulse);
  double       GetT2(int pulse);

  private:
  double _pulseThreshold;
  double _leadingEdgeThreshold;
  double _integralFactor;

  std::vector<int>    _PEs;
  std::vector<double> _leadingEdges;
  std::vector<double> _pulseHeights;
  std::vector<double> _integrals, _landauParams0, _landauParams1, _landauParams2, _T1s, _T2s;
};

