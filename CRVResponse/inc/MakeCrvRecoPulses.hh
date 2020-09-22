#ifndef MakeCrvRecoPulses_h
#define MakeCrvRecoPulses_h

#include <vector>
#include <TF1.h>

namespace mu2eCrv
{

class MakeCrvRecoPulses
{
  public:
  MakeCrvRecoPulses(double minADCdifference, double defaultBeta, double minBeta, double maxBeta,
                    double maxTimeDifference, double minPulseHeightRatio, double maxPulseHeightRatio,
                    double LEtimeFactor);
  void         SetWaveform(const std::vector<unsigned int> &waveform, unsigned int startTDC, 
                           double digitizationPeriod, double pedestal, double calibrationFactor, 
                           double calibrationFactorPulseHeight);
  unsigned int GetNPulses();
  int          GetPEs(int pulse);
  int          GetPEsPulseHeight(int pulse);
  double       GetPulseTime(int pulse);
  double       GetPulseHeight(int pulse);
  double       GetPulseBeta(int pulse);
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
  MakeCrvRecoPulses();

  TF1    _f;
  double _minADCdifference;
  double _defaultBeta;
  double _minBeta, _maxBeta;
  double _maxTimeDifference;
  double _minPulseHeightRatio, _maxPulseHeightRatio;
  double _LEtimeFactor;

  std::vector<int>    _PEs, _PEsPulseHeight;
  std::vector<double> _pulseTimes, _pulseHeights, _pulseBetas, _pulseFitChi2s;
  std::vector<double> _fitParams0, _fitParams1, _fitParams2, _t1s, _t2s;
  std::vector<double> _LEtimes, _LEfitChi2s;
  std::vector<int>    _peakBins;
};

}

#endif
