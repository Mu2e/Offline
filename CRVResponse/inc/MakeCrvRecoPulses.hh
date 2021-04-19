#ifndef MakeCrvRecoPulses_h
#define MakeCrvRecoPulses_h

#include <vector>
#include <TF1.h>
#include <TGraph.h>

namespace mu2eCrv
{

class MakeCrvRecoPulses
{
  public:
  MakeCrvRecoPulses(float minADCdifference, float defaultBeta, float minBeta, float maxBeta,
                    float maxTimeDifference, float minPulseHeightRatio, float maxPulseHeightRatio,
                    float LEtimeFactor, bool allowDoubleGumbel, float doubleGumbelThreshold);
  void         SetWaveform(const std::vector<unsigned int> &waveform, unsigned int startTDC, 
                           float digitizationPeriod, float pedestal, float calibrationFactor, 
                           float calibrationFactorPulseHeight);

  const std::vector<float>  &GetPEs() const            {return _PEs;}
  const std::vector<float>  &GetPEsPulseHeight() const {return _PEsPulseHeight;}
  const std::vector<double> &GetPulseTimes() const     {return _pulseTimes;}
  const std::vector<double> &GetLEtimes() const        {return _LEtimes;}
  const std::vector<float>  &GetPulseHeights() const   {return _pulseHeights;}
  const std::vector<float>  &GetPulseBetas() const     {return _pulseBetas;}
  const std::vector<float>  &GetPulseFitChi2s() const  {return _pulseFitChi2s;}
  const std::vector<bool>   &GetFailedFits() const     {return _failedFits;}

  private:
  MakeCrvRecoPulses();
  void FillGraphAndFindPeaks(const std::vector<unsigned int> &waveform, unsigned int startTDC,
                             float digitizationPeriod, float pedestal,
                             TGraph &g, std::vector<std::pair<size_t,size_t> > &peaks);
  void RangeFinderNarrow(const std::vector<unsigned int> &waveform, const size_t peakStart, const size_t peakEnd, size_t &start, size_t &end);
  void RangeFinder(const std::vector<unsigned int> &waveform, const size_t peakStart, const size_t peakEnd, size_t &start, size_t &end);
  bool FailedFit(TFitResultPtr fr, int paramStart, int paramEnd);
  double Chi2(TF1 &f, const TGraph &g);

  TF1    _f1, _f2;
  float  _minADCdifference;
  float  _defaultBeta;
  float  _minBeta, _maxBeta;
  float  _maxTimeDifference;
  float  _minPulseHeightRatio, _maxPulseHeightRatio;
  float  _LEtimeFactor;
  bool   _allowDoubleGumbel;
  float  _doubleGumbelThreshold;

  std::vector<float>  _PEs, _PEsPulseHeight;
  std::vector<double> _pulseTimes, _LEtimes; 
  std::vector<float>  _pulseHeights, _pulseBetas, _pulseFitChi2s;
  std::vector<bool>   _failedFits;

  public:
  const std::vector<float>  &GetPEsNoFit() const        {return _PEsNoFit;}
  const std::vector<double> &GetPulseTimesNoFit() const {return _pulseTimesNoFit;}
  const std::vector<double> &GetPulseStarts() const     {return _pulseStart;}
  const std::vector<double> &GetPulseEnds() const       {return _pulseEnd;}

  private:
  void NoFitOption(const std::vector<unsigned int> &waveform, float pedestal, 
                   size_t peakStart, float &sum, size_t &pulseStart, size_t &pulseEnd);
  std::vector<float>  _PEsNoFit;
  std::vector<double> _pulseTimesNoFit;
  std::vector<double> _pulseStart;
  std::vector<double> _pulseEnd;
};

}

#endif
