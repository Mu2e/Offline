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
                    float LEtimeFactor, float pulseThreshold, float pulseAreaThreshold, float doublePulseSeparation);
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
  const std::vector<bool>   &GetZeroNdfs() const       {return _zeroNdf;}
  const std::vector<bool>   &GetFailedFits() const     {return _failedFits;}

  private:
  MakeCrvRecoPulses();
  void FillGraphAndFindPeaks(const std::vector<unsigned int> &waveform, unsigned int startTDC,
                             float digitizationPeriod, float pedestal,
                             TGraph &g, std::vector<std::pair<size_t,size_t> > &peaks);
  void RangeFinder(const std::vector<unsigned int> &waveform, const size_t peakStart, const size_t peakEnd, size_t &start, size_t &end);
  bool FailedFit(TFitResultPtr fr);

  TF1    _f1;
  float  _minADCdifference;
  float  _defaultBeta;
  float  _minBeta, _maxBeta;
  float  _maxTimeDifference;
  float  _minPulseHeightRatio, _maxPulseHeightRatio;
  float  _LEtimeFactor;
  float  _pulseThreshold;
  float  _pulseAreaThreshold;
  float  _doublePulseSeparation;

  std::vector<float>  _PEs, _PEsPulseHeight;
  std::vector<double> _pulseTimes, _LEtimes;
  std::vector<float>  _pulseHeights, _pulseBetas, _pulseFitChi2s;
  std::vector<bool>   _zeroNdf, _failedFits, _duplicateNoFitPulses, _separatedDoublePulses;

  public:
  const std::vector<float>  &GetPEsNoFit() const        {return _PEsNoFit;}
  const std::vector<double> &GetPulseTimesNoFit() const {return _pulseTimesNoFit;}
  const std::vector<double> &GetPulseStarts() const     {return _pulseStart;}
  const std::vector<double> &GetPulseEnds() const       {return _pulseEnd;}
  const std::vector<bool>   &GetDuplicateNoFitPulses() const  {return _duplicateNoFitPulses;}
  const std::vector<bool>   &GetSeparatedDoublePulses() const {return _separatedDoublePulses;}

  private:
  void NoFitOption(const std::vector<unsigned int> &waveform, const std::vector<std::pair<size_t,size_t> > &peaks,
                   unsigned int startTDC, float digitizationPeriod, float pedestal, float calibrationFactor);
  std::vector<float>  _PEsNoFit;
  std::vector<double> _pulseTimesNoFit;
  std::vector<double> _pulseStart;
  std::vector<double> _pulseEnd;
};

}

#endif
