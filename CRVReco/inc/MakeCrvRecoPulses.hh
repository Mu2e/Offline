#ifndef MakeCrvRecoPulses_h
#define MakeCrvRecoPulses_h

#include <array>
#include <cstdint>
#include <vector>

namespace mu2eCrv
{

class MakeCrvRecoPulses
{
  public:
  MakeCrvRecoPulses(float minADCdifference, float defaultBeta, float minBeta, float maxBeta,
                    float maxTimeDifference, float minPulseHeightRatio, float maxPulseHeightRatio,
                    float LEtimeFactor, int samplesBefore, int samplesAfter);
  void         SetWaveform(const std::vector<int16_t> &waveform, uint16_t startTDC,
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
  typedef std::array<double,3> Vec3;

  struct FitResult
  {
    Vec3   par;
    double chi2;
    bool   converged;
    bool   hessianOk;
  };

  MakeCrvRecoPulses();
  bool FindNextPeak(const std::vector<double> &y, int start, int &peakStart, int &peakEnd, int &fitStart, int &fitEnd) const;
  FitResult FitGumbel(const std::vector<double> &t, const std::vector<double> &y,
                      const Vec3 &p0, const Vec3 &lower, const Vec3 &upper) const;

  float  _minADCdifference;
  float  _defaultBeta;
  float  _minBeta, _maxBeta;
  float  _maxTimeDifference;
  float  _minPulseHeightRatio, _maxPulseHeightRatio;
  float  _LEtimeFactor;
  int    _samplesBefore, _samplesAfter;

  std::vector<float>  _PEs, _PEsPulseHeight;
  std::vector<double> _pulseTimes, _LEtimes;
  std::vector<float>  _pulseHeights, _pulseBetas, _pulseFitChi2s;
  std::vector<bool>   _zeroNdf, _failedFits;
};

}

#endif
