#include "Offline/CRVResponse/inc/MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>

namespace
{
  double Gumbel(double* xs, double* par)
  {
    double const x = xs[0];
    return par[0]*(TMath::Exp(-(x-par[1])/par[2]-TMath::Exp(-(x-par[1])/par[2])));
  }
}

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses(float minADCdifference, float defaultBeta, float minBeta, float maxBeta,
                                     float maxTimeDifference, float minPulseHeightRatio, float maxPulseHeightRatio,
                                     float LEtimeFactor, float pulseThreshold, float pulseAreaThreshold, float doublePulseSeparation) :
                                     _f1("peakfitter",Gumbel,0,0,3),
                                     _minADCdifference(minADCdifference),
                                     _defaultBeta(defaultBeta), _minBeta(minBeta), _maxBeta(maxBeta),
                                     _maxTimeDifference(maxTimeDifference),
                                     _minPulseHeightRatio(minPulseHeightRatio),
                                     _maxPulseHeightRatio(maxPulseHeightRatio),
                                     _LEtimeFactor(LEtimeFactor),
                                     _pulseThreshold(pulseThreshold),
                                     _pulseAreaThreshold(pulseAreaThreshold),
                                     _doublePulseSeparation(doublePulseSeparation)
{}

void MakeCrvRecoPulses::FillGraphAndFindPeaks(const std::vector<unsigned int> &waveform, unsigned int startTDC,
                                              float digitizationPeriod, float pedestal,
                                              TGraph &g, std::vector<std::pair<size_t,size_t> > &peaks)
{
  size_t nBins = waveform.size();
  size_t peakStartBin=0;
  size_t peakEndBin=0;
  for(size_t bin=0; bin<nBins; ++bin)
  {
    g.SetPoint(bin,(startTDC+bin)*digitizationPeriod,waveform[bin]-pedestal);

    if(bin<1) continue; //don't search for peaks here

    if(waveform[bin-1]<waveform[bin]) //rising edge
    {
      peakStartBin=bin;
      peakEndBin=bin;
    }
    if(waveform[bin-1]==waveform[bin]) //potentially a peak where consecutive ADC values are equal
    {
      peakEndBin=bin;
    }
    if(waveform[bin-1]>waveform[bin]) //falling edge
    {
      if(peakStartBin>0)  //found a peak
      {
        if(waveform[peakEndBin]-pedestal>_minADCdifference) peaks.emplace_back(peakStartBin,peakEndBin);
        peakStartBin=0;  //so that the loop has to wait for the next rising edge
      }
    }
  }
}

void MakeCrvRecoPulses::RangeFinder(const std::vector<unsigned int> &waveform, const size_t peakStart, const size_t peakEnd, size_t &start, size_t &end)
{
  assert(peakStart>0);
  assert(peakEnd<waveform.size()-1);

  end=peakEnd+1;
  for(size_t i=peakEnd; i<waveform.size()-1 && i<peakEnd+4; ++i)
  {
    if(waveform[i+1]<=waveform[i]) end=i+1;
    else break;
  }
  start=peakStart-1;
  for(size_t i=peakStart; i>0  && i+4>peakStart; --i)
  {
    if(waveform[i-1]<=waveform[i]) start=i-1;
    else break;
  }
}

bool MakeCrvRecoPulses::FailedFit(TFitResultPtr fr)
{
  if(fr!=0) return true;
  if(!fr->IsValid()) return true;

  const double tolerance=0.01; //TODO: Try to ask the minimizer, if the parameter is at the limit
  for(int i=0; i<=2; ++i)
  {
    double v=fr->Parameter(i);
    double lower, upper;
    fr->ParameterBounds(i,lower,upper);
    if((v-lower)/(upper-lower)<tolerance) return true;
    if((upper-v)/(upper-lower)<tolerance) return true;
  }
  return false;
}

void MakeCrvRecoPulses::NoFitOption(const std::vector<unsigned int> &waveform, const std::vector<std::pair<size_t,size_t> > &peaks,
                                    unsigned int startTDC, float digitizationPeriod, float pedestal, float calibrationFactor)
{
  //find troughs between peaks, that may be used to separate double pusles
  std::vector<size_t> troughs;
  for(size_t i=1; i<peaks.size(); ++i)
  {
    size_t peak1=peaks[i-1].first;
    size_t peak2=peaks[i].first;
    auto troughIter=std::min_element(waveform.begin()+peak1,waveform.begin()+peak2+1);
    size_t trough=std::distance(waveform.begin(),troughIter);
    if(*troughIter-pedestal<_doublePulseSeparation*(waveform[peak1]-pedestal) && *troughIter-pedestal<_doublePulseSeparation*(waveform[peak2]-pedestal)) troughs.push_back(trough);
  }

  bool    aboveAreaThreshold=false;
  bool    pulseFound=false;
  bool    doublePulseThisPeak=false;
  bool    doublePulseNextPeak=false;
  size_t  pulseStart=0, pulseEnd=0;
  float   sum=0;
  for(size_t i=0; i<waveform.size(); ++i)
  {
    float ADC=waveform[i]-pedestal;
    if(ADC>_pulseAreaThreshold && !aboveAreaThreshold) //waveform rises above area threshold
    {
      aboveAreaThreshold=true;
      pulseStart=i;
    }
    if(ADC<=_pulseAreaThreshold && aboveAreaThreshold) //waveform falls below area threshold
    {
      aboveAreaThreshold=false;
      pulseFound=true;
      pulseEnd=i-1;
    }
    if(aboveAreaThreshold) sum+=ADC;
    if(i+1==waveform.size() && aboveAreaThreshold) //end of waveform, but still above area threshold
    {
      aboveAreaThreshold=false;
      pulseFound=true;
      pulseEnd=i;
    }
    if(find(troughs.begin(),troughs.end(),i)!=troughs.end() && aboveAreaThreshold) //found the lowest point between the peaks of a double pulse
    {
      aboveAreaThreshold=false;
      pulseFound=true;
      pulseEnd=i;
      doublePulseThisPeak=true;
    }
    if(pulseFound)  //a full waveform section about min threshold has been found
    {
      std::vector<float> peaksInCurrentPulse;
      for(auto ipeak=peaks.begin(); ipeak!=peaks.end(); ++ipeak)
      {
        if(ipeak->first>=pulseStart && ipeak->second<=pulseEnd)  //found a peak within this pulse
        {
          float peakADC=waveform[ipeak->first]-pedestal;
          peaksInCurrentPulse.push_back(peakADC);
          _pulseTimesNoFit.push_back((startTDC+0.5*(ipeak->first+ipeak->second))*digitizationPeriod);
          _PEsNoFit.push_back(sum*digitizationPeriod/calibrationFactor);  //every peak of this pulse gets the same PEs
          _separatedDoublePulses.push_back(doublePulseThisPeak || doublePulseNextPeak);
        }
      }

      pulseFound=false;
      sum=0;
      doublePulseNextPeak=false;
      if(doublePulseThisPeak)
      {
        doublePulseThisPeak=false;
        doublePulseNextPeak=true;
        sum=ADC; //the shared trough point gets added to both pulses
      }

      if(peaksInCurrentPulse.empty()) continue;

      //every peak of this pulse gets the same start/end pulse time:
      //the times when it crosses the threshold fraction of the lowest peak of this pulse
      float  minPeakInCurrentPulse=*std::min_element(peaksInCurrentPulse.begin(),peaksInCurrentPulse.end());
      float  rangeThreshold=_pulseThreshold*minPeakInCurrentPulse;
      size_t pulseRangeStart=pulseStart;
      size_t pulseRangeEnd=pulseEnd;
      for(size_t j=pulseStart; j<=pulseEnd; ++j)
      {
//        if(waveform[j]-pedestal>rangeThreshold) {pulseRangeStart=j; break;}
        if(waveform[j]-pedestal>rangeThreshold) {pulseRangeStart=(j>pulseStart?j-1:j); break;}  //include point before it crosses the threshold
      }
      for(size_t j=pulseEnd; j>=pulseStart; --j)
      {
//        if(waveform[j]-pedestal>rangeThreshold) {pulseRangeEnd=j; break;}
        if(waveform[j]-pedestal>rangeThreshold) {pulseRangeEnd=(j<pulseEnd?j+1:j); break;}  //include point after it crosses the threshold
        if(j==0) break;
      }
      for(size_t j=0; j<peaksInCurrentPulse.size(); ++j)
      {
        _pulseStart.push_back((startTDC+pulseRangeStart)*digitizationPeriod);
        _pulseEnd.push_back((startTDC+pulseRangeEnd)*digitizationPeriod);
        _duplicateNoFitPulses.push_back(j==0?false:true);
      }
    }
  }
}

void MakeCrvRecoPulses::SetWaveform(const std::vector<unsigned int> &waveform,
                                    unsigned int startTDC, float digitizationPeriod, float pedestal,
                                    float calibrationFactor, float calibrationFactorPulseHeight)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseBetas.clear();
  _pulseFitChi2s.clear();
  _PEs.clear();
  _PEsPulseHeight.clear();
  _LEtimes.clear();
  _failedFits.clear();
  _duplicateNoFitPulses.clear();
  _separatedDoublePulses.clear();
  _PEsNoFit.clear();
  _pulseTimesNoFit.clear();
  _pulseStart.clear();
  _pulseEnd.clear();

  //fill graph and find peaks
  std::vector<std::pair<size_t,size_t> > peaks;
  size_t nBins = waveform.size();
  TGraph g(nBins);
  FillGraphAndFindPeaks(waveform, startTDC, digitizationPeriod, pedestal, g, peaks);

  //loop through all peaks
  for(size_t ipeak=0; ipeak<peaks.size(); ++ipeak)
  {
    size_t peakStartBin=peaks[ipeak].first;
    size_t peakEndBin=peaks[ipeak].second;
    double peakStartTime=(startTDC+peakStartBin)*digitizationPeriod;
    double peakEndTime=(startTDC+peakEndBin)*digitizationPeriod;
    double peakTime=0.5*(peakStartTime+peakEndTime);

    _f1.SetParameter(0, (waveform[peakStartBin]-pedestal)*TMath::E());
    _f1.SetParameter(1, peakTime);
    _f1.SetParameter(2, _defaultBeta);
    _f1.SetParLimits(0,(waveform[peakStartBin]-pedestal)*TMath::E()*_minPulseHeightRatio,(waveform[peakStartBin]-pedestal)*TMath::E()*_maxPulseHeightRatio);
    _f1.SetParLimits(1,peakStartTime-_maxTimeDifference,peakEndTime+_maxTimeDifference);
    _f1.SetParLimits(2, _minBeta, _maxBeta);

    size_t fitStartBin, fitEndBin;
    RangeFinder(waveform, peakStartBin, peakEndBin, fitStartBin, fitEndBin);
    double fitStartTime=(startTDC+fitStartBin)*digitizationPeriod;
    double fitEndTime=(startTDC+fitEndBin)*digitizationPeriod;
    _f1.SetRange(fitStartTime,fitEndTime);

    //do the fit
    TFitResultPtr fr = g.Fit(&_f1,"NQSR");
    double fitParam0 = fr->Parameter(0);
    double fitParam1 = fr->Parameter(1);
    double fitParam2 = fr->Parameter(2);

    //collect fit information for the first peak
    float  PEs          = fitParam0*fitParam2 / calibrationFactor;
    double pulseTime    = fitParam1;
    float  pulseHeight  = fitParam0/TMath::E();
    float  pulseBeta    = fitParam2;
    float  pulseFitChi2 = (fr->Ndf()>0?fr->Chi2()/fr->Ndf():NAN);
    bool   failedFit    = FailedFit(fr);

    if(failedFit)
    {
      PEs          = (waveform[peakStartBin]-pedestal)*TMath::E() * _defaultBeta / calibrationFactor;
      pulseTime    = peakTime;
      pulseHeight  = waveform[peakStartBin]-pedestal;
      pulseBeta    = _defaultBeta;
      pulseFitChi2 = NAN;
    }

    double LEtime         = pulseTime-_LEtimeFactor*pulseBeta;  //50% pulse height is reached at -0.985*beta before the peak
    float  PEsPulseHeight = pulseHeight / calibrationFactorPulseHeight;

    _pulseTimes.push_back(pulseTime);
    _pulseHeights.push_back(pulseHeight);
    _pulseBetas.push_back(pulseBeta);
    _pulseFitChi2s.push_back(pulseFitChi2);
    _PEs.push_back(PEs);
    _PEsPulseHeight.push_back(PEsPulseHeight);
    _LEtimes.push_back(LEtime);
    _failedFits.push_back(failedFit);
  }

  NoFitOption(waveform, peaks, startTDC, digitizationPeriod, pedestal, calibrationFactor);

}

}

