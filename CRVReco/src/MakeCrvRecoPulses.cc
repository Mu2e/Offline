#include "Offline/CRVReco/inc/MakeCrvRecoPulses.hh"
#ifndef CRVStandalone
#include "canvas/Utilities/Exception.h"
#endif
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
                                     float LEtimeFactor, int samplesBefore, int samplesAfter) :
                                     _f1("peakfitter",Gumbel,0,0,3),
                                     _minADCdifference(minADCdifference),
                                     _defaultBeta(defaultBeta), _minBeta(minBeta), _maxBeta(maxBeta),
                                     _maxTimeDifference(maxTimeDifference),
                                     _minPulseHeightRatio(minPulseHeightRatio),
                                     _maxPulseHeightRatio(maxPulseHeightRatio),
                                     _LEtimeFactor(LEtimeFactor),
                                     _samplesBefore(samplesBefore), _samplesAfter(samplesAfter)
{
}

bool MakeCrvRecoPulses::FindNextPeak(const TGraph &g, int start, int &peakStart, int &peakEnd, int &fitStart, int &fitEnd)
{
  bool foundRisingEdge=false;
  bool foundPeak=false;
  for(int bin=start+1; bin<g.GetN(); ++bin)
  {
    if(g.GetPointY(bin-1)<g.GetPointY(bin)) //rising edge
    {
      foundRisingEdge=true;
      peakStart=bin;
      peakEnd=bin;
    }
    if(g.GetPointY(bin-1)==g.GetPointY(bin) && foundRisingEdge) //potentially a peak where consecutive ADC values are equal
    {
      peakEnd=bin;
    }
    if(g.GetPointY(bin-1)>g.GetPointY(bin) && foundRisingEdge) //falling edge
    {
      if(g.GetPointY(peakEnd)>_minADCdifference)
      {
        foundPeak=true;
        break; //found the peak. stop searching.
      }
      foundRisingEdge=false;  //so that the loop has to wait for the next rising edge
    }
  }

  if(!foundPeak) return false;

  //select a range of up to 3 points before and 2 after the peak (to avoid hidden reflection pulses)
  //-find up to 4/3 points before and after the peak for which the waveform is stricly decreasing
  //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)

  fitStart=peakStart-1;
  fitEnd=peakEnd+1;
  for(int i=peakStart-1; i>=peakStart-_samplesBefore; --i)
  {
    if(g.GetPointY(i)<=g.GetPointY(i+1)) fitStart=i;
    else break;
    if(i==0) break;
  }
  for(int i=peakEnd+1; i<g.GetN() && i<=peakEnd+_samplesAfter; ++i)
  {
    if(g.GetPointY(i-1)>=g.GetPointY(i)) fitEnd=i;
    else break;
  }
  if(peakStart-fitStart>1) fitStart++;
  if(fitEnd-peakEnd>1) fitEnd--;

  return true;
}

void MakeCrvRecoPulses::SubtractPulse(TGraph &g, uint16_t startTDC, float digitizationPeriod)
{
  for(int bin=0; bin<g.GetN(); ++bin)
  {
    double t=(startTDC+bin)*digitizationPeriod;
    double diff=std::round(_f1.Eval(t));

    g.SetPointY(bin,g.GetPointY(bin)-diff);
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

void MakeCrvRecoPulses::SetWaveform(const std::vector<int16_t> &waveform,
                                    uint16_t startTDC, float digitizationPeriod, float pedestal,
                                    float calibrationFactor, float calibrationFactorPulseHeight)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseBetas.clear();
  _pulseFitChi2s.clear();
  _PEs.clear();
  _PEsPulseHeight.clear();
  _LEtimes.clear();
  _zeroNdf.clear();
  _failedFits.clear();

  //fill the TGraph
  size_t nBins=waveform.size();
  TGraph g(nBins);
  for(size_t bin=0; bin<nBins; ++bin)
  {
    g.SetPoint(bin,(startTDC+bin)*digitizationPeriod,waveform[bin]-pedestal);
  }

  //loop through all peaks
  int start=0;
  int peakStart=0;
  int peakEnd=0;
  int fitStart=0;
  int fitEnd=0;
  while(FindNextPeak(g, start, peakStart, peakEnd, fitStart, fitEnd))
  {
    double peakStartTime=(startTDC+peakStart)*digitizationPeriod;
    double peakEndTime=(startTDC+peakEnd)*digitizationPeriod;
    double peakTime=0.5*(peakStartTime+peakEndTime);

    _f1.SetParameter(0, g.GetPointY(peakStart)*TMath::E());
    _f1.SetParameter(1, peakTime);
    _f1.SetParameter(2, _defaultBeta);
    _f1.SetParLimits(0, g.GetPointY(peakStart)*TMath::E()*_minPulseHeightRatio,g.GetPointY(peakStart)*TMath::E()*_maxPulseHeightRatio);
    _f1.SetParLimits(1, peakStartTime-_maxTimeDifference,peakEndTime+_maxTimeDifference);
    _f1.SetParLimits(2, _minBeta, _maxBeta);

    double fitStartTime=(startTDC+fitStart)*digitizationPeriod;
    double fitEndTime=(startTDC+fitEnd)*digitizationPeriod;
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
    float  pulseFitChi2 = (fr->Ndf()>0?fr->Chi2()/fr->Ndf():-1);
    bool   zeroNdf      = (fr->Ndf()>0?false:true);
    bool   failedFit    = FailedFit(fr);

    if(failedFit)
    {
      PEs          = g.GetPointY(peakStart) * TMath::E() * _defaultBeta / calibrationFactor;
      pulseTime    = peakTime;
      pulseHeight  = g.GetPointY(peakStart);
      pulseBeta    = _defaultBeta;
      pulseFitChi2 = -1;
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
    _zeroNdf.push_back(zeroNdf);
    //if one pulse reconstruction fails, all subsequent pulses should also get the error flag,
    //because a wrong pulse subtraction will create wrong secondary pulses
    if(_failedFits.size()>0)
    {
      if(_failedFits.back()) _failedFits.push_back(true);
      else _failedFits.push_back(failedFit);
    }
    else _failedFits.push_back(failedFit);

    SubtractPulse(g, startTDC, digitizationPeriod);
    start=fitEnd;
  }

}

}

