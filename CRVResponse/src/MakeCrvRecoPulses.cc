#include "CRVResponse/inc/MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
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

MakeCrvRecoPulses::MakeCrvRecoPulses(double minADCdifference, double defaultBeta, double minBeta, double maxBeta,
                                     double maxTimeDifference, double minPulseHeightRatio, double maxPulseHeightRatio,
                                     double LEtimeFactor) :
                                     _f("peakfitter",Gumbel,0,0,3), _minADCdifference(minADCdifference), 
                                     _defaultBeta(defaultBeta), _minBeta(minBeta), _maxBeta(maxBeta), 
                                     _maxTimeDifference(maxTimeDifference),
                                     _minPulseHeightRatio(minPulseHeightRatio),
                                     _maxPulseHeightRatio(maxPulseHeightRatio),
                                     _LEtimeFactor(LEtimeFactor)
{}

void MakeCrvRecoPulses::SetWaveform(const std::vector<unsigned int> &waveform, 
                                    unsigned int startTDC, double digitizationPeriod, double pedestal, 
                                    double calibrationFactor, double calibrationFactorPulseHeight)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseBetas.clear();
  _pulseFitChi2s.clear();
  _fitParams0.clear();
  _fitParams1.clear();
  _fitParams2.clear();
  _t1s.clear();
  _t2s.clear();
  _PEs.clear();
  _PEsPulseHeight.clear();
  _LEtimes.clear();
  _peakBins.clear();

  //find the maxima
  int nBins = static_cast<int>(waveform.size());
  std::vector<std::pair<int,int> > peaks;
  int peakStartBin=0;
  int peakEndBin=0;
  for(int bin=2; bin<nBins-2; bin++) 
  {
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
        peaks.emplace_back(peakStartBin,peakEndBin);
        peakStartBin=0;  //so that the loop as to wait for the next rising edge
      }
    }
  }

  for(size_t i=0; i<peaks.size(); i++)
  {
  //select a range of up to 4 points before and after the maximum point
  //-find up to 5 points before and after the maximum point for which the waveform is stricly decreasing
  //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)
    peakStartBin=peaks[i].first;
    peakEndBin=peaks[i].second;
    if(waveform[peakStartBin]-pedestal<_minADCdifference) continue;  //ignore peaks which are too small
                                                                     //(even smaller than 1PE dark counts)

    int startBin=peakStartBin;
    int endBin=peakEndBin;
    for(int bin=peakStartBin-1; bin>=0 && bin>=peakStartBin-5; bin--)
    {
      if(waveform[bin]<=waveform[bin+1]) startBin=bin;
      else break;
    }
    for(int bin=peakEndBin+1; bin<nBins && bin<=peakEndBin+5; bin++)
    {
      if(waveform[bin]<=waveform[bin-1]) endBin=bin;
      else break;
    }
    if(peakStartBin-startBin>1) startBin++;
    if(endBin-peakEndBin>1) endBin--;

    double t1=(startTDC+startBin)*digitizationPeriod;
    double t2=(startTDC+endBin)*digitizationPeriod;
    double peakTime=(startTDC+0.5*(peakStartBin+peakEndBin))*digitizationPeriod;

    //fill the graph
    std::vector<double> t,v;
    for(int bin=startBin; bin<=endBin; bin++) 
    {
      t.emplace_back((startTDC+bin)*digitizationPeriod);
      v.emplace_back(waveform[bin]-pedestal);
    }
    TGraph g(t.size(), t.data(), v.data());

    //set the fit function
    _f.SetParameter(0, (waveform[peakStartBin]-pedestal)*TMath::E());
    _f.SetParameter(1, peakTime);
    _f.SetParameter(2, _defaultBeta);

    //do the fit
    TFitResultPtr fr = g.Fit(&_f,"NQS");

    bool invalidFit=false;
    if(!fr->IsValid()) invalidFit=true;;

    double fitParam0 = fr->Parameter(0);
    double fitParam1 = fr->Parameter(1);
    double fitParam2 = fr->Parameter(2);

    //trying to identify fake pulse
    if(fitParam0<=0) invalidFit=true;
    if(fitParam2<_minBeta || fitParam2>_maxBeta) invalidFit=true;
    if(fabs(fitParam1-peakTime)>_maxTimeDifference) invalidFit=true;
    if(fitParam0/((waveform[peakStartBin]-pedestal)*TMath::E())>_maxPulseHeightRatio) invalidFit=true;
    if(fitParam0/((waveform[peakStartBin]-pedestal)*TMath::E())<_minPulseHeightRatio) invalidFit=true;

    int    PEs          = lrint(fitParam0*fitParam2 / calibrationFactor);
    double pulseTime    = fitParam1;
    double pulseHeight  = fitParam0/TMath::E();
    double pulseBeta    = fitParam2;
    double pulseFitChi2 = fr->Chi2();
    double LEtime       = 0;
    if(invalidFit)
    {
      PEs          = lrint((waveform[peakStartBin]-pedestal)*TMath::E() * _defaultBeta / calibrationFactor);
      pulseTime    = peakTime;
      pulseHeight  = waveform[peakStartBin]-pedestal;
      pulseBeta    = NAN;
      pulseFitChi2 = NAN;
      LEtime       = peakTime-_LEtimeFactor*_defaultBeta;  //50% pulse height is reached at -0.985*beta before the peak
    }
    else LEtime    = pulseTime-_LEtimeFactor*pulseBeta;  //50% pulse height is reached at -0.985*beta before the peak

    int  PEsPulseHeight = lrint(pulseHeight / calibrationFactorPulseHeight);

    _pulseTimes.push_back(pulseTime);
    _pulseHeights.push_back(pulseHeight);
    _pulseBetas.push_back(pulseBeta);
    _pulseFitChi2s.push_back(pulseFitChi2);
    _fitParams0.push_back(fitParam0);
    _fitParams1.push_back(fitParam1);
    _fitParams2.push_back(fitParam2);
    _t1s.push_back(t1);
    _t2s.push_back(t2);
    _PEs.push_back(PEs);
    _PEsPulseHeight.push_back(PEsPulseHeight);
    _LEtimes.push_back(LEtime);
    _peakBins.push_back(peakStartBin);
  }
}

unsigned int MakeCrvRecoPulses::GetNPulses()
{
  return _PEs.size();
}

int MakeCrvRecoPulses::GetPEs(int pulse)
{
  int n = _PEs.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _PEs[pulse];
}

int MakeCrvRecoPulses::GetPEsPulseHeight(int pulse)
{
  int n = _PEs.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _PEsPulseHeight[pulse];
}

double MakeCrvRecoPulses::GetPulseTime(int pulse)
{
  int n = _pulseTimes.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseTimes[pulse];
}

double MakeCrvRecoPulses::GetPulseHeight(int pulse)
{
  int n = _pulseHeights.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseHeights[pulse];
}

double MakeCrvRecoPulses::GetPulseBeta(int pulse)
{
  int n = _pulseBetas.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseBetas[pulse];
}

double MakeCrvRecoPulses::GetPulseFitChi2(int pulse)
{
  int n = _pulseFitChi2s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseFitChi2s[pulse];
}

double MakeCrvRecoPulses::GetFitParam0(int pulse)
{
  int n = _fitParams0.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _fitParams0[pulse];
}

double MakeCrvRecoPulses::GetFitParam1(int pulse)
{
  int n = _fitParams1.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _fitParams1[pulse];
}

double MakeCrvRecoPulses::GetFitParam2(int pulse)
{
  int n = _fitParams2.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _fitParams2[pulse];
}

double MakeCrvRecoPulses::GetT1(int pulse)
{
  int n = _t1s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _t1s[pulse];
}

double MakeCrvRecoPulses::GetT2(int pulse)
{
  int n = _t2s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _t2s[pulse];
}

double MakeCrvRecoPulses::GetLEtime(int pulse)
{
  int n = _LEtimes.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _LEtimes[pulse];
}

int MakeCrvRecoPulses::GetPeakBin(int pulse)
{
  int n = _peakBins.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _peakBins[pulse];
}

}

