#include "CRVResponse/inc/MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses()  
{}

void MakeCrvRecoPulses::SetWaveform(const std::vector<unsigned int> &waveform, unsigned int startTDC, double digitizationPeriod, 
                                    double pedestal, double calibrationFactor, double calibrationFactorPulseHeight, bool darkNoise)
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
  std::vector<std::pair<int,bool> > peaks;
  for(int bin=2; bin<nBins-2; bin++) 
  {
    if(waveform[bin-1]<waveform[bin] && waveform[bin]>waveform[bin+1]) peaks.emplace_back(bin,false);
    if(waveform[bin-1]<waveform[bin] && waveform[bin]==waveform[bin+1] && waveform[bin+1]>waveform[bin+2]) peaks.emplace_back(bin,true);
  }

  for(size_t i=0; i<peaks.size(); i++)
  {
  //select a range of up to 4 points before and after the maximum point
  //-find up to 5 points before and after the maximum point for which the waveform is stricly decreasing
  //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)
    int maxBin = peaks[i].first;
    if(waveform[maxBin]-pedestal<5) continue; //FIXME: need a better way to identify these fake pulse which are caused by electronic noise

    int startBin=maxBin;
    int endBin=maxBin;
    for(int bin=maxBin-1; bin>=0 && bin>=maxBin-5; bin--)
    {
      if(waveform[bin]<=waveform[bin+1]) startBin=bin;
      else break;
    }
    for(int bin=maxBin+1; bin<nBins && bin<=maxBin+5; bin++)
    {
      if(waveform[bin]<=waveform[bin-1]) endBin=bin;
      else break;
    }
    if(maxBin-startBin>1) startBin++;
    if(endBin-maxBin>1) endBin--;

    double t1=(startTDC+startBin)*digitizationPeriod;
    double t2=(startTDC+endBin)*digitizationPeriod;

    //fill the graph
    TGraph g;
    for(int bin=startBin; bin<=endBin; bin++) 
    {
      double t=(startTDC+bin)*digitizationPeriod;
      double v=waveform[bin]-pedestal;
      g.SetPoint(g.GetN(), t, v);
    }

    //set the fit function
    TF1 f("peakfitter","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))");
    f.SetParameter(0, (waveform[maxBin]-pedestal)*TMath::E());
    f.SetParameter(1, (startTDC+maxBin)*digitizationPeriod);
    f.SetParameter(2, darkNoise?12.6:19.0);
    if(peaks[i].second) f.SetParameter(1, (startTDC+maxBin+0.5)*digitizationPeriod);

    //do the fit
    TFitResultPtr fr = g.Fit(&f,"NQS");
    if(!fr->IsValid()) continue;

    double fitParam0 = fr->Parameter(0);
    double fitParam1 = fr->Parameter(1);
    double fitParam2 = fr->Parameter(2);
    if(fitParam0<=0 || fitParam2<=0) continue;
    if(fitParam2>50) continue; //FIXME: need a better way to identify these fake pulse which are caused by electronic noise
    if(fabs(fitParam1-(startTDC+maxBin)*digitizationPeriod)>30) continue; //FIXME
    if(fitParam0/((waveform[maxBin]-pedestal)*TMath::E())>2.0) continue; //FIXME

    int    PEs          = lrint(fitParam0*fitParam2 / calibrationFactor);
    double pulseTime    = fitParam1;
    double pulseHeight  = fitParam0/TMath::E();
    double pulseBeta    = fitParam2;
    double pulseFitChi2 = fr->Chi2();

    double LEtime=f.GetX(0.5*pulseHeight,pulseTime-50,pulseTime);   //i.e. at 50% of pulse height
    int    PEsPulseHeight = lrint(pulseHeight / calibrationFactorPulseHeight);

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
    _peakBins.push_back(maxBin);
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

