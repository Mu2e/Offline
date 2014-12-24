#include "CrvRecoPulseResponse.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

CrvRecoPulseResponse::CrvRecoPulseResponse(double pulseThreshold,       //in V
                                           double leadingEdgeThreshold, //in percent
                                           double integralFactor) : 
                                           _pulseThreshold(pulseThreshold), 
                                           _leadingEdgeThreshold(leadingEdgeThreshold), 
                                           _integralFactor(integralFactor)
{}

void CrvRecoPulseResponse::SetWaveform(const std::vector<double> &waveform, double startTime, double binWidth)
{
  _PEs.clear();
  _leadingEdges.clear();
  _pulseHeights.clear();

  _integrals.clear();
  _landauParams0.clear();
  _landauParams1.clear();
  _landauParams2.clear();
  _T1s.clear();
  _T2s.clear();

  unsigned int nBins = waveform.size();
  double time = startTime;
  for(unsigned bin=0; bin<nBins; bin++, time+=binWidth)
  {
    double voltage = waveform[bin];
    if(voltage>_pulseThreshold)
    {
      TGraph g;
      double T1 = time;
      double T2 = NAN;
      double integral = 0;
      double maxVoltage = NAN;
      bool   insideFitInterval = true;

      //add one more fit point before the waveform crosses the threshold
      if(bin>0) {g.SetPoint(0,time-binWidth,waveform[bin-1]); T1-=binWidth;}

      //loop over all points over the threshold
      for( ; bin<nBins; bin++, time+=binWidth)
      {
        voltage = waveform[bin];
        if(voltage<_pulseThreshold) break;

        if(insideFitInterval) g.SetPoint(g.GetN(),time,voltage);
        integral += voltage;
        if(voltage>maxVoltage || isnan(maxVoltage)) maxVoltage=voltage;
        else
        {
          insideFitInterval=false;  //a local maximum (and possibly the global maximum) has been reached
          if(isnan(T2)) T2=time;    //don't include subsequent points in fit
        }
      }
      if(isnan(T2)) T2=time;    //if T2 hasn't been found yet, take the last time
                                //this can happen e.g. if the voltage drops below the pulse threshold
                                //right after the maximum

      double PEs = static_cast<int>(integral*_integralFactor+0.5);
      double leadingEdge = T1;

      TF1 f("","landau");
      TFitResultPtr fr = g.Fit(&f,"NQS");
      int fitStatus = fr;
      if(fitStatus==0)
      {
//if fit is successfull, walk through the Landau graph to find where it becomes 
//20% (_leadingEdgeThreshold) of the maximum.
//this will replace the leading edge
        double param0 = fr->Parameter(0);
        double param1 = fr->Parameter(1);
        double param2 = fr->Parameter(2);
        _landauParams0.push_back(param0);
        _landauParams1.push_back(param1);
        _landauParams2.push_back(param2);
        for(double t=T1; t<T2; t+=1.0)
        {
          double v = param0*TMath::Landau(t, param1, param2);
          if(v>=_leadingEdgeThreshold*maxVoltage)
          {
            leadingEdge=t;
            break;
          }
        }
      }
      else
      {
        _landauParams0.push_back(NAN);
        _landauParams1.push_back(NAN);
        _landauParams2.push_back(NAN);
      }  

      _PEs.push_back(PEs);
      _leadingEdges.push_back(leadingEdge);
      _pulseHeights.push_back(maxVoltage);
      _integrals.push_back(integral);
      _T1s.push_back(T1);
      _T2s.push_back(T2);
    }
  }
}

unsigned int CrvRecoPulseResponse::GetNPulses()
{
  return _PEs.size();
}

double CrvRecoPulseResponse::GetPEs(int pulse)
{
  int n = _PEs.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _PEs[pulse];
}

double CrvRecoPulseResponse::GetLeadingEdge(int pulse)
{
  int n = _leadingEdges.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _leadingEdges[pulse];
}

double CrvRecoPulseResponse::GetPulseHeight(int pulse)
{
  int n = _pulseHeights.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseHeights[pulse];
}

double CrvRecoPulseResponse::GetIntegral(int pulse)
{
  int n = _integrals.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _integrals[pulse];
}

double CrvRecoPulseResponse::GetLandauParam0(int pulse)
{
  int n = _landauParams0.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _landauParams0[pulse];
}

double CrvRecoPulseResponse::GetLandauParam1(int pulse)
{
  int n = _landauParams1.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _landauParams1[pulse];
}

double CrvRecoPulseResponse::GetLandauParam2(int pulse)
{
  int n = _landauParams2.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _landauParams2[pulse];
}

double CrvRecoPulseResponse::GetT1(int pulse)
{
  int n = _T1s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _T1s[pulse];
}

double CrvRecoPulseResponse::GetT2(int pulse)
{
  int n = _T2s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _T2s[pulse];
}

