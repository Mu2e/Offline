#include "MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses(double pulseThreshold,       //in V
                                     double leadingEdgeThreshold, //in percent
                                     double param0, 
                                     double param1) : 
                                     _pulseThreshold(pulseThreshold), 
                                     _leadingEdgeThreshold(leadingEdgeThreshold), 
                                     _param0(param0),
                                     _param1(param1)
{}

void MakeCrvRecoPulses::SetWaveform(const std::vector<double> &waveform, double startTime, double binWidth)
{
  _PEs.clear();
  _leadingEdges.clear();
  _pulseHeights.clear();
  _pulseHeightsLandau.clear();
  _peakTimes.clear();
  _peakTimesLandau.clear();
  _integrals.clear();
  _landauParams0.clear();
  _landauParams1.clear();
  _landauParams2.clear();
  _T1s.clear();
  _T2s.clear();
  _TOTs.clear();

  unsigned int nBins = waveform.size();
  double time = startTime;
  for(unsigned bin=0; bin<nBins; bin++, time+=binWidth)
  {
    double voltage = waveform[bin];
    if(voltage>_pulseThreshold)
    {
      TGraph g, gFull;
      double T1 = time;  //start of Landau fit (one bin before 1st bin over threshold)
      double T2 = NAN;   //end of Landau fit (one bin after 1st maximum)
      double TOTstart = time; //start of pulse (1st bin over threshold)
      double TOTend = NAN;
      double integral = 0;
      double maxVoltage = NAN;  //global maximum
      double maxVoltageLandau = NAN; //maximum for Landau fit (1st peak)
      double peakTime = NAN;  //global maximum
      double peakTimeLandau = NAN; //maximum for Landau fit (1st peak)
      bool   insideFitInterval = true;

      //add one more fit point before the waveform crosses the threshold
      if(bin>0) 
      {
        g.SetPoint(0,time-binWidth,waveform[bin-1]); 
        gFull.SetPoint(0,time-binWidth,waveform[bin-1]); 
        T1-=binWidth;
      }

      //loop over all points over the threshold
      for( ; bin<nBins; bin++, time+=binWidth)
      {
        voltage = waveform[bin];
        gFull.SetPoint(g.GetN(),time,voltage);
        if(voltage<_pulseThreshold) break;

        integral += voltage;
        if(insideFitInterval) g.SetPoint(g.GetN(),time,voltage);
        if(voltage>maxVoltage || isnan(maxVoltage)) 
        {
          maxVoltage=voltage;
          peakTime=time;
        }
        else
        {
          insideFitInterval=false;  //the first local maximum (and possibly the global maximum) of this pulse 
          if(isnan(T2)) T2=time;    //has been reached; don't include subsequent points in fit
                                    //subsequent points of this pulse are still used to determine the integral
                                    //and the pulse height
        }
      }
      TOTend=time-binWidth;
      if(isnan(T2)) T2=time;    //if T2 hasn't been found yet, take the last time
                                //this can happen e.g. if the voltage drops below the pulse threshold
                                //right after the maximum

      double leadingEdge = T1;

      TF1 f("","landau");
      TFitResultPtr fr = g.Fit(&f,"NQS");
      int fitStatus = fr;
      if(fitStatus==0)
      {
//if fit is successfull, walk through the Landau graph to find where it becomes 
//20% (_leadingEdgeThreshold) of the maximum.
//this will replace the leading edge and max voltage
        double param0 = fr->Parameter(0);
        double param1 = fr->Parameter(1);
        double param2 = fr->Parameter(2);
        _landauParams0.push_back(param0);
        _landauParams1.push_back(param1);
        _landauParams2.push_back(param2);
        bool foundLeadingEdge=false;
        for(double t=T1; t<T2; t+=1.0)
        {
          double v = param0*TMath::Landau(t, param1, param2);
          if(v>=maxVoltageLandau || isnan(maxVoltageLandau)) 
          {
            maxVoltageLandau=v;
            peakTimeLandau=t;
          }
          if(v>=_leadingEdgeThreshold*maxVoltageLandau && !foundLeadingEdge)
          {
            leadingEdge=t;
            foundLeadingEdge=true;
          }
        }
      }
      else
      {
        _landauParams0.push_back(NAN);
        _landauParams1.push_back(NAN);
        _landauParams2.push_back(NAN);
      }  

//      int PEs = static_cast<int>(integral*_param1+_param0+0.5);  //the 0.5 is used to properly round the doubles
//      int PEs = static_cast<int>(maxVoltageLandau*_param1+_param0+0.5);  //the 0.5 is used to properly round the doubles
      int PEs = static_cast<int>(maxVoltage*_param1+_param0+0.5);  //the 0.5 is used to properly round the doubles  <-- that's what's currently done with the test beam data

      _PEs.push_back(PEs);
      _leadingEdges.push_back(leadingEdge);
      _pulseHeightsLandau.push_back(maxVoltageLandau);
      _pulseHeights.push_back(maxVoltage);
      _peakTimesLandau.push_back(peakTimeLandau);
      _peakTimes.push_back(peakTime);
      _integrals.push_back(integral);
      _T1s.push_back(T1);
      _T2s.push_back(T2);
      _TOTs.push_back(TOTend-TOTstart);
    }
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

double MakeCrvRecoPulses::GetLeadingEdge(int pulse)
{
  int n = _leadingEdges.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _leadingEdges[pulse];
}

double MakeCrvRecoPulses::GetPulseHeight(int pulse)
{
  int n = _pulseHeights.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseHeights[pulse];
}

double MakeCrvRecoPulses::GetPulseHeightLandau(int pulse)
{
  int n = _pulseHeightsLandau.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseHeightsLandau[pulse];
}

double MakeCrvRecoPulses::GetPeakTime(int pulse)
{
  int n = _peakTimes.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _peakTimes[pulse];
}

double MakeCrvRecoPulses::GetPeakTimeLandau(int pulse)
{
  int n = _peakTimesLandau.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _peakTimesLandau[pulse];
}

double MakeCrvRecoPulses::GetIntegral(int pulse)
{
  int n = _integrals.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _integrals[pulse];
}

double MakeCrvRecoPulses::GetLandauParam0(int pulse)
{
  int n = _landauParams0.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _landauParams0[pulse];
}

double MakeCrvRecoPulses::GetLandauParam1(int pulse)
{
  int n = _landauParams1.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _landauParams1[pulse];
}

double MakeCrvRecoPulses::GetLandauParam2(int pulse)
{
  int n = _landauParams2.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _landauParams2[pulse];
}

double MakeCrvRecoPulses::GetT1(int pulse)
{
  int n = _T1s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _T1s[pulse];
}

double MakeCrvRecoPulses::GetT2(int pulse)
{
  int n = _T2s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _T2s[pulse];
}

double MakeCrvRecoPulses::GetTimeOverThreshold(int pulse)
{
  int n = _TOTs.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _TOTs[pulse];
}

}
