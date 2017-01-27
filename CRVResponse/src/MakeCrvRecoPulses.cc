#include "MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses(double scale, bool useFittedPulseHeight, bool useFittedPulseTime) : 
                                     _scale(scale),             //conversion between pulse height and PE
                                     _useFittedPulseHeight(useFittedPulseHeight),
                                     _useFittedPulseTime(useFittedPulseTime)
{}

void MakeCrvRecoPulses::SetWaveform(const std::vector<double> &waveform, double startTime, double binWidth)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseWidths.clear();
  _fitParams0.clear();
  _fitParams1.clear();
  _fitParams2.clear();
  _t1s.clear();
  _t2s.clear();
  _PEs.clear();

  unsigned int nBins = waveform.size();
  bool increasing=true;
  int pulseStartBin=0;
  int pulseMaxBin=0;
  int pulseEndBin=0;
  double prevVoltage=0;
  bool doFit=false;
  for(unsigned bin=0; bin<nBins; bin++)
  {
    double voltage = waveform[bin];
    if(increasing)
    {
      if(voltage<prevVoltage && bin>0) //found maximum
      {
        increasing=false;
        pulseMaxBin=bin-1;
      }
      if(bin==nBins-1)  //reached last bin while still increasing, use this last bin as maximum
      {
        pulseMaxBin=bin;
        pulseEndBin=bin;
        doFit=true;
      }
    }
    else
    {
      if(voltage>prevVoltage) //found minimum = border between two pulses
      {
        increasing=true;
        pulseEndBin=bin-1;
        doFit=true;
      }
      if(bin==nBins-1) //reach last bin while decreasing
      {
        pulseEndBin=bin;
        doFit=true;
      }
    }

    if(doFit)
    {
      //find pulse parameters without fit. 
      //will be used, if fit does not succeed or is not possible.
      double pulseTime=startTime+pulseMaxBin*binWidth;
      double pulseHeight=waveform[pulseMaxBin];
      double pulseWidth=binWidth;
      double fitParam0=NAN;
      double fitParam1=NAN;
      double fitParam2=NAN;
      double t1=pulseTime;
      double t2=pulseTime;

      if(pulseStartBin!=pulseMaxBin && pulseEndBin!=pulseMaxBin &&  //pulse maximum is not at the ends of the waveform
         (_useFittedPulseHeight || _useFittedPulseTime))   //only do a fit, if at least one of the fitted values is needed
      {
        int fitStartBin=pulseMaxBin-2;
        int fitEndBin=pulseMaxBin+2;
        if(fitStartBin<pulseStartBin) fitStartBin=pulseStartBin;
        if(fitEndBin>pulseEndBin) fitEndBin=pulseEndBin;

        t1=startTime+fitStartBin*binWidth;
        t2=startTime+fitEndBin*binWidth;

        //fill the graph
        TGraph g;
        for(int i=fitStartBin; i<fitEndBin; i++) 
        {
          double t=startTime+i*binWidth;
          double v=waveform[i];
          g.SetPoint(g.GetN(), t, v);
        }

        //set the fit function
        TF1 f("peakfinder","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))");
        f.SetParameter(0, pulseHeight*2.718);
        f.SetParameter(1, pulseTime);
        f.SetParameter(2, 8);  //TODO: need to check

        //do the fit
        TFitResultPtr fr = g.Fit(&f,"NQS");
        if(fr->IsValid())
        {
          fitParam0 = fr->Parameter(0);
          fitParam1 = fr->Parameter(1);
          fitParam2 = fr->Parameter(2);
          if(_useFittedPulseTime) pulseTime   = fitParam1;
          if(_useFittedPulseHeight) pulseHeight = fitParam0/2.718;    //=fitParam0/e
          pulseWidth  = fitParam2*1.283;    //=fitParam2*pi/sqrt(6)  // =standard deviation of the Gumpel distribution
        }
      }    

      _pulseTimes.push_back(pulseTime);
      _pulseHeights.push_back(pulseHeight);
      _pulseWidths.push_back(pulseWidth);
      _fitParams0.push_back(fitParam0);
      _fitParams1.push_back(fitParam1);
      _fitParams2.push_back(fitParam2);
      _t1s.push_back(t1);
      _t2s.push_back(t2);

      int PE = static_cast<int>(pulseHeight*_scale+0.5);  //the 0.5 is used to properly round the doubles
      _PEs.push_back(PE);

      doFit=false;
      pulseStartBin=pulseEndBin;
    }

    prevVoltage=voltage;
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

double MakeCrvRecoPulses::GetPulseWidth(int pulse)
{
  int n = _pulseWidths.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseWidths[pulse];
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

}
