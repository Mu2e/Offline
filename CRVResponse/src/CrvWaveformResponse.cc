#include "CrvWaveformResponse.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>

void CrvWaveformResponse::LoadSinglePEWaveform(const std::string &filename, double binWidth, double maxTime) 
{
  _binWidth=binWidth;  //in ns

  std::ifstream f(filename.c_str());
  double currentTime=0, currentVoltage=0;
  double previousTime=NAN, previousVoltage=NAN;
  unsigned int bin=0;
  while(f >> currentTime >> currentVoltage)
  {
    currentTime*=1e9; //convert s into ns, since the text file is in s.
    if(currentTime>maxTime) break;
    if(!isnan(previousTime))
    {
      double t=bin*_binWidth;  
      while(currentTime>=t)  //there can be several bins between two time entries in the text file
      {
        double fraction=(t-previousTime)/(currentTime-previousTime);
        double voltage=(currentVoltage-previousVoltage)*fraction+previousVoltage;
        _waveformSinglePE.push_back(voltage);
        bin++;
        t=bin*_binWidth;  
      }
    }
    else
    {
      _waveformSinglePE.push_back(currentVoltage);
      bin++;
    }
    previousTime=currentTime;
    previousVoltage=currentVoltage;
  }
  f.close();
}

void CrvWaveformResponse::makeWaveforms(const std::vector<double> &arrivalTimes, 
                                        std::vector<double> &waveform,
                                        double &startTime) 
{
  waveform.clear();
  startTime=NAN;

  if(arrivalTimes.size()==0) return;

  startTime=*std::min_element(arrivalTimes.begin(),arrivalTimes.end());

  std::vector<double>::const_iterator iter;
  for(iter=arrivalTimes.begin(); iter!=arrivalTimes.end(); iter++)
  {
    double t=*iter;
    unsigned int waveformIndexStart=static_cast<unsigned int>((t-startTime)/_binWidth);
    for(unsigned int bin=0; bin<_waveformSinglePE.size(); bin++)
    {
      unsigned int waveformIndex=waveformIndexStart+bin;
      if(waveform.size()<waveformIndex+1) waveform.resize(waveformIndex+1);
      waveform[waveformIndex]+=_waveformSinglePE[bin];
    }
  }
}

