#include "CrvWaveformResponse.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>

void CrvWaveformResponse::LoadSinglePEWaveform(const std::string &filename, double binWidth, unsigned int nBins) 
{
  _singlePEbinWidth = binWidth;
  std::ifstream f(filename.c_str());
  double currentTime=0, currentVoltage=0;
  double previousTime=NAN, previousVoltage=NAN;
  unsigned int bin=0;
  while(f >> currentTime >> currentVoltage)
  {
    currentTime*=1e9; //convert s into ns, since the text file is in s.
    if(!isnan(previousTime))
    {
      double t=bin*binWidth;  
      while(currentTime>=t && bin<nBins)  //there can be several bins between two time entries in the text file
      {
        double fraction=(t-previousTime)/(currentTime-previousTime);
        double voltage=(currentVoltage-previousVoltage)*fraction+previousVoltage;
        _waveformSinglePE.push_back(voltage);
        bin++;
        t=bin*binWidth;  
      }
      if(bin>=nBins) break;
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

void CrvWaveformResponse::MakeWaveforms(const std::vector<double> &arrivalTimes, 
                                        std::vector<double> &waveform,
                                        double startTime, double binWidth) 
{
  waveform.clear();

  if(arrivalTimes.size()==0) return;

  std::vector<double>::const_iterator iter;
  for(iter=arrivalTimes.begin(); iter!=arrivalTimes.end(); iter++)
  {
    double arrivalTime=*iter;
    unsigned int waveformIndex=ceil((arrivalTime-startTime)/binWidth);  //first available time index
    double digiTime = waveformIndex*binWidth + startTime;  //the time of this time index

    for(; ; waveformIndex++, digiTime+=binWidth)
    {
      double singlePEtime = digiTime - arrivalTime;
      unsigned int singlePEwaveformIndex=static_cast<unsigned int>(singlePEtime/_singlePEbinWidth + 0.5);
      if(singlePEwaveformIndex>=_waveformSinglePE.size()) break; 

      if(waveform.size()<waveformIndex+1) waveform.resize(waveformIndex+1);
      waveform[waveformIndex]+=_waveformSinglePE[singlePEwaveformIndex];
    }
  }
}

