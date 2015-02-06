#include "MakeCrvWaveforms.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>

void MakeCrvWaveforms::LoadSinglePEWaveform(const std::string &filename, double binWidth, unsigned int nBins) 
{
  _singlePEbinWidth = binWidth;
  std::ifstream f(filename.c_str());
  double currentTime=0, currentVoltage=0;
  double previousTime=NAN, previousVoltage=NAN;
  unsigned int bin=0;
  while(f >> currentTime >> currentVoltage)
  {
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

void MakeCrvWaveforms::MakeWaveform(const std::vector<double> &times, 
                                    const std::vector<double> &charges, 
                                    std::vector<double> &waveform,
                                    double startTime, double binWidth) 
{
  waveform.clear();

  if(times.size()==0) return;

  std::vector<double>::const_iterator iterTime=times.begin();
  std::vector<double>::const_iterator iterCharge=charges.begin();
  for(; iterTime!=times.end() && iterCharge!=charges.end(); iterTime++, iterCharge++)
  {
    double time=*iterTime;
    double charge=*iterCharge;
    unsigned int waveformIndex=ceil((time-startTime)/binWidth);  //first available time index
    double digiTime = waveformIndex*binWidth + startTime;  //the time of this time index

    for(; ; waveformIndex++, digiTime+=binWidth)
    {
      double singlePEtime = digiTime - time;
      unsigned int singlePEwaveformIndex=static_cast<unsigned int>(singlePEtime/_singlePEbinWidth + 0.5);
      if(singlePEwaveformIndex>=_waveformSinglePE.size()) break; 

      if(waveform.size()<waveformIndex+1) waveform.resize(waveformIndex+1);
      waveform[waveformIndex]+=_waveformSinglePE[singlePEwaveformIndex]*charge; 
    }
  }
}

