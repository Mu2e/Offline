#include "CRVResponse/inc/MakeCrvWaveforms.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace mu2eCrv
{

void MakeCrvWaveforms::LoadSinglePEWaveform(const std::string &filename, double singlePEWaveformPrecision, double singlePEWaveformStretchFactor,
                                            double singlePEWaveformMaxTime, double singlePEReferenceCharge) 
{
  _singlePEWaveformPrecision = singlePEWaveformPrecision;
  _singlePEReferenceCharge = singlePEReferenceCharge;
  std::ifstream f(filename.c_str());
  if(!f.good()) throw std::logic_error("Could not open single PE waveform file. "+filename);

  double currentTime=0, currentVoltage=0;
  double previousTime=NAN, previousVoltage=NAN;
  unsigned int index=0;
  while(f >> currentTime >> currentVoltage)
  {
    currentTime*=singlePEWaveformStretchFactor;
    if(!std::isnan(previousTime))
    {
      double t=index*singlePEWaveformPrecision;  
      while(currentTime>=t && index*singlePEWaveformPrecision<singlePEWaveformMaxTime)
      {
        double fraction=(t-previousTime)/(currentTime-previousTime);
        double voltage=(currentVoltage-previousVoltage)*fraction+previousVoltage;
        _singlePEWaveform.push_back(voltage);
        index++;
        t=index*singlePEWaveformPrecision;  
      }
      if(index*singlePEWaveformPrecision>=singlePEWaveformMaxTime) break;
    }
    else
    {
      _singlePEWaveform.push_back(currentVoltage);
      index++;
    }
    previousTime=currentTime;
    previousVoltage=currentVoltage;
  }
  f.close();
}

void MakeCrvWaveforms::MakeWaveform(const std::vector<double> &times, 
                                    const std::vector<double> &charges, 
                                    std::vector<double> &waveform,
                                    double startTime, double digitizationPrecision) 
{
  waveform.clear();

  if(times.size()==0) return;

  std::vector<double>::const_iterator iterTime=times.begin();
  std::vector<double>::const_iterator iterCharge=charges.begin();
  for(; iterTime!=times.end() && iterCharge!=charges.end(); iterTime++, iterCharge++)
  {
    double timeOfCharge=*iterTime;  //the time when the charge happened
    double charge=*iterCharge/_singlePEReferenceCharge;  //scale it to the 1PE reference charge used for the single PE waveform
    unsigned int waveformIndex = (timeOfCharge>startTime ? ceil((timeOfCharge-startTime)/digitizationPrecision) : 0);  //waveform index of the first digitization point for this particular charge
    double waveformTime = waveformIndex*digitizationPrecision + startTime;  //the time for this waveform index

    for(; ; waveformIndex++, waveformTime+=digitizationPrecision)
    {
      double singlePEWaveformTime = waveformTime - timeOfCharge;
      if(singlePEWaveformTime<0) continue;  //this shouldn't happen
      unsigned int singlePEwaveformIndex=static_cast<unsigned int>(lrint(singlePEWaveformTime/_singlePEWaveformPrecision));
      if(singlePEwaveformIndex>=_singlePEWaveform.size()) break; 

      if(waveform.size()<waveformIndex+1) waveform.resize(waveformIndex+1,0);  //new vector elements are set to 0
      waveform[waveformIndex]+=_singlePEWaveform[singlePEwaveformIndex]*charge; 
    }
  }
}

void MakeCrvWaveforms::AddElectronicNoise(std::vector<double> &waveform, double noise, CLHEP::RandGaussQ &randGaussQ)
{
  std::vector<double>::iterator iter;
  for(iter=waveform.begin(); iter!=waveform.end(); iter++)
  {
    double n = randGaussQ.fire(0, noise);
    *iter+=n;
  }
}

}
