#include "Offline/CRVResponse/inc/MakeCrvWaveforms.hh"

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
  _singlePEWaveformMaxTime = singlePEWaveformMaxTime;
  _singlePEReferenceCharge = singlePEReferenceCharge;
  std::ifstream f(filename.c_str());
  if(!f.good()) throw std::logic_error("Could not open single PE waveform file. "+filename);

  double currentTime=0, currentVoltage=0;
  double previousTime=0, previousVoltage=0;
  unsigned int index=0;
  while(f >> currentTime >> currentVoltage)
  {
    currentTime*=singlePEWaveformStretchFactor;
    if(index!=0)
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

  _singlePEMaxVoltage = *std::max_element(_singlePEWaveform.begin(), _singlePEWaveform.end());
}

void MakeCrvWaveforms::MakeWaveform(const std::vector<std::pair<double,double> > &timesAndCharges,
                                    std::vector<double> &waveform,
                                    double startTime, double digitizationPrecision)
{
  waveform.clear();

  if(timesAndCharges.size()==0) return;
  size_t estimatedNumberOfSamples=(timesAndCharges.back().first-timesAndCharges.front().first+_singlePEWaveformMaxTime)/digitizationPrecision;
  waveform.resize(estimatedNumberOfSamples);

  std::vector<std::pair<double,double> >::const_iterator iterTimesAndCharges=timesAndCharges.begin();
  for(; iterTimesAndCharges!=timesAndCharges.end(); ++iterTimesAndCharges)
  {
    double timeOfCharge=iterTimesAndCharges->first;  //the time when the charge happened
    double charge=iterTimesAndCharges->second/_singlePEReferenceCharge;  //scale it to the 1PE reference charge used for the single PE waveform
    double waveformIndexTmp = ceil((timeOfCharge-startTime)/digitizationPrecision);
    if(waveformIndexTmp<0) waveformIndexTmp=0;
    size_t waveformIndex = static_cast<size_t>(lrint(waveformIndexTmp));  //waveform index of the first digitization point for this particular charge
    double waveformTime = waveformIndex*digitizationPrecision + startTime;  //the time for this waveform index

    for(; ; waveformIndex++, waveformTime+=digitizationPrecision)
    {
      double singlePEWaveformTime = waveformTime - timeOfCharge;
      if(singlePEWaveformTime<0) continue;
      size_t singlePEwaveformIndex=static_cast<size_t>(lrint(singlePEWaveformTime/_singlePEWaveformPrecision));
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
