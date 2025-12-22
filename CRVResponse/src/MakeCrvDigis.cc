#include "Offline/CRVResponse/inc/MakeCrvDigis.hh"

namespace mu2eCrv
{

void MakeCrvDigis::SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal, double startTime, double digitizationPrecision, int minADC, int maxADC)
{
  _ADCs.clear();
  _ADCs.resize(waveform.size());
  for(size_t i=0; i<waveform.size(); i++)
  {
    int16_t ADC = static_cast<int16_t>(waveform[i]*ADCconversionFactor+pedestal+0.5);
    if(ADC<minADC) ADC=minADC;
    if(ADC>maxADC) ADC=maxADC;
    _ADCs.at(i)=ADC;
  }

  int TDCtmp=lrint(startTime/digitizationPrecision);
  if(TDCtmp<0) throw std::logic_error("ERROR: found a waveform start time (relative to the event marker) < 0"); //this shouldn't happen
  _TDC=static_cast<uint16_t>(TDCtmp);
}

}
