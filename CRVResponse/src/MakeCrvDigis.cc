#include "Offline/CRVResponse/inc/MakeCrvDigis.hh"

namespace mu2eCrv
{

void MakeCrvDigis::SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal, double startTime, double digitizationPrecision, int minADC, int maxADC)
{
  _ADCs.clear();
  _ADCs.resize(waveform.size());
  for(size_t i=0; i<waveform.size(); i++)
  {
    // clamp before narrowing: an extreme sample cast to int16_t first would wrap around and
    // land back inside [minADC,maxADC], so the saturation clamp would never see it
    double ADCd = waveform[i]*ADCconversionFactor+pedestal+0.5;
    if(ADCd<minADC) ADCd=minADC;
    if(ADCd>maxADC) ADCd=maxADC;
    _ADCs.at(i)=static_cast<int16_t>(ADCd);
  }

  int TDCtmp=lrint(startTime/digitizationPrecision);
  if(TDCtmp<0) throw std::logic_error("ERROR: found a waveform start time (relative to the event marker) < 0"); //this shouldn't happen
  _TDC=static_cast<uint16_t>(TDCtmp);
}

}
