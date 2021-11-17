#include "Offline/CRVResponse/inc/MakeCrvDigis.hh"

namespace mu2eCrv
{

void MakeCrvDigis::SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal, double startTime, double digitizationPrecision)
{
  _ADCs.clear();
  _ADCs.reserve(waveform.size());
  for(size_t i=0; i<waveform.size(); i++)
  {
    _ADCs.push_back(static_cast<int16_t>(waveform[i]*ADCconversionFactor+pedestal+0.5));
  }

  int TDCtmp=lrint(startTime/digitizationPrecision);
  if(TDCtmp<0) throw std::logic_error("ERROR: found a waveform start time (relative to the event marker) < 0"); //this shouldn't happen
  _TDC=static_cast<uint16_t>(TDCtmp);
}

}
