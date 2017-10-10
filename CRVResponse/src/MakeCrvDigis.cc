#include "MakeCrvDigis.hh"

namespace mu2eCrv
{

void MakeCrvDigis::SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal)
{
  _ADCs.clear();
  for(size_t i=0; i<waveform.size(); i++)
  {
    _ADCs.push_back(waveform[i]*ADCconversionFactor+pedestal);
  }
}

}
