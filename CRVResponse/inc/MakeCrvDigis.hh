#ifndef MakeCrvDigis_h
#define MakeCrvDigis_h

#include <string>
#include <vector>
#include "CLHEP/Random/Randomize.h"

namespace mu2eCrv
{

class MakeCrvDigis
{
  public:
    MakeCrvDigis() {}
    ~MakeCrvDigis() {}
   
    void SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal, double startTime, double digitizationPrecision);

    std::vector<unsigned int> GetADCs() {return _ADCs;}
    const std::vector<unsigned int> &GetADCs() const {return _ADCs;}

    unsigned int GetTDC() {return _TDC;}

  private:
    std::vector<unsigned int> _ADCs;
    unsigned int _TDC;
};

}

#endif
