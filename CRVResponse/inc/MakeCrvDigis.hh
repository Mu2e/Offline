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

    void SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal, double startTime, double digitizationPrecision, int minADC, int maxADC);

    std::vector<int16_t> GetADCs() {return _ADCs;}
    const std::vector<int16_t> &GetADCs() const {return _ADCs;}

    uint16_t GetTDC() {return _TDC;}

  private:
    std::vector<int16_t> _ADCs;
    uint16_t _TDC;
};

}

#endif
