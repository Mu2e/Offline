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
   
    void SetWaveform(const std::vector<double> &waveform, double ADCconversionFactor, int pedestal);

    std::vector<int> GetADCs() {return _ADCs;}
    const std::vector<int> &GetADCs() const {return _ADCs;}

  private:
    std::vector<int> _ADCs;
};

}

#endif
