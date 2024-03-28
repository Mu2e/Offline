#ifndef MakeCrvWaveforms_h
#define MakeCrvWaveforms_h

#include <string>
#include <vector>
#include "CLHEP/Random/Randomize.h"

namespace mu2eCrv
{

class MakeCrvWaveforms
{
  public:
    MakeCrvWaveforms() {}
    ~MakeCrvWaveforms() {}

//the precision which is used when the single PE waveform is stored must be more precise than
//the digitizaionInterval of the final waveform
    void LoadSinglePEWaveform(const std::string &filename, double singlePEWaveformPrecision, double singlePEWaveformStrechFactor,
                              double singlePEWaveformMaxTime, double singlePEReferenceCharge);
    void MakeWaveform(const std::vector<std::pair<double,double> > &timesAndCharges,
                      std::vector<double> &waveform,
                      double startTime, double digitizationInterval);
    void AddElectronicNoise(std::vector<double> &waveform, double noise, CLHEP::RandGaussQ &randGaussQ);
    double GetSinglePEMaxVoltage() {return _singlePEMaxVoltage;}

  private:
    std::vector<double> _singlePEWaveform;
    double _singlePEWaveformPrecision{0.0};
    double _singlePEWaveformMaxTime{0.0};
    double _singlePEReferenceCharge{0.0};
    double _singlePEMaxVoltage{0.0};
};

}

#endif
