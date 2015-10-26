#ifndef MakeCrvWaveforms_h
#define MakeCrvWaveforms_h

#include <string>
#include <vector>

namespace mu2eCrv
{

class MakeCrvWaveforms
{
  public:
    MakeCrvWaveforms() {}
    ~MakeCrvWaveforms() {}
   
//the bin width which is used when the single PE waveform is stored must be more precise than
//the bin width of the final waveform
    void LoadSinglePEWaveform(const std::string &filename, double binWidth, unsigned int nBins);
    void MakeWaveform(const std::vector<double> &times, 
                      const std::vector<double> &charges, 
                      std::vector<double> &waveform,
                      double startTime, double binWidth, double timeShift);

  private:
    std::vector<double> _waveformSinglePE;
    double _singlePEbinWidth;
};

}

#endif
