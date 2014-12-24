#ifndef CrvWaveformResponse_h
#define CrvWaveformResponse_h

#include <string>
#include <vector>

class CrvWaveformResponse
{
  public:
    CrvWaveformResponse() {}
    ~CrvWaveformResponse() {}
   
//the bin width which is used when the single PE waveform is stored must be more precise than
//the bin width of the final waveform
    void LoadSinglePEWaveform(const std::string &filename, double binWidth, unsigned int nBins);
    void MakeWaveforms(const std::vector<double> &arrivalTimes, 
                       std::vector<double> &waveform,
                       double startTime, double binWidth);

  private:
    std::vector<double> _waveformSinglePE;
    double _singlePEbinWidth;
};

#endif
