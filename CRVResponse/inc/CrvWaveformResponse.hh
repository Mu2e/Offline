#ifndef CrvWaveformResponse_h
#define CrvWaveformResponse_h

#include <string>
#include <vector>

class CrvWaveformResponse
{
  public:
    CrvWaveformResponse() {}
    ~CrvWaveformResponse() {}
   
    void LoadSinglePEWaveform(const std::string &filename, double binWidth, double maxTime);
    void makeWaveforms(const std::vector<double> &arrivalTimes, 
                       std::vector<double> &waveform,
                       double &startTime);

  private:
    double _binWidth;
    std::vector<double> _waveformSinglePE;
};

#endif
