#ifndef MCDataProducts_CrvWaveforms_hh
#define MCDataProducts_CrvWaveforms_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CrvWaveforms
  {
    public:

    CrvWaveforms() {}

    std::vector<double> &GetWaveform(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _waveform[SiPM];
    }

    std::vector<double> &GetWaveform(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _waveform[SiPMNumber];
    }

    const std::vector<double> &GetWaveform(int fiberNumber, int side) const
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _waveform[SiPM];
    }

    const std::vector<double> &GetWaveform(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _waveform[SiPMNumber];
    }

    double GetStartTime(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _startTime[SiPM];
    }

    double GetStartTime(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _startTime[SiPMNumber];
    }

    void SetStartTime(int fiberNumber, int side, double startTime) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      _startTime[SiPM]=startTime;
    }

    void SetStartTime(int SiPMNumber, double startTime) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      _startTime[SiPMNumber]=startTime;
    }

    double GetBinWidth() const 
    {
      return _binWidth;
    }

    void SetBinWidth(double binWidth) 
    {
      _binWidth=binWidth;
    }

    private:

    std::vector<double> _waveform[4];
    double _startTime[4];
    double _binWidth;
  };
}

#endif /* MCDataProducts_CrvWaveforms_hh */
