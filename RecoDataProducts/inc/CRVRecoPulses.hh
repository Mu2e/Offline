#ifndef RecoDataProducts_CRVRecoPulses_hh
#define RecoDataProducts_CRVRecoPulses_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CRVRecoPulses
  {
    public:

    CRVRecoPulses() {}

    struct CRVSingleRecoPulse
    {
      int    _PEs;
      double _leadingEdge;
      double _timeOverThreshold;
      double _pulseHeight;
    };

    std::vector<CRVSingleRecoPulse> &GetRecoPulses(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _crvPulses[SiPM];
    }

    std::vector<CRVSingleRecoPulse> &GetRecoPulses(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _crvPulses[SiPMNumber];
    }

    const std::vector<CRVSingleRecoPulse> &GetRecoPulses(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _crvPulses[SiPM];
    }

    const std::vector<CRVSingleRecoPulse> &GetRecoPulses(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _crvPulses[SiPMNumber];
    }

    unsigned int GetNumberOfRecoPulses(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _crvPulses[SiPM].size();
    }

    unsigned int GetNumberOfRecoPulses(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _crvPulses[SiPMNumber].size();
    }

    private:

    std::vector<CRVSingleRecoPulse> _crvPulses[4];
  };
}

#endif /* RecoDataProducts_CRVRecoPulses_hh */
