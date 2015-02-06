#ifndef MCDataProducts_CrvSiPMRespones_hh
#define MCDataProducts_CrvSiPMRespones_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CrvSiPMResponses
  {
    public:

    CrvSiPMResponses() {}

    struct CrvSingleSiPMResponse
    {
      double _time;
      double _charge;
    };

    std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _crvSiPMResponses[SiPM];
    }

    std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _crvSiPMResponses[SiPMNumber];
    }

    const std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _crvSiPMResponses[SiPM];
    }

    const std::vector<CrvSingleSiPMResponse> &GetSiPMResponses(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _crvSiPMResponses[SiPMNumber];
    }

    unsigned int GetNumberOfSiPMResponses(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _crvSiPMResponses[SiPM].size();
    }

    unsigned int GetNumberOfSiPMResponses(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _crvSiPMResponses[SiPMNumber].size();
    }

    double GetFirstSiPMResponseTime() const
    {
      double firstTime = NAN;
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        if(_crvSiPMResponses[SiPM].size()==0) continue;
        double t = _crvSiPMResponses[SiPM].front()._time;
        if(isnan(firstTime) || t<firstTime) firstTime=t;
      }
      return firstTime;
    }

    double GetLastSiPMResponseTime() const
    {
      double lastTime = NAN;
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        if(_crvSiPMResponses[SiPM].size()==0) continue;
        double t = _crvSiPMResponses[SiPM].back()._time;
        if(isnan(lastTime) || t<lastTime) lastTime=t;
      }
      return lastTime;
    }

    private:

    std::vector<CrvSingleSiPMResponse> _crvSiPMResponses[4];
  };
}

#endif /* MCDataProducts_CrvSiPMRespones_hh */
