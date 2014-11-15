#ifndef MCDataProducts_CRVPEs_hh
#define MCDataProducts_CRVPEs_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//


namespace mu2e 
{
  class CRVPEs
  {
    public:

    CRVPEs() {}

    std::vector<double> &GetPEtimes(int fiberNumber, int side) 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _PEtimes[SiPM];
    }

    std::vector<double> &GetPEtimes(int SiPMNumber) 
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _PEtimes[SiPMNumber];
    }

    const std::vector<double> &GetPEtimes(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _PEtimes[SiPM];
    }

    const std::vector<double> &GetPEtimes(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _PEtimes[SiPMNumber];
    }

    unsigned int GetNumberOfPEs(int fiberNumber, int side) const 
    {
      if(fiberNumber<0 || fiberNumber>1) throw std::logic_error("Wrong CRV fiber number.");
      if(side<0 || side>1) throw std::logic_error("Wrong CRV side.");
      int SiPM = 2*fiberNumber + side;
      return _PEtimes[SiPM].size();
    }

    unsigned int GetNumberOfPEs(int SiPMNumber) const
    {
      if(SiPMNumber<0 || SiPMNumber>3) throw std::logic_error("Wrong CRV SiPM number.");
      return _PEtimes[SiPMNumber].size();
    }

    double GetFirstPEtime() const
    {
      double firstTime = NAN;
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        if(_PEtimes[SiPM].size()==0) continue;
        double t = *std::min_element(_PEtimes[SiPM].begin(),_PEtimes[SiPM].end());
        if(isnan(firstTime) || t<firstTime) firstTime=t;
      }
      return firstTime;
    }

    private:

    std::vector<double> _PEtimes[4];
  };
}

#endif /* MCDataProducts_CRVPEs_hh */
