#ifndef RecoDataProducts_CrvRecoPulses_hh
#define RecoDataProducts_CrvRecoPulses_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include <vector>
#include <cmath>
#include <cstddef>

namespace mu2e 
{
  class CrvRecoPulses
  {
    public:

    CrvRecoPulses() {}

    struct CrvSingleRecoPulse
    {
      int    _PEs;
      double _pulseTime;
      double _pulseHeight;
      double _pulseWidth;
      double _pulseFitChi2;
      double _LEtime;
      size_t _singleWaveformIndex;   //index in the vector of singleWaveforms in the CrvDigi (which is the same as the index in the CrvDigiMC)
      CrvSingleRecoPulse(int PEs, double pulseTime, double pulseHeight, double pulseWidth, 
                         double pulseFitChi2, double LEtime, size_t singleWaveformIndex) : 
                                                                            _PEs(PEs), 
                                                                            _pulseTime(pulseTime), 
                                                                            _pulseHeight(pulseHeight),
                                                                            _pulseWidth(pulseWidth),
                                                                            _pulseFitChi2(pulseFitChi2),
                                                                            _LEtime(LEtime),
                                                                            _singleWaveformIndex(singleWaveformIndex) {}
      CrvSingleRecoPulse() : _PEs(0), _pulseTime(NAN), _pulseHeight(NAN), _pulseWidth(NAN), 
                             _pulseFitChi2(NAN), _LEtime(NAN), _singleWaveformIndex(0) {}  //to make ROOT happy
    };

    std::vector<CrvSingleRecoPulse> &GetRecoPulses(int fiberNumber, int side);
    std::vector<CrvSingleRecoPulse> &GetRecoPulses(int SiPMNumber); 

    const std::vector<CrvSingleRecoPulse> &GetRecoPulses(int fiberNumber, int side) const;
    const std::vector<CrvSingleRecoPulse> &GetRecoPulses(int SiPMNumber) const;

    unsigned int GetNumberOfRecoPulses(int fiberNumber, int side) const;
    unsigned int GetNumberOfRecoPulses(int SiPMNumber) const;

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<CrvSingleRecoPulse> _crvPulses[4];
  };
}

#endif /* RecoDataProducts_CrvRecoPulses_hh */
