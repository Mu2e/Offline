#ifndef RecoDataProducts_CrvRecoPulse_hh
#define RecoDataProducts_CrvRecoPulse_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include <vector>

namespace mu2e 
{
  class CrvRecoPulse
  {
    public:

    CrvRecoPulse() : _PEs(0), _pulseTime(0), _pulseHeight(0), _pulseWidth(0), _pulseFitChi2(0), _LEtime(0), 
                     _scintillatorBarIndex(0), _SiPMNumber(0) {}

    CrvRecoPulse(int PEs, double pulseTime, double pulseHeight, double pulseWidth, double pulseFitChi2, double LEtime, 
                 const std::vector<size_t> &waveformIndices, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) : 
                                                                            _PEs(PEs), 
                                                                            _pulseTime(pulseTime), 
                                                                            _pulseHeight(pulseHeight),
                                                                            _pulseWidth(pulseWidth),
                                                                            _pulseFitChi2(pulseFitChi2),
                                                                            _LEtime(LEtime),
                                                                            _waveformIndices(waveformIndices),
                                                                            _scintillatorBarIndex(scintillatorBarIndex),
                                                                            _SiPMNumber(SiPMNumber)
                                                                             {}

    int    GetPEs() const          {return _PEs;}
    double GetPulseTime() const    {return _pulseTime;}
    double GetPulseHeight() const  {return _pulseHeight;}
    double GetPulseWidth() const   {return _pulseWidth;}
    double GetPulseFitChi2() const {return _pulseFitChi2;}
    double GetLEtime() const       {return _LEtime;}

    const std::vector<size_t>    &GetWaveformIndices() const      {return _waveformIndices;}
    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    private:

    int    _PEs;
    double _pulseTime;
    double _pulseHeight;
    double _pulseWidth;
    double _pulseFitChi2;
    double _LEtime;

    std::vector<size_t>            _waveformIndices;  //indices in the vector of the CrvDigiCollection (which is the same as the index in the CrvDigiMCCollection)
    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
  };
}

#endif /* RecoDataProducts_CrvRecoPulse_hh */
