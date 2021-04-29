#ifndef RecoDataProducts_CrvRecoPulse_hh
#define RecoDataProducts_CrvRecoPulse_hh
//
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "RecoDataProducts/inc/CrvRecoPulseFlags.hh"

#include <vector>

namespace mu2e 
{
  class CrvRecoPulse
  {
    public:

    CrvRecoPulse() : _PEs(0), _PEsPulseHeight(0), _pulseTime(0), _pulseHeight(0), _pulseBeta(0), _pulseFitChi2(0), _LEtime(0), 
                     _flags(0), _PEsNoFit(0), _pulseTimeNoFit(0), _pulseStart(0), _pulseEnd(0), _scintillatorBarIndex(0), _SiPMNumber(0) {}

    CrvRecoPulse(float PEs, float PEsPulseHeight, double pulseTime, float pulseHeight, float pulseBeta, float pulseFitChi2, double LEtime, 
                 const CrvRecoPulseFlags &flags,
                 float PEsNoFit, double pulseTimeNoFit, double pulseStart, double pulseEnd, 
                 const std::vector<size_t> &waveformIndices, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) : 
                                                                            _PEs(PEs), 
                                                                            _PEsPulseHeight(PEsPulseHeight), 
                                                                            _pulseTime(pulseTime), 
                                                                            _pulseHeight(pulseHeight),
                                                                            _pulseBeta(pulseBeta),
                                                                            _pulseFitChi2(pulseFitChi2),
                                                                            _LEtime(LEtime),
                                                                            _flags(flags),
                                                                            _PEsNoFit(PEsNoFit),
                                                                            _pulseTimeNoFit(pulseTimeNoFit),
                                                                            _pulseStart(pulseStart),
                                                                            _pulseEnd(pulseEnd),
                                                                            _waveformIndices(waveformIndices),
                                                                            _scintillatorBarIndex(scintillatorBarIndex),
                                                                            _SiPMNumber(SiPMNumber)
                                                                             {}

    float  GetPEs() const            {return _PEs;}
    float  GetPEsPulseHeight() const {return _PEsPulseHeight;}
    double GetPulseTime() const      {return _pulseTime;}
    float  GetPulseHeight() const    {return _pulseHeight;}
    float  GetPulseBeta() const      {return _pulseBeta;}
    float  GetPulseFitChi2() const   {return _pulseFitChi2;}
    double GetLEtime() const         {return _LEtime;}
    const  CrvRecoPulseFlags &GetRecoPulseFlags() const {return _flags;}

    float  GetPEsNoFit() const       {return _PEsNoFit;}
    double GetPulseTimeNoFit() const {return _pulseTimeNoFit;}
    double GetPulseStart() const     {return _pulseStart;}
    double GetPulseEnd() const       {return _pulseEnd;}

    const std::vector<size_t>    &GetWaveformIndices() const      {return _waveformIndices;}
    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    private:

    float  _PEs;
    float  _PEsPulseHeight;  //used for PEs which were calculated using the pulse height and the pulse height calibration factor
    double _pulseTime;
    float  _pulseHeight;
    float  _pulseBeta;
    float  _pulseFitChi2;
    double _LEtime;
    CrvRecoPulseFlags  _flags;

    float   _PEsNoFit;        //based on the sum of the pedestal-subtracted ADC values of the pulse.
    double  _pulseTimeNoFit;  //time of largest ADC value. 
    double  _pulseStart;      //based on the time when the pulse starts to be above a threshold (FWHM). 
    double  _pulseEnd;

    std::vector<size_t>            _waveformIndices;  //indices in the vector of the CrvDigiCollection (which is the same as the index in the CrvDigiMCCollection)
    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
  };
  typedef std::vector<mu2e::CrvRecoPulse> CrvRecoPulseCollection;
}

#endif /* RecoDataProducts_CrvRecoPulse_hh */
