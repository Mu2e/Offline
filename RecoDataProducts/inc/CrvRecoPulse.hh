#ifndef RecoDataProducts_CrvRecoPulse_hh
#define RecoDataProducts_CrvRecoPulse_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulseFlags.hh"

#include <vector>

namespace mu2e
{
  class CrvRecoPulse
  {
    public:

    CrvRecoPulse() {}

    CrvRecoPulse(float PEs, float PEsPulseHeight, double pulseTime, float pulseHeight, float pulseBeta, float pulseFitChi2, double LEtime,
                 const CrvRecoPulseFlags &flags,
                 float PEsNoFit, double pulseTimeNoFit, double pulseStart, double pulseEnd,
                 const std::vector<size_t> &waveformIndices, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber,
		 float pedestal, bool pedestalFromDB) :
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
                                                                            _SiPMNumber(SiPMNumber),
                                                                            _pedestal(pedestal),
                                                                            _pedestalFromDB(pedestalFromDB)
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
    std::vector<size_t>    &GetWaveformIndices() {return _waveformIndices;} // used in reco compression
    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}
    float                         GetPedestal() const             {return _pedestal;}
    bool                          IsPedestalFromDB() const        {return _pedestalFromDB;}

    private:

    float  _PEs{0};
    float  _PEsPulseHeight{0};  //used for PEs which were calculated using the pulse height and the pulse height calibration factor
    double _pulseTime{0};
    float  _pulseHeight{0};
    float  _pulseBeta{0};
    float  _pulseFitChi2{0};
    double _LEtime{0};
    CrvRecoPulseFlags  _flags;

    float   _PEsNoFit{0};        //based on the sum of the pedestal-subtracted ADC values of the pulse.
    double  _pulseTimeNoFit{0};  //time of largest ADC value.
    double  _pulseStart{0};      //based on the time when the pulse starts to be above a threshold (FWHM).
    double  _pulseEnd{0};

    std::vector<size_t>            _waveformIndices;  //indices in the vector of the CrvDigiCollection (which is the same as the index in the CrvDigiMCCollection)
    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber{0};
    float                          _pedestal;
    bool                           _pedestalFromDB;
  };
  typedef std::vector<mu2e::CrvRecoPulse> CrvRecoPulseCollection;
}

#endif /* RecoDataProducts_CrvRecoPulse_hh */
