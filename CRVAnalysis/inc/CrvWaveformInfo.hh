#ifndef CrvWaveformInfo_hh
#define CrvWaveformInfo_hh

#include <vector>
#include "Rtypes.h"

namespace mu2e
{

  struct CrvWaveformInfo   //information about CRV waveforms
  {
    Float_t             _adc;
    Float_t             _time;
    Int_t               _SiPMId;
    CrvWaveformInfo(float adc, float time, int SiPMId) :
                _adc(adc),
                _time(time),
                _SiPMId(SiPMId)
                {}
    CrvWaveformInfo() :
                _adc(0),
                _time(0),
                _SiPMId(0)
                {}
  };

  typedef std::vector<CrvWaveformInfo> CrvWaveformInfoCollection;  //this is the reco vector which will be stored in the main TTree

}
#endif
