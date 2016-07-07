#ifndef TrkDiag_TimePeakHitInfo_hh
#define TrkDiag_TimePeakHitInfo_hh
#include "Rtypes.h"

using namespace std; 

namespace mu2e 
{
// struct for time peak hit diagnostics
  struct TimePeakHitInfo {
    Float_t _dt; // time relative to the peak
    Float_t _dphi; // resolved phi relative to the peak
    Float_t _rho; // transverse radius of the peak
    Float_t _mva; // peak MVA output
    Int_t _mcpdg, _mcgen, _mcproc; // MC truth info
  };
}
