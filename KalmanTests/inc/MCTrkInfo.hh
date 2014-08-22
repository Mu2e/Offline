#ifndef MCTrkInfo_HH
#define MCTrkInfo_HH
#include "KalmanTests/inc/threevec.hh"
#include "KalmanTests/inc/helixpar.hh"
#include "Rtypes.h"
namespace mu2e
{
  struct MCTrkInfo {
    Float_t _time;
    Float_t _mom;
    threevec _pos;
    helixpar _hpar;
    MCTrkInfo() :  _time(-1.0),_mom(-1.0) {}
    void reset() { _time = _mom = -1; _pos.reset(); _hpar.reset(); }
  };
}
#endif

