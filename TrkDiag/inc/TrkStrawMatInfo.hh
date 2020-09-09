//
// Struct to describe straw materials on the Kalman fit
//
#ifndef TrkStrawMatInfo_HH
#define TrkStrawMatInfo_HH
#include "Rtypes.h"
namespace mu2e 
{
  struct TrkStrawMatInfo {
    Bool_t _active;	    // was this material was used in the Kalman fit or not
    Int_t _plane, _panel, _layer, _straw; // StrawId fields of the straw
    Float_t _doca;    // DOCA between the track fit and the wire
    Float_t _tlen;   // length along the helix of the POCA
    Float_t _dp;      // momentum (energy) loss induced by this straw's material, including both entry and exit wall and the gas
    Float_t _radlen;  // radiation length of this straw's material seen by the track (including angular effects)
    TrkStrawMatInfo() : _active(false),
    _plane(-1), _panel(-1), _layer(-1), _straw(-1),
    _doca(-1000.0), _tlen(-1000.0),
    _dp(-1000.0), _radlen(-1000.0) {}
  };
}
#endif
