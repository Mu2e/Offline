//
// Struct to describe straw materials on the Kalman fit
// $Id: TrkStrawMatInfo.hh,v 1.2 2014/09/22 12:13:17 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/22 12:13:17 $
//
#ifndef TrkStrawMatInfo_HH
#define TrkStrawMatInfo_HH
#include "Rtypes.h"
namespace mu2e 
{
  struct TrkStrawMatInfo {
    Bool_t _active;	    // was this material was used in the Kalman fit or not
    Bool_t _thit, _thita;   // Is this material associated with a hit (or active hit)?
    Int_t _plane, _panel, _layer, _straw; // StrawId fields of the straw
    Float_t _doca;    // DOCA between the track fit and the wire
    Float_t _tlen;   // length along the helix of the POCA
    Float_t _dp;      // momentum (energy) loss induced by this straw's material, including both entry and exit wall and the gas
    Float_t _radlen;  // radiation length of this straw's material seen by the track (including angular effects)
    Float_t _sigMS;   // RMS of the (1-dimensional) scattering angle induced by this straw's material
    TrkStrawMatInfo() : _active(false), _thit(false), _thita(false),
    _plane(-1), _panel(-1), _layer(-1), _straw(-1),
    _doca(-1000.0), _tlen(-1000.0),
    _dp(-1000.0), _radlen(-1000.0), _sigMS(-1000.0) {}
  };
}
#endif
