//
// Struct describing MC truth of single straw hit assigned to a track, for use in TTree diagnostics
// $Id: TrkStrawHitInfo.hh,v 1.2 2014/09/22 12:13:17 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/22 12:13:17 $
//
#ifndef TrkStrawHitInfoMC_HH
#define TrkStrawHitInfoMC_HH
#include "Rtypes.h"
namespace mu2e 
{
   struct TrkStrawHitInfoMC {
    Int_t _pdg, _gen, _proc; // PDG particle code, generator code and process code of the particle which caused the electronics to cross threshold in simulation
    Int_t _rel; // relationship (same, mother, daughter, sibling, unrelated) of this particle to the particle generating most of the hits on this track
    Float_t _t0;  // true time this particle passed closest to this wire
    Float_t _ht;  // true time this particle's signal reach the electronics (includes true drift time)
    Float_t _dist;  // true transverse distance between the particle trajectory and the wire
    Float_t _len;   // true distance from the straw center to the Point of Minimum Transverse Distance (PMTD)
    Float_t _edep;  // true energy deposited by this particle in the straw gas
    Float_t _mom;   // true particle momentum at the PMTD
    Float_t _r, _phi; // polar coordinates of the PMTD
    Int_t _ambig;   // true left-right ambiguity = true angular momentum sign of the particle WRT the wire
    Bool_t _xtalk;  // whether or not this hit was caused by cross-talk of a real signal.
    TrkStrawHitInfoMC() : _pdg(-1), _gen(-1), _proc(-1), 
    _t0(-1000.0), _ht(-1000.0), _dist(-1000.0), _len(-1000.0),
    _edep(-1000.0),_r(-1000.0),_phi(-1000.0),
    _ambig(-100), _xtalk(false) {}
  };
}
#endif
