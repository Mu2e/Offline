//
// Struct describing a calorimeter cluster assigned to a track, for use in TTree diagnostics
// original author: Dave Brown (LBNL), Jan 2019
//
#ifndef TrkCaloHitInfo_HH
#define TrkCaloHitInfo_HH
#include "Rtypes.h"
#include "DataProducts/inc/XYZVec.hh"
namespace mu2e 
{
  struct TrkCaloHitInfo {
    Int_t _active;   // is this hit used in the track?
    Int_t _did; // disk ID
    XYZVec _poca; // Position of Point Of Closest Approach (POCA)
    XYZVec _mom; // Track momentum vector at Point Of Closest Approach (POCA)
    Float_t _trklen;	// distance along the helix of the POCA for this hit
    Float_t _clen;    // length along the crystal from the front face
    Float_t _doca;	// DOCA of this hit
    Float_t _t0, _t0err;  // Estimated time the particle passed the POCA of this hit (and error).  Note this is mass hypothesis dependent
    Float_t _ct, _cterr;    // reconstructed time (including propagation time)
    Float_t _edep;    // reconstructed crystal energy deposition
    TrkCaloHitInfo() : _active(false),_did(-1),
    _trklen(-1000.0), _clen(-1000.0),
    _doca(-1000.0), _t0(-1000.0), _t0err(-1000.0),
    _ct(-1000.0), _cterr(-1000.0), _edep(-1000.0)  {}
    void reset() { *this = TrkCaloHitInfo(); }
    static std::string const& leafnames() { 
      static const std::string leaves = 
	std::string("active/I:disk/I:POCAx/F:POCAy/F:POCAz/F:momx/F:momy/F:momz/F:") +
	std::string("trklen/F:clen/F:doca/F:t0/F:t0err/F:ctime/F:cterr/F:edep/F");
      return leaves;
    }
  };
}
#endif
