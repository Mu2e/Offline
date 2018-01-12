#ifndef RecoDataProducts_ComboHit_hh
#define RecoDataProducts_ComboHit_hh
//
// Class to describe a combination of several StrawHits, either as adjacent
// hits in a panel or stereo hits (or a combination of these)
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/XYZVec.hh"
#include <stdint.h>
// root includes
#include "Rtypes.h"
// C++ includes
#include <array>
#include <vector>
namespace mu2e {

  struct ComboHit {
    enum edir {wire=0, trans };
    ComboHit();
    // accessors
    XYZVec const& pos() const { return _pos; }
    XYZVec const& wdir() const { return _wdir; }
    Float_t posRes(edir dir) const;
    Float_t energyDep() const { return _edep; }
    Float_t time() const { return _time; }
    Float_t qual() const { return _qual; }
    StrawHitFlag const& flag() const { return _flag; }
    uint16_t nCombo() const { return _nsh; }
    StrawHitIndex shIndex(uint16_t ish) const {
      if(ish < _nsh)
	return _sh[ish];
      else
	return 0; // should throw here FIXME!
    }
    //
    XYZVec _pos; // average position
    XYZVec _wdir; // direction at this position (typically the wire direction)
    Float_t _wdist; // distance from wire center along this direction
    Float_t _wres; // resolution along this direction
    Float_t _tres; // resolution perpendicular to this direction
    uint16_t _nsh; // number of associated straw hits
    constexpr static size_t MaxNStraws = 8;
    typedef std::array<StrawHitIndex,MaxNStraws> SHIArray;
    SHIArray _sh; // Indices back to straw hits
    Float_t _time; // Average time for these
    Float_t _edep; // average energy deposit for these
    Float_t _qual; // quality of combination
    StrawHitFlag _flag; // bit flags for this combo
  };
  typedef std::vector<mu2e::ComboHit> ComboHitCollection;
}
#endif


