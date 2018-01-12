#ifndef RecoDataProducts_StrawHitPosition_hh
#define RecoDataProducts_StrawHitPosition_hh
//
// Class to describe derived information from a StrawHit, in particular pos().
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitFlag.hh"
// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
// root includes
#include "Rtypes.h"
// C++ includes
#include <vector>
namespace mu2e {

  struct StrawHitPosition {
    enum edir {wire=0, trans };
    StrawHitPosition();
    // accessors
    CLHEP::Hep3Vector const& pos() const { return _pos; }
    CLHEP::Hep3Vector const& wdir() const { return _wdir; }
    // center of the wire
    CLHEP::Hep3Vector centerPos() const;
    Float_t wireDist() const { return _wdist; }
    Float_t posRes(edir dir) const;
    Float_t phi() const { return _phi;}
    Int_t stereoHitIndex() const { return _stindex; } // negative if there's no stereo hit
    StrawHitFlag const& flag() const { return _flag; }
    CLHEP::Hep3Vector _pos; // pos() of the hit.  This is on the wire
    CLHEP::Hep3Vector _wdir; // wire direction at this position (unit vector)
    Float_t _phi; //cache phi angle to speed up processing
    Float_t _wdist; // distance along the wire from the center
    Float_t _wres; // resolution along the wire
    Float_t _tres; // resolution perpendicular to the wire (= straw radius)
    Int_t _stindex; // index into stereo hit collection (-1 if not based on stereo)
    StrawHitFlag _flag; // bit flags for this hit position
  };
  typedef std::vector<mu2e::StrawHitPosition> StrawHitPositionCollection;
}
#endif


