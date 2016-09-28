//
// Class to describe derived information from a StrawHit, in particular position.
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitPosition.hh"

namespace mu2e {

  Float_t StrawHitPosition::posRes(edir dir) const {
    switch ( dir ) {
      case StrawHitPosition::wire : {
	return _wres;
      }
      case StrawHitPosition::trans : {
	return _tres;
      }
      default : {
	return -1.0;
      }
    }
  }

  StrawHitPosition::StrawHitPosition() : _wdist(0.0), _wres(-1.0),_tres(-1.0), _stindex(-1) {}

}
