//
// Class to describe a combination of Straw Hits
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {

  Float_t ComboHit::posRes(edir dir) const {
    switch ( dir ) {
      case ComboHit::wire : {
	return _wres;
      }
      case ComboHit::trans : {
	return _tres;
      }
      default : {
	return -1.0;
      }
    }
  }

  ComboHit::ComboHit() : _wres(-1.0),_tres(-1.0), _nsh(0), _sh{0}, _time(0.0), _edep(0.0), _qual(0.0) {}

}
