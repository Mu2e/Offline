#ifndef DataProducts_StrawEnd_hh
#define DataProducts_StrawEnd_hh
#include <iostream>
//
// define the wire ends
namespace mu2e {
  struct StrawEnd {
    enum strawend{unknown=-1,minus=0,plus=1,nends}; // azimuth of straw end WRT center
    strawend _end;
    StrawEnd(strawend end) : _end(end) {}
    bool operator ==(strawend end) const { return _end == end; }
    bool operator ==(StrawEnd const& other) const { return other._end == _end; }
    bool operator !=(StrawEnd const& other) const { return other._end != _end; }
  
    friend std::ostream& operator << (std::ostream& os, StrawEnd const& strawend) {
      switch (strawend._end) {
      case unknown: default:
	os << "Unknown";
	break;
      case plus:
	os << "Plus";
	break;
      case minus:
	os << "Minus";
	break;
      }
      return os;
    }
  };
}
#endif

