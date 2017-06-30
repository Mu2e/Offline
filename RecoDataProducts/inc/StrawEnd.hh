#ifndef DataProducts_StrawEnd_hh
#define DataProducts_StrawEnd_hh
#include <iostream>
#include "TrackerConditions/inc/Types.hh"
//
// define the wire ends
namespace mu2e {
  struct StrawEnd {
    size_t _end;
    operator size_t() const { return _end; }

    StrawEnd(TrkTypes::End end=TrkTypes::cal) : _end(static_cast<size_t>(end)) {}
    bool operator == (TrkTypes::End end) const { return _end == static_cast<size_t>(end); }
    bool operator != (TrkTypes::End end) const { return _end != static_cast<size_t>(end); }
    bool operator == (StrawEnd const& other) const { return other._end == _end; }
    bool operator != (StrawEnd const& other) const { return other._end != _end; }

    friend std::ostream& operator << (std::ostream& os, StrawEnd const& end) {
      switch ( end._end ) {
	default:
	  os << "Unknown";
	  break;
	case TrkTypes::cal:
	  os << "Cal";
	  break;
	case TrkTypes::hv:
	  os << "HV";
	  break;
      }
      return os;
    }
  };
}
#endif

