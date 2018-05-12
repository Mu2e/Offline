#ifndef DataProducts_StrawEnd_hh
#define DataProducts_StrawEnd_hh
#include <ostream>
#include <Rtypes.h>
//
// define the ends of a straw
namespace mu2e {
  struct StrawEnd {
    enum End{unknown=-1,cal=0,hv, nends};
    int8_t _end;
    operator int8_t() const { return _end; }
    End end() const { return static_cast<End>(_end); }
    StrawEnd(End end=cal) : _end(static_cast<int8_t>(end)) {}
    bool operator == (End end) const { return _end == end; }
    bool operator != (End end) const { return _end != end; }
    bool operator == (StrawEnd const& other) const { return other._end == _end; }
    bool operator != (StrawEnd const& other) const { return other._end != _end; }

    friend std::ostream& operator << (std::ostream& os, StrawEnd const& end) {
      switch ( end._end ) {
	case cal:
	  os << "Cal";
	  break;
	case hv:
	  os << "HV";
	  break;
      	default:
	  os << "Unknown";
	  break;
      }
      return os;
    }
  };
}
#endif

