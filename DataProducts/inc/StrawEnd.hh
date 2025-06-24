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
    End otherEnd() const { return (_end == cal) ? hv : cal; }
    StrawEnd(End end=cal) : _end(static_cast<int8_t>(end)) {}
    // convention for signing directions according to the end
    double endSign() const { return 2.0*(static_cast<double>(_end) - 0.5); }
    bool operator == (End end) const { return _end == end; }
    bool operator != (End end) const { return _end != end; }
    bool operator == (StrawEnd const& other) const { return other._end == _end; }
    bool operator != (StrawEnd const& other) const { return other._end != _end; }

  };
}
std::ostream& operator << (std::ostream& os, mu2e::StrawEnd const& end);
#endif

