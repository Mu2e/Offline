//
// Hold information about one Straw.
//
// Original author Rob Kutschke
//

#include <sstream>
#include "Offline/TrackerGeom/inc/Straw.hh"

using std::vector;

namespace mu2e {

  // unaligned construct
  Straw::Straw( const StrawId& id,
      const xyzVec& c,
      double halflen,
      double wtx,
      double wty
      ):
    Straw(id,c,xyzVec(wtx,wty,1.),halflen) {}

  Straw::Straw( const StrawId& id,
      const xyzVec& c,
      const xyzVec& w, double halfLength):
    _id(id),
    _wmid(c), _smid(c),_wdir(w.unit()), _sdir(w.unit()), _hlen(halfLength) {}


  // aligned constructor
  Straw::Straw (const StrawId& id,
      const xyzVec& calwireend, const xyzVec& hvwireend,
      const xyzVec& calstrawend, const xyzVec& hvstrawend) :
    _id(id),
    _wmid(0.5*(hvwireend + calwireend)),
    _smid(0.5*(hvstrawend + calstrawend)),
    _wdir((calwireend - hvwireend).unit()),  // convention is U points from HV to cal
    _sdir((calstrawend - hvstrawend).unit()),
    _hlen(0.5*(hvwireend - calwireend).mag()) {}

  std::string Straw::name( std::string const& base ) const{
    std::ostringstream os;
    os << base << _id;
    return os.str();
  }


} // namespace mu2e
