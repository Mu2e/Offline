//
// Hold information about one Straw.
//
// Original author Rob Kutschke
//

#include <sstream>
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"

#ifndef __CINT__

using std::vector;
using CLHEP::Hep3Vector;

namespace mu2e {

  Straw::Straw():
    _id(StrawId()),
    _c(CLHEP::Hep3Vector(0.,0.,0.)),
    _detail(0),
    _detailIndex(0),
    _w(CLHEP::Hep3Vector(0.,0.,1.))
  {
  }

  Straw::Straw( const StrawId& id,
                CLHEP::Hep3Vector const& c,
                const StrawDetail* detail,
                int detailIndex,
                double wtx,
                double wty
                ):
    _id(id),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex)
  {
    _w = CLHEP::Hep3Vector(wtx,wty,1.).unit();
  }

  Straw::Straw( const StrawId& id,
                CLHEP::Hep3Vector const& c,
                const StrawDetail* detail,
                int detailIndex,
                CLHEP::Hep3Vector const& w
                ):
    _id(id),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex),
    _w(w.unit())
  {
  }

  void Straw::fillPointers ( const Tracker& tracker ) const{
    _detail = &tracker.getStrawDetails().at(_detailIndex);
  }

  std::string Straw::name( std::string const& base ) const{
    std::ostringstream os;

    os << base
       << _id.getPlane() << "_"
       << _id.getPanel() << "_"
       << _id.getLayer()  << "_"
       << _id.getStraw();
    return os.str();
  }

} // namespace mu2e
#endif
