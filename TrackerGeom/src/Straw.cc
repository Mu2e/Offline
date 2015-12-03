//
// Hold information about one Straw.
//
// $Id: Straw.cc,v 1.5 2011/06/01 16:02:58 mf Exp $
// $Author: mf $
// $Date: 2011/06/01 16:02:58 $
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
    _index(StrawIndex(0)),
    _c(CLHEP::Hep3Vector(0.,0.,0.)),
    _detail(0),
    _detailIndex(0),
    _w(CLHEP::Hep3Vector(0.,0.,1.)),
    _nextOuterL(StrawIndex::NO_STRAW),
    _nextInnerL(StrawIndex::NO_STRAW),
    _nextOuterP(StrawIndex::NO_STRAW),
    _nextInnerP(StrawIndex::NO_STRAW)
  {
  }

  Straw::Straw( const StrawId& id,
                StrawIndex  index,
                CLHEP::Hep3Vector const& c,
                const StrawDetail* detail,
                int detailIndex,
                double wtx,
                double wty
                ):
    _id(id),
    _index(index),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex),
    _nextOuterL(StrawIndex::NO_STRAW),
    _nextInnerL(StrawIndex::NO_STRAW),
    _nextOuterP(StrawIndex::NO_STRAW),
    _nextInnerP(StrawIndex::NO_STRAW)
  {
    _w = CLHEP::Hep3Vector(wtx,wty,1.).unit();
  }

  Straw::Straw( const StrawId& id,
                StrawIndex index,
                CLHEP::Hep3Vector const& c,
                const StrawDetail* detail,
                int detailIndex,
                CLHEP::Hep3Vector const& w
                ):
    _id(id),
    _index(index),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex),
    _w(w.unit()),
    _nextOuterL(StrawIndex::NO_STRAW),
    _nextInnerL(StrawIndex::NO_STRAW),
    _nextOuterP(StrawIndex::NO_STRAW),
    _nextInnerP(StrawIndex::NO_STRAW)
  {
  }

  void Straw::fillPointers ( const Tracker& tracker ) const{
    _detail = &tracker.getStrawDetails().at(_detailIndex);
  }

  bool Straw::isNearestNeighbour( StrawIndex idx ) const{

    for ( vector<StrawIndex>::const_iterator i=_nearestByIndex.begin(),
            e=_nearestByIndex.end();
          i<e; ++i ){
      if ( *i == idx ) return true;
    }

    return false;
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
