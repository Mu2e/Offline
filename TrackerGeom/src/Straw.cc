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
    _id2(StrawId2()),
    _index(StrawIndex(0)),
    _index2(StrawIndex2(0)),
    _c(CLHEP::Hep3Vector(0.,0.,0.)),
    _detail(0),
    _detailIndex(0),
    _w(CLHEP::Hep3Vector(0.,0.,1.))
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
    _id2(StrawId2()),
    _index(index),
    _index2(StrawIndex2(0)),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex)
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
    _id2(StrawId2()),
    _index(index),
    _index2(StrawIndex2(0)),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex),
    _w(w.unit())
  {
  }

  Straw::Straw( const StrawId& id, const StrawId2& id2,
                StrawIndex index, StrawIndex2 index2,
                CLHEP::Hep3Vector const& c,
                const StrawDetail* detail,
                int detailIndex,
                CLHEP::Hep3Vector const& w
                ):
    _id(id),
    _id2(id2),
    _index(index),
    _index2(index2),
    _c(c),
    _detail(detail),
    _detailIndex(detailIndex),
    _w(w.unit())
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

  bool Straw::isSamePreamp( StrawIndex idx ) const{
   for ( vector<StrawIndex>::const_iterator i=_preampByIndex.begin(),
            e=_preampByIndex.end();
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
