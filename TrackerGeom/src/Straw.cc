//
// Hold information about one Straw.
//
//
// $Id: Straw.cc,v 1.1 2010/02/07 00:29:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/07 00:29:41 $
//
// Original author Rob Kutschke
//

#include <sstream>
#include "TrackerGeom/inc/Straw.hh"

#ifndef __CINT__ 

using std::vector;
using CLHEP::Hep3Vector;

namespace mu2e {

Straw::Straw():
  _id(StrawId()),
  _index(StrawIndex(0)),
  _c(Hep3Vector(0.,0.,0.)),
  _detail(0),
  _detailIndex(0),
  _w(Hep3Vector(0.,0.,1.))
{
}


Straw::Straw( const StrawId& id,
	      StrawIndex  index,
	      Hep3Vector const& c,
	      const StrawDetail* detail,
	      int detailIndex,
	      double wtx,
	      double wty
	      ):
  _id(id),
  _index(index),
  _c(c),
  _detail(detail),
  _detailIndex(detailIndex)
{
  _w = Hep3Vector(wtx,wty,1.).unit();
}

Straw::Straw( const StrawId& id,
	      StrawIndex index,
	      Hep3Vector const& c,
	      const StrawDetail* detail,
	      int detailIndex,
	      Hep3Vector const& w
	      ):
  _id(id),
  _index(index),
  _c(c),
  _detail(detail),
  _detailIndex(detailIndex),
  _w(w.unit())
{
}

Straw::~Straw (){
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
     << _id.getDevice() << "_"
     << _id.getSector() << "_"
     << _id.getLayer()  << "_"
     << _id.getStraw();
  return os.str();
}

} // namespace mu2e

#endif
