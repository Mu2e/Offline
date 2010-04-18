//
// Hold information about one Layer in a tracker.
//
//
// $Id: Layer.cc,v 1.2 2010/04/18 00:31:56 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:31:56 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Layer.hh"
#include "TrackerGeom/inc/Tracker.hh"

#ifndef __CINT__ 

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  Layer::Layer():
    _id(LayerId()),
    _nStraws(0),
    _orig(Hep3Vector(0.,0.,0.)),
    _delta(Hep3Vector(0.,0.,0.)){
  }

  Layer::Layer(const LayerId& id,
               int      nStraws,
               const Hep3Vector& origin,
               const Hep3Vector& delta
               ):
    _id(id),
    _nStraws(nStraws),
    _orig(origin),
    _delta(delta){
  }

  Layer::Layer(const LayerId& id ):
    _id(id){
  }

  string Layer::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getDevice() << "_"
       << _id.getSector() << "_"
       << _id.getLayer();
    return os.str();
  }

  void Layer::fillPointers ( const Tracker& tracker ) const{
    _straws.clear();
    for ( size_t i=0; i<_indices.size(); ++i ){
      StrawIndex idx = _indices[i];
      const Straw* straw =  &tracker.getStraw(idx);
      _straws.push_back(straw);
      straw->fillPointers(tracker);
    }
  }

} // namespace mu2e 
#endif
  
