//
// Hold information about one Layer in a tracker.
//
//
// $Id: Layer.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "LTrackerGeom/inc/Layer.hh"

#ifndef __CINT__ 

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  Layer::Layer():
    _id(LayerId()),
    _nStraws(0),
    _orig(Hep3Vector(0.,0.,0.)),
    _delta(Hep3Vector(0.,0.,0.))
  {
  }
  
  Layer::Layer(const LayerId& id,
	       int      nStraws,
	       const Hep3Vector& origin,
	       const Hep3Vector& delta
	       ):
    _id(id),
    _nStraws(nStraws),
    _orig(origin),
    _delta(delta)
  {
  }
  
  Layer::Layer(const LayerId& id ):
    _id(id)
  {
  }
  
  Layer::~Layer(){}
  
} // namespace mu2e 

#endif
  
