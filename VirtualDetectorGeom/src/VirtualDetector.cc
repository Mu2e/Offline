//
// $Id: VirtualDetector.cc,v 1.5 2011/12/10 00:16:34 genser Exp $
// $Author: genser $
//

#include "VirtualDetectorGeom/inc/VirtualDetector.hh"

namespace mu2e {

  VirtualDetector::VirtualDetector():
    _halfLength(0.01)
  {}

  void VirtualDetector::addVirtualDetector( int id, std::string name,
                                            CLHEP::Hep3Vector posParent,
                                            CLHEP::HepRotation *rotParent,
                                            CLHEP::Hep3Vector posLocal) {
    _name[id] = name;
    _local[id] = posLocal;
    if( rotParent==0 ) {
      _global[id]   = posParent+posLocal;
      _rotation[id] = 0;
    } else {
      _rotation[id] = rotParent;
      _global[id]   = posParent+rotParent->inverse()*posLocal;
    }
  }
}

