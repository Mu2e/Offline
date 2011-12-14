//
// $Id: VirtualDetector.cc,v 1.7 2011/12/14 00:30:16 gandr Exp $
// $Author: gandr $
//

#include "VirtualDetectorGeom/inc/VirtualDetector.hh"

namespace mu2e {

  VirtualDetector::VirtualDetector():
    _halfLength(0.01)
  {}

  void VirtualDetector::addVirtualDetector( int id, const std::string& name,
                                            const CLHEP::Hep3Vector& posParent,
                                            const CLHEP::HepRotation *rotParent,
                                            const CLHEP::Hep3Vector& posLocal) {
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

