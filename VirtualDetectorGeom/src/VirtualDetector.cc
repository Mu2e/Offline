//
// $Id: VirtualDetector.cc,v 1.10 2011/12/14 22:20:47 gandr Exp $
// $Author: gandr $
//

#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

namespace mu2e {

  VirtualDetector::VirtualDetector():
    _halfLength(0.01)
  {}

  const std::string VirtualDetector::_baseName("VirtualDetector");

  std::string VirtualDetector::volumeName(int i) { return _baseName + "_" + VirtualDetectorId(i).name(); }

  void VirtualDetector::addVirtualDetector( int id,
                                            const CLHEP::Hep3Vector& posParent,
                                            const CLHEP::HepRotation *rotParent,
                                            const CLHEP::Hep3Vector& posLocal) {
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

