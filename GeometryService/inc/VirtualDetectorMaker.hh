#ifndef VirtualDetectorGeom_VirtualDetectorMaker_hh
#define VirtualDetectorGeom_VirtualDetectorMaker_hh
//
// Construct and return an VirtualDetector.
//
// $Id: VirtualDetectorMaker.hh,v 1.1 2012/02/27 06:05:35 gandr Exp $
// $Author: gandr $
//
// Original author Peter Shanahan
//

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class SimpleConfig;
class VirtualDetector;

class VirtualDetectorMaker {

public:

  VirtualDetectorMaker( SimpleConfig const& config );

  ~VirtualDetectorMaker ();

  // This is depracted and will go away soon.
  // Still needed for root graphics version.
  const VirtualDetector& getVirtualDetector() const { return *_vd;}

  // This is the accessor that will remain.
  std::auto_ptr<VirtualDetector> getVirtualDetectorPtr() { return _vd; }

private:

  void BuildVirtualDetector(SimpleConfig const&);

  std::auto_ptr<VirtualDetector> _vd;

};

}  //namespace mu2e

#endif /* VirtualDetectorGeom_VirtualDetectorMaker_hh */
