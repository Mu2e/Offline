#ifndef VIRTUALDETECTORMAKER_HH
#define VIRTUALDETECTORMAKER_HH
//
// Construct and return an VirtualDetector.
//
// Original author Peter Shanahan
//

#include <vector>
#include <string>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class VirtualDetector;
class SimpleConfig;

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

#endif 
