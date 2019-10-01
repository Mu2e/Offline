#ifndef TrackerConditions_Mu2eDetectorMaker_hh
#define TrackerConditions_Mu2eDetectorMaker_hh
//
// Make Mu2eDetector from fcl or (eventually) database
//

#include "TrackerConditions/inc/Mu2eMaterial.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
#include "TrackerConditions/inc/Mu2eDetectorConfig.hh"

namespace mu2e {

  class Mu2eDetectorMaker {

  public:
    Mu2eDetectorMaker(Mu2eDetectorConfig const& config):_config(config) {}
    Mu2eDetector::ptr_t fromFcl(Mu2eMaterial::cptr_t material, 
				Tracker::cptr_t det);

  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const Mu2eDetectorConfig _config;
  };


} // namespace mu2e

#endif /* TrackerConditions_Mu2eDetectorMaker_hh */
