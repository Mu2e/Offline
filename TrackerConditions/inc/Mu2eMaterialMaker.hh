#ifndef TrackerConditions_Mu2eMaterialMaker_hh
#define TrackerConditions_Mu2eMaterialMaker_hh
//
// Make Mu2eMaterial from fcl
// Since Mu2eMaterial holds pointers to singletons
// inside of BTrk, there should only ever be one made
//

#include "TrackerConditions/inc/Mu2eMaterial.hh"
#include "TrackerConditions/inc/Mu2eMaterialConfig.hh"

namespace mu2e {

  class Mu2eMaterialMaker {

  public:
    Mu2eMaterialMaker(Mu2eMaterialConfig const& config):_config(config) {}
    Mu2eMaterial::ptr_t fromFcl();

  private:

    // _config should only be initialized once
    const Mu2eMaterialConfig _config;
  };


} // namespace mu2e

#endif /* TrackerConditions_Mu2eMaterialMaker_hh */
