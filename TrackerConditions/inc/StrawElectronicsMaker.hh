#ifndef TrackerConditions_StrawElectronicsMaker_hh
#define TrackerConditions_StrawElectronicsMaker_hh

//
// construct a StrawElectronics conditions entity
// from fcl or database
//

#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawElectronicsConfig.hh"
#include "DbTables/inc/TrkDelayPanel.hh"
#include "DbTables/inc/TrkPreampRStraw.hh"
#include "DbTables/inc/TrkPreampStraw.hh"
#include "DbTables/inc/TrkThresholdRStraw.hh"


namespace mu2e {

  class StrawElectronicsMaker {
  public:
    StrawElectronicsMaker(StrawElectronicsConfig const& config):_config(config) {}
    StrawElectronics::ptr_t fromFcl();
    StrawElectronics::ptr_t fromDb(TrkDelayPanel::cptr_t tdp,
				   TrkPreampRStraw::cptr_t tprs,
				   TrkPreampStraw::cptr_t tps,
				   TrkThresholdRStraw::cptr_t ttrs );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const StrawElectronicsConfig _config;

  };
}


#endif

