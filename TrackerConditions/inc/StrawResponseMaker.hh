#ifndef TrackerConditions_StrawResponseMaker_hh
#define TrackerConditions_StrawResponseMaker_hh

//
// construct a StrawResponse conditions entity
// from fcl or database
//

#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConfig/inc/StrawResponseConfig.hh"


namespace mu2e {

  class StrawResponseMaker {
  public:
    StrawResponseMaker(StrawResponseConfig const& config):_config(config) {}
    StrawResponse::ptr_t fromFcl(StrawDrift::cptr_t strawDrift,
				 StrawElectronics::cptr_t strawElectronics,
				 StrawPhysics::cptr_t strawPhysics);
    StrawResponse::ptr_t fromDb( /* db tables will go here*/ );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const StrawResponseConfig _config;

  };
}

#endif

