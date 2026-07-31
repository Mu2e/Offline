#ifndef CaloConditions_CalSimParamsMaker_hh
#define CaloConditions_CalSimParamsMaker_hh

//
// construct a CalSimParams conditions entity
// from fcl or database
//

#include "Offline/CaloConditions/inc/CalSimParams.hh"
#include "Offline/CaloConfig/inc/CalSimParamsConfig.hh"
#include "Offline/DbTables/inc/CalSimCrystals.hh"

namespace mu2e {

  class CalSimParamsMaker {
    typedef std::shared_ptr<CalSimParams> ptr_t;

    public:
      CalSimParamsMaker(CalSimParamsConfig const& config):_config(config) {}
      ptr_t fromFcl();
      ptr_t fromDb(CalSimCrystals::cptr_t cch_p);

    private:
      // this object needs to be thread safe,
      // _config should only be initialized once
      const CalSimParamsConfig _config;

  };
}


#endif
