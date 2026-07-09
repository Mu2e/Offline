#ifndef CaloConditions_CalSimCrystalMaker_hh
#define CaloConditions_CalSimCrystalMaker_hh

//
// construct a CalSimCrystal conditions entity
// from fcl or database
//

#include "Offline/CaloConditions/inc/CalSimCrystal.hh"
#include "Offline/CaloConfig/inc/CalSimCrystalConfig.hh"
#include "Offline/DbTables/inc/CalSimCrystals.hh"

namespace mu2e {

  class CalSimCrystalMaker {
    typedef std::shared_ptr<CalSimCrystal> ptr_t;

    public:
      CalSimCrystalMaker(CalSimCrystalConfig const& config):_config(config) {}
      ptr_t fromFcl();
      ptr_t fromDb(CalSimCrystals::cptr_t cch_p);

    private:
      // this object needs to be thread safe,
      // _config should only be initialized once
      const CalSimCrystalConfig _config;

  };
}


#endif
