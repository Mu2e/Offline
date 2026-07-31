#ifndef CaloConditions_CalDAQMapMaker_hh
#define CaloConditions_CalDAQMapMaker_hh

//
// construct a CalDAQMap conditions entity
// from fcl or database
//

#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/CaloConfig/inc/CalDAQMapConfig.hh"
#include "Offline/DbTables/inc/CalChannels.hh"

namespace mu2e {

  class CalDAQMapMaker {
  typedef std::shared_ptr<CalDAQMap> ptr_t;

  public:
    CalDAQMapMaker(CalDAQMapConfig const& config):_config(config) {}
    ptr_t fromFcl();
    ptr_t fromDb(CalChannels::cptr_t cch_p);

  private:

    // this object needs to be thread safe,
    // _config should only be initialized once
    const CalDAQMapConfig _config;

  };
}


#endif
