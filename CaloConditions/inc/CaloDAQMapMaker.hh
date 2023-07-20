#ifndef CaloConditions_CaloDAQMapMaker_hh
#define CaloConditions_CaloDAQMapMaker_hh

//
// construct a CaloDAQMap conditions entity
// from fcl or database
//

#include "Offline/CaloConditions/inc/CaloDAQMap.hh"
#include "Offline/CaloConfig/inc/CaloDAQMapConfig.hh"


namespace mu2e {

  class CaloDAQMapMaker {
  typedef std::shared_ptr<CaloDAQMap> ptr_t;

  public:
    CaloDAQMapMaker(CaloDAQMapConfig const& config):_config(config) {}
    ptr_t fromFcl();
    ptr_t fromDb();

  private:

    // this object needs to be thread safe,
    // _config should only be initialized once
    const CaloDAQMapConfig _config;

  };
}


#endif
