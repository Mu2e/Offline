#ifndef CaloConditions_CaloDAQMapMaker_hh
#define CaloConditions_CaloDAQMapMaker_hh

//
// construct a CaloDAQMap conditions entity
// from fcl or database
//

#include "CaloConditions/inc/CaloDAQMap.hh"
#include "CaloConfig/inc/CaloDAQMapConfig.hh"
#include "DbTables/inc/CalRoIDMapDIRACToOffline.hh"
#include "DbTables/inc/CalRoIDMapOfflineToDIRAC.hh"


namespace mu2e {

  class CaloDAQMapMaker {
  typedef std::shared_ptr<CaloDAQMap> ptr_t;

  public:
    CaloDAQMapMaker(CaloDAQMapConfig const& config):_config(config) {}
    ptr_t fromFcl();
    ptr_t fromDb(CalRoIDMapDIRACToOffline::cptr_t tdtc,
		 CalRoIDMapOfflineToDIRAC::cptr_t tctd );
  
  private:

    // this object needs to be thread safe, 
    // _config should only be initialized once
    const CaloDAQMapConfig _config;

  };
}


#endif

