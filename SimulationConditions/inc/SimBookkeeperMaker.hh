#ifndef SimulationConditions_SimBookkeeperMaker_hh
#define SimulationConditions_SimBookkeeperMaker_hh

//
// Makes the SimBookkeeper ProditionsEntitiy
//

#include "Offline/SimulationConditions/inc/SimBookkeeper.hh"
#include "Offline/SimulationConfig/inc/SimBookkeeperConfig.hh"
#include "Offline/DbTables/inc/SimEfficiencies2.hh"

namespace mu2e {

  class SimBookkeeperMaker {
  public:
    SimBookkeeperMaker(SimBookkeeperConfig const& config):_config(config) {}

    SimBookkeeper::ptr_t fromFcl();
    SimBookkeeper::ptr_t fromDb(SimEfficiencies2::cptr_t effDb);

  private:
    // this object needs to be thread safe,
    // _config should only be initialized once
    const SimBookkeeperConfig _config;
  };
}

#endif
