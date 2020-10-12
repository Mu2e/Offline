#ifndef SimulationConditions_SimBookkeeperMaker_hh
#define SimulationConditions_SimBookkeeperMaker_hh

//
// Makes the SimBookkeeper ProditionsEntitiy
//

#include "SimulationConditions/inc/SimBookkeeper.hh"
#include "SimulationConfig/inc/SimBookkeeperConfig.hh"
#include "DbTables/inc/SimEfficiencies.hh"

namespace mu2e {

  class SimBookkeeperMaker {
  public:
    SimBookkeeperMaker(SimBookkeeperConfig const& config):_config(config) {}

    SimBookkeeper::ptr_t fromFcl() {
      auto ptr = std::make_shared<SimBookkeeper>();
      for (const auto& i_effConf : _config.simStageEfficiencies()) {
        ptr->addEff(i_effConf.tag(), i_effConf.eff());
      }
      return ptr;
    }

    SimBookkeeper::ptr_t fromDb(SimEfficiencies::cptr_t effDb) {
      // fill the SimBookkeeper with initial values
      auto ptr = fromFcl();
      // now overwrite with values from database
      for (const auto& i_row : effDb->rows()) {
        ptr->addEff(i_row.tag(), i_row.eff());
      }
      return ptr;
    }

  private:
    // this object needs to be thread safe,
    // _config should only be initialized once
    const SimBookkeeperConfig _config;
  };
}

#endif
