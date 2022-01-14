//
// Makes the SimBookkeeper ProditionsEntitiy
//

#include "Offline/SimulationConditions/inc/SimBookkeeperMaker.hh"

namespace mu2e {

  SimBookkeeper::ptr_t SimBookkeeperMaker::fromFcl() {
    auto ptr = std::make_shared<SimBookkeeper>();
    for (const auto& i_effConf : _config.simStageEfficiencies()) {
      ptr->addEff(i_effConf.tag(), i_effConf.eff());
    }
    return ptr;
  }

  SimBookkeeper::ptr_t SimBookkeeperMaker::fromDb(SimEfficiencies2::cptr_t effDb) {
    // fill the SimBookkeeper with initial values
    auto ptr = fromFcl();
    // now overwrite with values from database
    for (const auto& i_row : effDb->rows()) {
      ptr->addEff(i_row.tag(), i_row.eff());
    }
    return ptr;
  }
}
