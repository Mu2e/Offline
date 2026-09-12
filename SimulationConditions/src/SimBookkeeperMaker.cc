//
// Makes the SimBookkeeper ProditionsEntitiy
//

#include "Offline/SimulationConditions/inc/SimBookkeeperMaker.hh"
#include "cetlib_except/exception.h"

#include <set>

namespace mu2e {

  SimBookkeeper::ptr_t SimBookkeeperMaker::fromFcl() {
    auto ptr = std::make_shared<SimBookkeeper>();
    std::set<std::string> fclTags;
    for (const auto& i_effConf : _config.simStageEfficiencies()) {
      if (!fclTags.insert(i_effConf.tag()).second) {
        throw cet::exception("SIMBOOKKEEPER_DUPLICATE_FCL_TAG")
          << "Efficiency tag \"" << i_effConf.tag()
          << "\" appears more than once in simStageEfficiencies\n";
      }
      ptr->addEff(i_effConf.tag(), i_effConf.eff());
    }
    return ptr;
  }

  SimBookkeeper::ptr_t SimBookkeeperMaker::fromDb(SimEfficiencies2::cptr_t effDb) {
    // fill the SimBookkeeper with initial values
    auto ptr = fromFcl();
    // now overwrite with values from database
    std::set<std::string> dbTags;
    for (const auto& i_row : effDb->rows()) {
      if (!dbTags.insert(i_row.tag()).second) {
        throw cet::exception("SIMBOOKKEEPER_DUPLICATE_DB_TAG")
          << "Efficiency tag \"" << i_row.tag()
          << "\" appears more than once in table " << effDb->name() << "\n";
      }
      ptr->addEff(i_row.tag(), i_row.eff());
    }
    return ptr;
  }
}
