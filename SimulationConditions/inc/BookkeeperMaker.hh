#ifndef SimulationConditions_BookkeeperMaker_hh
#define SimulationConditions_BookkeeperMaker_hh

//
// Makes the Bookkeeper ProditionsEntitiy
//

#include "SimulationConditions/inc/Bookkeeper.hh"
#include "SimulationConfig/inc/BookkeeperConfig.hh"
#include "DbTables/inc/SimEfficiencies.hh"

namespace mu2e {

  class BookkeeperMaker {
  public:
    BookkeeperMaker(BookkeeperConfig const& config):_config(config) {}

    typename Bookkeeper::ptr_t fromFcl() {
      auto ptr = std::make_shared<Bookkeeper>();
      for (const auto& i_effConf : _config.simStageEfficiencies()) {
	ptr->addEff(i_effConf.tag(), i_effConf.eff());
      }
      return ptr;
    }

    typename Bookkeeper::ptr_t fromDb(typename SimEfficiencies::cptr_t effDb) {
      // fill the Bookkeeper with initial values
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
    const BookkeeperConfig _config;
  };
}

#endif
