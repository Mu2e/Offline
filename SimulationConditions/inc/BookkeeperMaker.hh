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
      ptr->addEff("muBeamPerPOT", _config.muBeamPerPOT());
      ptr->addEff("flashPerMuBeam", _config.flashPerMuBeam());
      return ptr;
    }

    typename Bookkeeper::ptr_t fromDb(typename SimEfficiencies::cptr_t effDb) {
      // fill the Bookkeeper with initial values
      auto ptr = fromFcl();
      // now overwrite with values from database
      double eff;
      effDb->findEff("muBeamPerPOT", eff);
      ptr->setEffVal("muBeamPerPOT", eff);
      effDb->findEff("flashPerMuBeam", eff);
      ptr->setEffVal("flashPerMuBeam", eff);
      return ptr;
    }

  private:
    // this object needs to be thread safe, 
    // _config should only be initialized once
    const BookkeeperConfig _config;
  };
}

#endif
