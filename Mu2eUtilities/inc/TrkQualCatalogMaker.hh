#ifndef Mu2eUtilities_TrkQualCatalogMaker_hh
#define Mu2eUtilities_TrkQualCatalogMaker_hh

//
// Makes the TrkQualCatalog ProditionsEntitiy
//

#include "Mu2eUtilities/inc/TrkQualCatalog.hh"
#include "Mu2eUtilities/inc/TrkQualCatalogConfig.hh"
#include "DbTables/inc/TrkQualDb.hh"

namespace mu2e {
  class TrkQualCatalogMaker {
  public:
    TrkQualCatalogMaker(TrkQualCatalogConfig const& config):_config(config) {}

    TrkQualCatalog::ptr_t fillEntries();
    void initializeMVAs(TrkQualCatalog::ptr_t ptr);
    TrkQualCatalog::ptr_t fromFcl();
    TrkQualCatalog::ptr_t fromDb(TrkQualDb::cptr_t tqDb);
  private:
    // this object needs to be thread safe, 
    // _config should only be initialized once
    const TrkQualCatalogConfig _config;
  };
}

#endif
