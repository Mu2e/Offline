#ifndef TrkDiag_TrkQualCatalogMaker_hh
#define TrkDiag_TrkQualCatalogMaker_hh

//
// Makes the TrkQualCatalog ProditionsEntitiy
//

#include "TrkDiag/inc/TrkQualCatalog.hh"
#include "TrkDiag/inc/TrkQualCatalogConfig.hh"
#include "DbTables/inc/TrkQualDb.hh"

namespace mu2e {
  class TrkQualCatalogMaker {
  public:
    TrkQualCatalogMaker(TrkQualCatalogConfig const& config):_config(config) {}
    TrkQualCatalog::ptr_t fromFcl();
    TrkQualCatalog::ptr_t fromDb(TrkQualDb::cptr_t tqDb);
  private:
    // this object needs to be thread safe, 
    // _config should only be initialized once
    const TrkQualCatalogConfig _config;
  };
}

#endif
