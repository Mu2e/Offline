#include "TrkDiag/inc/TrkQualCatalogMaker.hh"

namespace mu2e {

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fromFcl() {
    auto ptr = std::make_shared<TrkQualCatalog>();//_config.deadTimeAnalog(), 
    return ptr;
  }

  TrkQualCatalog::ptr_t TrkQualCatalogMaker::fromDb(TrkQualDb::cptr_t tqDb) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // now overwrite with db values
    //tqDb....

    return ptr;
  }

};
