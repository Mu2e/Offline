#ifndef AnalysisConditions_TrkQualCatalogCache_hh
#define AnalysisConditions_TrkQualCatalogCache_hh

//
// TrkQualCatalogCache is the ProditionsCache for TrkQualDb
//

// Mu2e includes
#include "AnalysisConditions/inc/MVACatalogCache.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "DbTables/inc/TrkQualDb.hh"

namespace mu2e {

  typedef MVACatalogCache<TrkQual, TrkQualDb> TrkQualCatalogCache;
}

#endif
