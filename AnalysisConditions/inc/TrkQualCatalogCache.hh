#ifndef AnalysisConditions_TrkQualCatalogCache_hh
#define AnalysisConditions_TrkQualCatalogCache_hh

//
// TrkQualCatalogCache is the ProditionsCache for AnaTrkQualDb
//

// Mu2e includes
#include "AnalysisConditions/inc/MVACatalogCache.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "DbTables/inc/AnaTrkQualDb.hh"

namespace mu2e {

  typedef MVACatalogCache<TrkQual, AnaTrkQualDb> TrkQualCatalogCache;
}

#endif
