#ifndef AnalysisConditions_TrkQualCatalogCache_hh
#define AnalysisConditions_TrkQualCatalogCache_hh

//
// TrkQualCatalogCache is the ProditionsCache for AnaTrkQualDb
//

// Mu2e includes
#include "Offline/AnalysisConditions/inc/MVACatalogCache.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/DbTables/inc/AnaTrkQualDb.hh"

namespace mu2e {

  typedef MVACatalogCache<TrkQual, AnaTrkQualDb> TrkQualCatalogCache;
}

#endif
