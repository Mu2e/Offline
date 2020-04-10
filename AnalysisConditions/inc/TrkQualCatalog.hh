#ifndef AnalysisConditions_TrkQualCatalog_hh
#define AnalysisConditions_TrkQualCatalog_hh

//
// TrkQualCatalog is the ProditionsEntry for TrkQualDb
//

// Mu2e includes
#include "AnalysisConditions/inc/MVACatalog.hh"
#include "RecoDataProducts/inc/TrkQual.hh"

namespace mu2e {

  typedef MVAEntry<TrkQual> TrkQualEntry;
  typedef MVAEntries<TrkQual> TrkQualEntries;
  typedef MVACatalog<TrkQual> TrkQualCatalog;

  template<>
  MVACatalog<TrkQual>::MVACatalog() : _name("TrkQualCatalog") { }

  template<>
  MVACatalog<TrkQual>::MVACatalog(TrkQualEntries entries) : _name("TrkQualCatalog"), _entries(entries) { }
}

#endif
