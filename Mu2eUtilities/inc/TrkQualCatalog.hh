#ifndef Mu2eUtilities_TrkQualCatalog_hh
#define Mu2eUtilities_TrkQualCatalog_hh

//
// TrkQualCatalog is the ProditionsEntry for TrkQualDb
//

// Mu2e includes
#include "Mu2eUtilities/inc/MVACatalog.hh"
#include "RecoDataProducts/inc/TrkQual.hh"

namespace mu2e {

  typedef MVAEntry<TrkQual> TrkQualEntry;
  typedef MVACatalog<TrkQual> TrkQualCatalog;
}

#endif
