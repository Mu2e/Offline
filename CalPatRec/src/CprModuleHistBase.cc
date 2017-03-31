///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

#include "CalPatRec/inc/CprModuleHistBase.hh"

namespace mu2e {


  // ModuleHist::ModuleHist() {
  // }

  // ModuleHist::~ModuleHist() {
  // }
//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  int CprModuleHistBase::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs, TObject* Hist) {
    return 0;
  }


//-----------------------------------------------------------------------------
// Mode = 0: event-level histograms
//      = 1: track-level histograms ... etc
//-----------------------------------------------------------------------------
  int CprModuleHistBase::fillHistograms(int Mode, const TObject* Data, TObject* Hist) {
    return 0;
  }

}
