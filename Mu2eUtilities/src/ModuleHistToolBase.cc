///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

// #include "art/Framework/Services/Registry/ServiceHandle.h"
// #include "art/Framework/Services/Optional/TFileService.h"
// #include "art/Framework/Principal/Handle.h"

#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"

namespace art { class TFileService; }

namespace mu2e {

  // ModuleHist::ModuleHist() {
  // }

  // ModuleHist::~ModuleHist() {
  // }
//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  int ModuleHistToolBase::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
    return 0;
  }


//-----------------------------------------------------------------------------
//  fillHistograms keeps full knowledge of what and how has to be filled
//-----------------------------------------------------------------------------
  int ModuleHistToolBase::fillHistograms(void* Data, int Module) {
    return 0;
  }

//-----------------------------------------------------------------------------
//  hook for debug printout
//-----------------------------------------------------------------------------
  int ModuleHistToolBase::debug(void* Data, int Mode) {
    return 0;
  }

}
