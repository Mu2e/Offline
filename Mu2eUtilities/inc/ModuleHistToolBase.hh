#ifndef __Mu2eUtilities_ModuleHistToolBase_hh__
#define __Mu2eUtilities_ModuleHistToolBase_hh__

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

namespace art { class TFileService; }
// #include "art/Framework/Principal/Handle.h"

// #include "art/Framework/Principal/Event.h"
// #include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  
  class ModuleHistToolBase {
  public:

    ModuleHistToolBase() noexcept = default ;
    virtual ~ModuleHistToolBase()  noexcept = default ;
    
    virtual int bookHistograms(art::ServiceHandle<art::TFileService> & Tfs) ; // = 0 ; 
    virtual int fillHistograms(void* Data, int Mode = -1) ;
    virtual int debug         (void* Data, int Mode = -1);
  };
}

#endif

