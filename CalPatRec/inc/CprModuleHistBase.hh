#ifndef __CalPatRec_ModuleHist_hh__
#define __CalPatRec_ModuleHist_hh__

#include "TObject.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  
  class CprModuleHistBase {
  public:

    CprModuleHistBase() noexcept = default ;
    virtual ~CprModuleHistBase()  noexcept = default ;
    
    virtual int bookHistograms(art::ServiceHandle<art::TFileService> & Tfs, TObject* Hist) ; // = 0 ; 
    virtual int fillHistograms(int Mode, const TObject* Data, TObject* Hist) ;
  };
}

#endif

