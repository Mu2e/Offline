#ifndef GlobalConstantsService_GlobalConstantsHandle_hh
#define GlobalConstantsService_GlobalConstantsHandle_hh

#include <string>

#include "GlobalConstantsService/inc/GlobalConstantsService.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace mu2e {
  template <typename ENTITY>
  class GlobalConstantsHandle {
  public:

    explicit GlobalConstantsHandle() {

      art::ServiceHandle<GlobalConstantsService> sg;

      // Key and version are being ignored by the implementation.
      // Not exposing them via handle constructor args makes it easier
      // to give them meaning in the future, if needed.
      const std::string& key="";
      const std::string& version="";

      entity_ = sg->getElement<ENTITY>(key,version);
    }

    const ENTITY* operator->() const { return entity_;}
    const ENTITY& operator*()  const { return *entity_;}

  private:
    const ENTITY* entity_;
  };
}

#endif /* GlobalConstantsService_GlobalConstantsHandle_hh */
