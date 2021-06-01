#ifndef Sandbox_Bug01Service_hh
#define Sandbox_Bug01Service_hh

//
// A test service that causes an exception to be thrown during
// the begin run callback.
//
// Contact person Rob Kutschke
//

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace mu2e {

  class Bug01Service {
  public:
    Bug01Service(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~Bug01Service();

    void preBeginRun( art::Run const &run);

  };

}

DECLARE_ART_SERVICE(mu2e::Bug01Service, LEGACY)
#endif /* Sandbox_Bug01Service_hh */
