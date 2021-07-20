#ifndef Sandbox_FooService_hh
#define Sandbox_FooService_hh

//
// A test service that prints tracer ouptut. Used to study order of calls.
//
//
// Contact person Rob Kutschke
//

// Framework include files
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

namespace mu2e {

  class FooService {
  public:
    FooService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~FooService();

    void preBeginRun( art::Run const &run);

    void poke() const;

  private:

  };

}

DECLARE_ART_SERVICE(mu2e::FooService, LEGACY)
#endif /* Sandbox_FooService_hh */
