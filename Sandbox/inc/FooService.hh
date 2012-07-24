#ifndef Sandbox_FooService_hh
#define Sandbox_FooService_hh

//
// A test service that prints tracer ouptut. Used to study order of calls.
//
// $Id: FooService.hh,v 1.1 2012/07/24 20:00:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/24 20:00:28 $
//
// Contact person Rob Kutschke
//

// Framework include files
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Run.h"

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

#endif /* Sandbox_FooService_hh */
