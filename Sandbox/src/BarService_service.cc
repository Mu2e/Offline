//
// A test service that prints tracer ouptut. Used to study order of calls when
// one service uses another service.
//
//
// Contact person Rob Kutschke
//
// Notes:
// 1) This works with art v1_00_11 and greater.
//
// 2) The point of this example is to show how to control the order
//    that art calls the callbacks of services.  For example, at beginRun,
//    will art call FooService::beginRun before or after it calls
//    BarService:: beginRun?  This example forces art to, at each begin run,
//    call FooService::beginRun before calling BarService::beginRun.
//    This mocks up the following idea: if the GeometryService (Bar) wants
//    to produce a fully aligned geometry, then we need to make sure that
//    ConditionsService::beginRun (Foo) gets called
//    before GeometryService::beginRun.
//
// 3) The rule used by art is that callbacks are executed in the order in
//    which they are registered in the registry.  In this example, the
//    call to register the callback BarService::beginRun is done AFTER getting
//    a handle to FooService. Suppose that FooService were already
//    instantiated before the call to the BarService c'tor ; in this case
//    Foo's callbacks would already be registered and they would be called
//    first.  The other alternative is that FooService is not yet instantiated
//    at the call to the c'tor of BarService. In this case, the action of
//    getting the handle to FooService will trigger its c'tor; this will
//    cause its callbacks to be registered. Only after this, does this
//    example register the callbacks for BarService. Therefore we have
//    enforced that Foo's callbacks are called before those of Bar.
//
// 4) It is possible to specify a long chain of services that depend on other services.
//    this technique will recursively to ensure that all callbacks are done in the correct order.
//
// 5) It is not possible to have both ServiceA call ServiceB and ServiceB to call ServiceA; the
//    dependency graph must be acyclic.
//

// C++ include files
#include <iostream>

// Mu2e include files.
#include "Sandbox/inc/BarService.hh"
#include "Sandbox/inc/FooService.hh"

// Framework include files
#include "canvas/Persistency/Provenance/RunID.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

using namespace std;

namespace mu2e {

  BarService::BarService(fhicl::ParameterSet const& pset,
                           art::ActivityRegistry&iRegistry){
    bool doPoke(pset.get<bool>("pokeFoo"));

    cout << "BarService:: constructor: " << doPoke << endl;

    if ( doPoke ) {
      art::ServiceHandle<FooService> foo;
      cout << "     ";
      foo->poke();
    }

    // Magic alert: see notes 2 and 3.
    iRegistry.sPreBeginRun.watch(this, &BarService::preBeginRun);
    cout << "BarService:: done constructor: " << endl;
  }

  BarService::~BarService(){
    cout << "BarService::destructor" << endl;
  }

  void  BarService::preBeginRun(art::Run const& run) {
    cout << "BarService::preBeginRun: " << run.id() << endl;
    art::ServiceHandle<FooService> foo;
    cout << "     ";
    foo->poke();
    cout << "BarService::preBeginRun done " << endl;
  }

  void BarService::poke() const{
    cout << "BarService::you poked me." << endl;
    art::ServiceHandle<FooService> foo;
    cout << "     ";
    foo->poke();
    cout << "BarService::you poked me. done" << endl;
  }

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::BarService);
