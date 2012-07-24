//
// A test service that prints tracer ouptut. Used to study order of calls.
// The BarSerivce depends on the FooService.
//
// $Id: BarService_service.cc,v 1.1 2012/07/24 20:00:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/24 20:00:28 $
//
// Contact person Rob Kutschke
//

// C++ include files
#include <iostream>

#include "Sandbox/inc/BarService.hh"
#include "Sandbox/inc/FooService.hh"

// Framework include files
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

using namespace std;

namespace mu2e {

  BarService::BarService(fhicl::ParameterSet const& pset,
                                   art::ActivityRegistry&iRegistry){
    cout << "BarService:: constructor: " << endl;
    iRegistry.watchPreBeginRun(this, &BarService::preBeginRun);

    // This fails since foo is constructed last.
    //art::ServiceHandle<FooService> foo;
    //foo->poke();
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
