//
// A test service that prints tracer ouptut. Used to study order of calls.
// The XBarSerivce depends on the FooService.  The XBarService is just
// a copy of the BarService but it comes after Foo in an alphabetic sort.
//
// $Id: XBarService_service.cc,v 1.1 2012/07/24 20:00:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/24 20:00:28 $
//
// Contact person Rob Kutschke
//

// C++ include files
#include <iostream>

#include "Sandbox/inc/XBarService.hh"
#include "Sandbox/inc/FooService.hh"

// Framework include files
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

using namespace std;

namespace mu2e {

  XBarService::XBarService(fhicl::ParameterSet const& pset,
                           art::ActivityRegistry&iRegistry){
    cout << "XBarService:: constructor: " << endl;
    iRegistry.watchPreBeginRun(this, &XBarService::preBeginRun);

    if ( pset.get<bool>("pokeFoo",false) ) {
      art::ServiceHandle<FooService> foo;
      cout << "     ";
      foo->poke();
    }
    cout << "XBarService:: done constructor: " << endl;
  }

  XBarService::~XBarService(){
    cout << "XBarService::destructor" << endl;
  }

  void  XBarService::preBeginRun(art::Run const& run) {
    cout << "XBarService::preBeginRun: " << run.id() << endl;
    art::ServiceHandle<FooService> foo;
    cout << "     ";
    foo->poke();
    cout << "XBarService::preBeginRun done " << endl;
  }

  void XBarService::poke() const{
    cout << "XBarService::you poked me." << endl;
    art::ServiceHandle<FooService> foo;
    cout << "     ";
    foo->poke();
    cout << "XBarService::you poked me. done" << endl;
  }


} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::XBarService);
