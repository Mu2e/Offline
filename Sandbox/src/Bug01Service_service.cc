//
//
// A test service that throws an exception to be thrown during
// the begin run callback.
//
// Contact person Rob Kutschke
//

#include "ConfigTools/inc/SimpleConfig.hh"
#include "Sandbox/inc/Bug01Service.hh"

#include "canvas/Persistency/Provenance/RunID.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "cetlib_except/exception.h"

#include <iostream>

mu2e::Bug01Service::Bug01Service(fhicl::ParameterSet const& pset,
                                 art::ActivityRegistry&     registry){
  std::cout << "Bug01Service::constructor" << std::endl;
  registry.sPreBeginRun.watch(this, &Bug01Service::preBeginRun);
}

mu2e::Bug01Service::~Bug01Service(){}

void  mu2e::Bug01Service::preBeginRun(art::Run const& ) {
  std::cout << "Bug01Service::preBeginRun: throwing now ... " << std::endl;
  throw cet::exception("FOO") << "Throwing an exception, just because ... "<< "\n";
}

DEFINE_ART_SERVICE(mu2e::Bug01Service);
