//
//
// A test service that causes an exception to be thrown during
// the begin run callback.
//
// Contact person Rob Kutschke
//

#include "ConfigTools/inc/SimpleConfig.hh"
#include "Sandbox/inc/Bug01Service.hh"

#include "art/Persistency/Provenance/RunID.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <iostream>
#include <fstream>

namespace mu2e {

  Bug01Service::Bug01Service(fhicl::ParameterSet const& pset,
                             art::ActivityRegistry&     registry):
    fileName_(pset.get<std::string>("fileName")){
    std::cout << "Bug01Service::constructor" << std::endl;
    registry.sPreBeginRun.watch(this, &Bug01Service::preBeginRun);
  }

  Bug01Service::~Bug01Service(){
    std::cout << "Bug01Service::destructor" << std::endl;
  }

  void  Bug01Service::preBeginRun(art::Run const& run) {
    std::cout << "Bug01Service::preBeginRun: "
         << run.id() << " "
         << fileName_
         << std::endl;

    SimpleConfig config( fileName_ );
    config.print();
  }

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::Bug01Service);
