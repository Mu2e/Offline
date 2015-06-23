#ifndef BTrkHelper_BTrkHelper_hh
#define BTrkHelper_BTrkHelper_hh

//
// The btrk code needs some information and tools that are supplied externally.
// This service is the agent that supplies the that information and those tools.
//
// Original author Rob Kutschke
//

#include "BTrkHelper/inc/FileFinder.hh"
#include "BTrkHelper/inc/ParticleInfo.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include <string>
#include <memory>

namespace mu2e {

  class BTrkHelper {

public:
    BTrkHelper(const fhicl::ParameterSet&, art::ActivityRegistry&);

    // Functions registered for callbacks.
    void beginRun   ( art::Run    const &run    );
    void beginSubRun( art::SubRun const &subRun );

private:

    // Access to particle data table via TrkParticle::type
    ParticleInfo particleInfo_;

    // Tool to resolve filenames from FHICL_FILE_PATH
    FileFinder   fileFinder_;

  };

}

DECLARE_ART_SERVICE(mu2e::BTrkHelper, LEGACY)
#endif /* BTrkHelper_BTrkHelper_hh */
