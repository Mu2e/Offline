#ifndef BtrkHelper_BtrkHelper_hh
#define BtrkHelper_BtrkHelper_hh

//
// The btrk code needs some information and tools that are supplied externally.
// This service is the agent that supplies the that information and those tools.
//
// Original author Rob Kutschke
//

#include "BtrkHelper/inc/FileFinder.hh"
#include "BtrkHelper/inc/ParticleInfo.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include <string>
#include <memory>

namespace mu2e {

  class BtrkHelper {

public:
    BtrkHelper(const fhicl::ParameterSet&, art::ActivityRegistry&);

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

DECLARE_ART_SERVICE(mu2e::BtrkHelper, LEGACY)
#endif /* BtrkHelper_BtrkHelper_hh */
