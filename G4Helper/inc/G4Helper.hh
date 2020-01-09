#ifndef G4Helper_G4Helper_hh
#define G4Helper_G4Helper_hh
//
// The design of G4 requires that users new many objects and then delete
// them at the appropriate time, usually the end of the G4 run.  This
// Service exists to manage the delete automatically.  It is also available
// as a place to create any other required singleton-like behaviour for
// support of G4.  For technical reasons, this cannot be done by making
// Mu2eG4RunManager a singleton.
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "G4Helper/inc/VolumeInfo.hh"

// Framework include files
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "boost/regex.hpp"

// C++ includes
#include <map>
#include <string>
#include <vector>

namespace mu2e {

  class G4Helper {
  public:
    G4Helper(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~G4Helper();

    AntiLeakRegistry& antiLeakRegistry(){ return _antiLeakRegistry; }

    // Versions of the map [] operator that check for errors.
    VolumeInfo& locateVolInfo( const std::string key);
    void addVolInfo( const VolumeInfo& info );

    // Find all VolumeInfo objects whose name matches a regex.
    std::vector<VolumeInfo const*> locateVolInfo( boost::regex const& re ) const;

  private:

    AntiLeakRegistry _antiLeakRegistry;

    std::map<std::string,VolumeInfo> _volumeInfoList;

  };

}

DECLARE_ART_SERVICE(mu2e::G4Helper, SHARED)
#endif /* G4Helper_G4Helper_hh */
