#ifndef Mu2eG4Helper_Mu2eG4Helper_hh
#define Mu2eG4Helper_Mu2eG4Helper_hh
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
#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

// Framework include files
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"

#include "boost/regex.hpp"

// C++ includes
#include <map>
#include <string>
#include <vector>

namespace mu2e {

  class Mu2eG4Helper {
  public:
    Mu2eG4Helper(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~Mu2eG4Helper();

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

DECLARE_ART_SERVICE(mu2e::Mu2eG4Helper, SHARED)
#endif /* Mu2eG4Helper_Mu2eG4Helper_hh */
