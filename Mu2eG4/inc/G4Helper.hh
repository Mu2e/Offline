#ifndef Mu2eG4_G4Helper_hh
#define Mu2eG4_G4Helper_hh
//
// The design of G4 requires that users new many objects and then delete
// them at the appropriate time, usually the end of the G4 run.  This 
// Service exists to manage the delete automatically.  It is also available
// as a place to create any other required singleton-like behaviour for
// support of G4.  For technical reasons, this cannot be done by making 
// Mu2eG4RunManager a singleton.
//
// $Id: G4Helper.hh,v 1.2 2010/11/16 14:43:11 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/16 14:43:11 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <map>
#include <string>

// Framework include files
#include "FWCore/ServiceRegistry/interface/Service.h"

// Mu2e includes
#include "Mu2eG4/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/AntiLeakRegistry.hh"

namespace mu2e {

  class G4Helper {
  public:
    G4Helper(const edm::ParameterSet&, edm::ActivityRegistry&);
    ~G4Helper();
    
    AntiLeakRegistry& antiLeakRegistry(){ return _antiLeakRegistry; }

    // Versions of the map [] operator that check for errors.
    VolumeInfo& locateVolInfo( const std::string key);
    void addVolInfo( const VolumeInfo& info );
    
  private:

    AntiLeakRegistry _antiLeakRegistry;

    std::map<std::string,VolumeInfo> _volumeInfoList;

  };

}

#endif
