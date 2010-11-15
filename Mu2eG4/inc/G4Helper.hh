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
// $Id: G4Helper.hh,v 1.1 2010/11/15 23:20:06 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/11/15 23:20:06 $
//
// Original author Rob Kutschke
//

// Framework include files
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Mu2eG4/inc/AntiLeakRegistry.hh"

namespace mu2e {

  class G4Helper {
  public:
    G4Helper(const edm::ParameterSet&, edm::ActivityRegistry&);
    ~G4Helper();
    
    AntiLeakRegistry& antiLeakRegistry(){ return _antiLeakRegistry; }
    
  private:

    AntiLeakRegistry _antiLeakRegistry;

  };

}

#endif
