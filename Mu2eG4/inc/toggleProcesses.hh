#ifndef Mu2eG4_toggleProcesses_hh
#define Mu2eG4_toggleProcesses_hh
//
// Function that handles the switching on and off of G4 processes.  This
// is handled through the configuration files
//
//-----------------------------------------------------------------------------

#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"

namespace mu2e{

  void switchDecayOff( const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

  void switchCaptureDModel( const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

  void addUserProcesses( const Mu2eG4Config::Physics& phys,
                         const Mu2eG4Config::Debug& debug,
                         const Mu2eG4ResourceLimits& lim
                         );

}  // end namespace mu2e

#endif /* Mu2eG4_toggleProcesses_hh */
