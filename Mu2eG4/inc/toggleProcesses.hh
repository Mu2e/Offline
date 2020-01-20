#ifndef Mu2eG4_toggleProcesses_hh
#define Mu2eG4_toggleProcesses_hh
//
// Function that handles the switching on and off of G4 processes.  This
// is handled through the configuration files and includes the following
// commands:
//
// g4.noDecay - turns off decays of specified particles
//
//-----------------------------------------------------------------------------

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace mu2e{

  void switchDecayOff( const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

  void switchCaptureDModel( const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

  void addUserProcesses( const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

}  // end namespace mu2e

#endif /* Mu2eG4_toggleProcesses_hh */
