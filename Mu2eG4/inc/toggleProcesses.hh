#ifndef Mu2eG4_toggleProcesses_hh
#define Mu2eG4_toggleProcesses_hh
//
// Function that handles the switching on and off of G4 processes.  This
// is handled through the configuration files and includes the following 
// commands:
//
// g4.noDecay - turns off decays of specified particles
// muMinusConversionAtRest.do - turns on the at rest G4 process 
// MuonMinusConversionAtRest and turns off MuonMinusCaptureAtRest
//
// $Id: toggleProcesses.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $ 
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
//-----------------------------------------------------------------------------

namespace mu2e{

  // Specializations for particular trackers.  Called by the public entry point.
  void switchDecayOff( const SimpleConfig& config  );
  void addUserProcesses( const SimpleConfig& config  );

}  // end namespace mu2e

#endif /* Mu2eG4_toggleProcesses_hh */


