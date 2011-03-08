#ifndef ToggleProcesses_HH
#define ToggleProcesses_HH
//
// Function that handles the switching on and off of G4 processes.  This
// is handled through the configuration files and includes the following 
// commands:
//
// g4.noDecay - turns off decays of specified particles
// muMinusConversionAtRest.do - turns on the at rest G4 process 
// MuonMinusConversionAtRest and turns off MuonMinusCaptureAtRest
//
// $Id: toggleProcesses.hh,v 1.1 2011/03/08 14:15:00 ayarritu Exp $ 
// $Author: ayarritu $
// $Date: 2011/03/08 14:15:00 $
//
//-----------------------------------------------------------------------------

namespace mu2e{

  // Specializations for particular trackers.  Called by the public entry point.
  void switchDecayOff( const SimpleConfig& config  );
  void addUserProcesses( const SimpleConfig& config  );

}  // end namespace mu2e

#endif


