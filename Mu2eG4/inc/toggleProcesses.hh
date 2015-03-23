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
// $Id: toggleProcesses.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
//-----------------------------------------------------------------------------

namespace fhicl { class ParameterSet; }

namespace mu2e{
  class SimpleConfig;

  // Specializations for particular trackers.  Called by the public entry point.
  void switchDecayOff( const SimpleConfig& config  );
  void addUserProcesses( const SimpleConfig& config  );

  void switchDecayOff(const fhicl::ParameterSet& pset);
  void addUserProcesses(const fhicl::ParameterSet& pset);

}  // end namespace mu2e

#endif /* Mu2eG4_toggleProcesses_hh */
