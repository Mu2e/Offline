#ifndef SteppingAction_h
#define SteppingAction_h 1
//
// Called at every G4 step.
//
// $Id: SteppingAction.hh,v 1.7 2010/12/11 00:42:51 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/12/11 00:42:51 $
//
// Original author Rob Kutschke
//
#include <vector>

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"

// G4 includes
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"

// Forward declarations outside of mu2e namespace.
class G4VPhysicalVolume;
class G4Track;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class SteppingAction : public G4UserSteppingAction{

  public:
    SteppingAction( const SimpleConfig& config );
    ~SteppingAction(){};
    
    void UserSteppingAction(const G4Step*);

    void BeginOfEvent();
    void EndOfEvent();

    void BeginOfTrack();
    void EndOfTrack();

    int nKilledStepLimit() const { return _nKilledStepLimit; }

    // Called by G4_plugin.
    void beginRun();
    
    G4ThreeVector const& lastPosition() const { return _lastPosition; }
    G4ThreeVector const& lastMomentum() const { return _lastMomentum; }

    void setZRef( G4double zref){ 
      _zref=zref;
    }
  
  private:

    // Start: information from the run time configuration.

    // Which killers will be enabled?
    bool _doKillLowEKine;
    bool _doKillInHallAir;
    bool _killerVerbose;

    // Minimum energy cut.
    double _eKineMin;

    // List of particles to remove (others will be kept)
    std::vector<int> _pdgToDrop;
    // Maximum allowed number of steps per event
    int _maxSteps;
    int _nSteps;
    int _nKilledStepLimit;

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugEventList;
    EventNumberList _debugTrackList;

    // End: information from the run time configuration.

    // Address of the physical volume object for the hall air.
    G4VPhysicalVolume* _hallAirPhysVol;

    // Used in the turn around code.
    G4ThreeVector _lastPosition;
    G4ThreeVector _lastMomentum;
    G4double      _zref;

    // Two functions to decide whether or not to kill tracks.
    bool killLowEKine ( const G4Track* );
    bool killInHallAir( const G4Track* );

  };
  
} // end namespace mu2e
#endif
