#ifndef TrackingAction_H
#define TrackingAction_H 1
//
// Steering routine for user tracking actions. 
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.hh,v 1.6 2010/11/08 23:52:33 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/08 23:52:33 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <string>
#include <map>

// Framework includes
#include "FWCore/Utilities/interface/CPUTimer.h"

// Mu2e includes
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "ToyDP/inc/SimParticleCollection.hh"

// G4 includes.
#include "G4UserTrackingAction.hh"

// Other includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class SteppingAction;
  //class SimParticleCollection;

  class TrackingAction: public G4UserTrackingAction{

  public:

    TrackingAction( const SimpleConfig& config, SteppingAction *);
    virtual ~TrackingAction();

    // These methods are required by G4
    virtual void PreUserTrackingAction (const G4Track* trk);
    virtual void PostUserTrackingAction(const G4Track* trk);

    // All methods after here are Mu2e specific.

    // Do all things that need to be done at the beginning/end of an event.
    void beginEvent( SimParticleCollection& simParticles );
    void endEvent();

    // Record start and end points of each track created by G4.
    // Copy to the event data.
    void saveSimParticleStart(const G4Track* trk);
    void saveSimParticleEnd  (const G4Track* trk);

    // Receive persistent volume information and save it for the duration of the run.
    void beginRun( const PhysicalVolumeHelper& physVolHelper, 
                   CLHEP::Hep3Vector const& mu2eOrigin ){
      _physVolHelper = &physVolHelper;
      _mu2eOrigin    =  mu2eOrigin;
    }

    // Check consistency of mother-daughter pointers.
    bool checkCrossReferences( bool doPrint, bool doThrow);

  private:

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugList;

    // Non-owning pointer to the collection of information about SimParticles.
    // G4_plugin owns the collection but this class is responsible for filling it.
    // G4_plugin sets this pointer at the start of each event.
    SimParticleCollection* _spmap;

    // Utility to translate between transient and persistent representations.
    const PhysicalVolumeHelper* _physVolHelper;

    CLHEP::Hep3Vector _mu2eOrigin;

    edm::CPUTimer _timer;

    // Debug printout.
    void printInfo(const G4Track* trk, const std::string& text, bool isEnd=false);

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // Non-owning pointer to stepping action.
    SteppingAction * _stepping; 

  };

} // end namespace mu2e

#endif

