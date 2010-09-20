#ifndef TrackingAction_H
#define TrackingAction_H 1
//
// Steering routine for user tracking actions. 
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.hh,v 1.4 2010/09/20 02:57:05 logash Exp $
// $Author: logash $
// $Date: 2010/09/20 02:57:05 $
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

  class TrackingAction: public G4UserTrackingAction{

  public:

    TrackingAction( const SimpleConfig& config);
    virtual ~TrackingAction();

    // These methods are required by G4
    virtual void PreUserTrackingAction (const G4Track* trk);
    virtual void PostUserTrackingAction(const G4Track* trk);

    // All methods after here are Mu2e specific.

    // Do all things that need to be done at the beginning/end of an event.
    void endEvent( SimParticleCollection& simParticles );
    void beginEvent();

    // Record start and end points of each track created by G4.
    // Copy to the event data.
    void saveSimParticleStart(const G4Track* trk);
    void saveSimParticleEnd  (const G4Track* trk);
    void saveSimParticleCopy (SimParticleCollection& simParticles);

    // Receive persistent volume information and save it for the duration of the run.
    void beginRun( const PhysicalVolumeHelper& physVolHelper, 
                   CLHEP::Hep3Vector const& mu2eOrigin ){
      _physVolHelper = &physVolHelper;
      _mu2eOrigin    =  mu2eOrigin;
    }

    static void setSizeLimit(int sizeLimit) {
      _sizeLimit = sizeLimit;
    }


  private:

    // Lists of events and tracks for which to enable debug printout.
    EventNumberList _debugList;

    // Transient information collected during the event.
    // Tracks are not processed in order, hence the map.
    // Used by saveSimParticle*
    std::map<uint32_t,SimParticle> _spmap;

    // Utility to translate between transient and persistent representations.
    const PhysicalVolumeHelper* _physVolHelper;

    CLHEP::Hep3Vector _mu2eOrigin;

    edm::CPUTimer _timer;

    // Debug printout.
    void printInfo(const G4Track* trk, const std::string& text, bool isEnd=false);

    // Limit maximum size of the steps collection
    static int _sizeLimit;
    int _currentSize;

  };

} // end namespace mu2e

#endif

