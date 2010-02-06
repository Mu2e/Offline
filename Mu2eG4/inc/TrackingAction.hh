#ifndef TrackingAction_H
#define TrackingAction_H 1
//
// Steering routine for user tracking actions. 
// If Mu2e needs many different user tracking actions, they
// should be called from this class.
//
// $Id: TrackingAction.hh,v 1.1 2010/02/06 19:39:09 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/06 19:39:09 $
//
// Original author Rob Kutschke
//

#include <string>

#include "Mu2eG4/inc/EventNumberList.hh"

#include "globals.hh"
#include "G4UserTrackingAction.hh"

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class TrackingAction: public G4UserTrackingAction{

  public:
    TrackingAction( const SimpleConfig& config);
    virtual ~TrackingAction();

    virtual void PreUserTrackingAction(const G4Track* trk);
    virtual void PostUserTrackingAction(const G4Track* trk);

  private:

    // Count number of calls to ClassifyNewTrack this event.
    int ncalls;
    EventNumberList debugList;

    void printInfo(const G4Track* trk, const std::string& text);

  };

} // end namespace mu2e

#endif

