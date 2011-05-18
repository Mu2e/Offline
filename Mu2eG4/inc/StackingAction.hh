#ifndef Mu2eG4_StackingAction_hh
#define Mu2eG4_StackingAction_hh
//
// Steering routine for user stacking actions.
// If Mu2e needs many different user stacking actions, they
// should be called from this class.
//
// $Id: StackingAction.hh,v 1.14 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//
// Notes
// 1) There is an option to kill particles that have a kinetic energy
//    below a given value.  In the long run we want to properly tune
//    physics lists.  The option provided here is a sledgehammer for
//    use when we need it during development.
//

#include <vector>

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4VPhysicalVolume;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class StackingAction: public G4UserStackingAction{

  public:
    StackingAction( const SimpleConfig& config);
    ~StackingAction();

  public:

    // Specified by G4UserStackingAction.
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    void NewStage();
    void PrepareNewEvent();

    // Methods specific to Mu2e.
    void beginRun( double dirtG4YMin, double dirtG4YMax );

  private:

    // Count number of calls to ClassifyNewTrack this event.
    int _ncalls;

    // Count events to limit printout.
    int _nevents;

    // Do we run the cosmic killer?
    bool _doCosmicKiller;
    int  _killLevel;
    double _cosmicpcut;
    double _yaboveDirtYmin;

    // Only stack primary particles.
    bool _primaryOnly;

    // Enable the code to kill particles with kinetic energy below threshold.
    // See Note 1).
    bool _killLowKineticEnergy;
    double _eKineMin;

    // Only one of these may be non-empty.
    // List of particles to remove (others will be kept).
    // List of particles to keep (others will be dropped).
    std::vector<int> _pdgToDrop;
    std::vector<int> _pdgToKeep;

    // Pointers to some physical volumes of interest.
    G4VPhysicalVolume * _dirtBodyPhysVol;
    G4VPhysicalVolume * _dirtCapPhysVol;

    // Y limits of the dirt overburden.
    double _dirtG4Ymin, _dirtG4Ymax;

    // Check to see if we kill low p tracks in the dirt or other shielding.
    bool cosmicKiller( const G4Track* aTrack);

  };

} // end namespace mu2e

#endif /* Mu2eG4_StackingAction_hh */
