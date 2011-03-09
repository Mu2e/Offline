#ifndef StackingAction_H
#define StackingAction_H 1
//
// Steering routine for user stacking actions. 
// If Mu2e needs many different user stacking actions, they
// should be called from this class.
//
// $Id: StackingAction.hh,v 1.6 2011/03/09 21:43:35 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/03/09 21:43:35 $
//
// Original author Rob Kutschke
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
    int ncalls;

    // Count events to limit printout.
    int nevents;

    // Do we run the cosmic killer?
    bool doCosmicKiller;
    int  killLevel;

    // List of particles to remove (others will be kept)
    std::vector<int> _pdgToDrop;

    // Pointers to some physical volumes of interest.
    G4VPhysicalVolume * dirtBodyPhysVol;
    G4VPhysicalVolume * dirtCapPhysVol;

    // Y limits of the dirt overburden.
    double _dirtG4Ymin, _dirtG4Ymax;

    // Check to see if we kill low p tracks in the dirt or other shielding.
    bool cosmicKiller( const G4Track* aTrack);

    // Drop tracks from a list of PDG Id's.
    bool dropByPDGId( G4Track const *);

  };

} // end namespace mu2e

#endif

