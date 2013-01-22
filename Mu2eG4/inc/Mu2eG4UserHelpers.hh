#ifndef Mu2eG4_Mu2eG4UserHelpers_hh
#define Mu2eG4_Mu2eG4UserHelpers_hh
//
// A collection of Geant4 user helper functions
// initially extracted from the TrackingAction
// $Id: Mu2eG4UserHelpers.hh,v 1.4 2013/01/22 21:42:06 mjlee Exp $
// $Author: mjlee $
// $Date: 2013/01/22 21:42:06 $
//
// Original author KLG based on Rob's TrackingAction
//


// C++ includes
#include <map>
#include <string>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "art/Utilities/CPUTimer.h"

// G4 includes
#include "globals.hh"

// Mu2e includes

#include "MCDataProducts/inc/SimParticleCollection.hh"

class G4Track;

namespace mu2e {

  namespace Mu2eG4UserHelpers {

    typedef SimParticleCollection::key_type    key_type;
    typedef SimParticleCollection::mapped_type mapped_type;
    typedef std::map<key_type,mapped_type>     map_type;

    // Check consistency of mother-daughter pointers.
    bool checkCrossReferences( bool doPrint, bool doThrow, map_type const& transientMap);

    // Debug printout.
    void printTrackInfo(G4Track const* trk, std::string const& text,
                        map_type const& transientMap,
                        art::CPUTimer const& timer,
                        CLHEP::Hep3Vector const& mu2eOrigin,
                        bool isEnd=false);

    G4String findStoppingProcessName(G4Track const* trk);
    ProcessCode findCreationCode(G4Track const* trk);

    double getPreLastStepKE(G4Track const* trk);

    int  getNSteps(G4Track const* trk);  

    // Control the saving of trajectories.
    // The first method does the big picture bookkeeping.
    // The second method decides yes/no for storing the trajectory of one track.
    void controlTrajectorySaving(G4Track const* trk, int sizeLimit, int currentSize, double pointTrajectoryMomentumCut = 50.);
    bool saveThisTrajectory( const G4Track* trk, double pointTrajectoryMomentumCut );

  }

}

#endif /* Mu2eG4_Mu2eG4UserHelpers_hh */
