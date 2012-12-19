#ifndef Mu2eG4_Mu2eG4UserHelpers_hh
#define Mu2eG4_Mu2eG4UserHelpers_hh
//
// A collection of Geant4 user helper functions
// initially extracted from the TrackingAction
// $Id: Mu2eG4UserHelpers.hh,v 1.1 2012/12/19 23:36:57 genser Exp $
// $Author: genser $
// $Date: 2012/12/19 23:36:57 $
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

    G4String findStoppingProcessName(G4Track const* trk);

    bool checkCrossReferences( bool doPrint, bool doThrow, map_type const& _transientMap);

    void printTrackInfo(G4Track const* trk, std::string const& text,
                        map_type const& _transientMap,
                        art::CPUTimer const& _timer,
                        CLHEP::Hep3Vector const& _mu2eOrigin,
                        bool isEnd=false);

    ProcessCode findCreationCode(G4Track const* trk);

    double getPreLastStepKE(G4Track const* trk);

    int  getNSteps(G4Track const* trk);  

    void controlTrajectorySaving(G4Track const* trk, int _sizeLimit, int _currentSize);

    bool saveThisTrajectory( const G4Track* trk );

  }

}

#endif /* Mu2eG4_Mu2eG4UserHelpers_hh */
