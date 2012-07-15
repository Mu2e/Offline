//
// $Id: ExtMonFNAL_SD.cc,v 1.7 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//

#include <iostream>

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4ios.hh"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ExtMonFNAL_SD.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

//#define AGDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

using namespace std;

namespace mu2e {

  ExtMonFNAL_SD::ExtMonFNAL_SD(G4String const name, SimpleConfig const & config )
    : Mu2eSensitiveDetector(name,config)
  {
    verboseLevel = config.getInt("extMonFNAL.sd.verbosityLevel");
  }

  G4bool ExtMonFNAL_SD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    // Calculate energy deposition
    G4double edep = aStep->GetTotalEnergyDeposit();

    AGDEBUG("AG: edep = "<<edep);

    if( edep<=0 ) return false;

    _currentSize += 1;

    if ( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in "
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    // These are in G4 world coordinates
    const G4ThreeVector prePosG4 = aStep->GetPreStepPoint()->GetPosition();
    const G4ThreeVector preMomG4 = aStep->GetPreStepPoint()->GetMomentum();
    AGDEBUG("AG: got prePosG4     = "<<prePosG4<<    ", preMomG4     = "<<preMomG4);

    // Transform the global coordinates into the ExtMon frame.
    // First do G4->Mu2e by subtracting _mu2eOrigin (defined in the base class)
    // then convert Mu2e->ExtMon using the ExtMon class

    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    const G4ThreeVector prePosExtMon = extmon->mu2eToExtMon_position(prePosG4 - _mu2eOrigin);
    const G4ThreeVector preMomExtMon = extmon->mu2eToExtMon_momentum(preMomG4);
    AGDEBUG("AG: got prePosExtMon = "<<prePosExtMon<<", preMomExtMon = "<<preMomExtMon);

    // Which process caused this step to end?
    const G4String& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    AGDEBUG("pname = "<<pname);
    const ProcessCode endCode(_processInfo->findAndCount(pname));

    // Add the hit to the framework collection.
    _collection->
      push_back(StepPointMC(art::Ptr<SimParticle>( *_simID, aStep->GetTrack()->GetTrackID(), _event->productGetter(*_simID) ),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                            edep,
                            aStep->GetNonIonizingEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            prePosExtMon,
                            preMomExtMon,
                            aStep->GetStepLength(),
                            endCode
                            ));

    return true;

  }

} //namespace mu2e
