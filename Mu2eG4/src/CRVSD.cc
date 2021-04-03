//
// Sensitive detector for CRV
//
// Added to fill the new visible energy deposit used in scintillators
//
// Original author KLG
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Mu2eG4/inc/CRVSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4Step.hh"

using namespace std;

namespace mu2e {

  CRVSD::CRVSD(G4String name, SimpleConfig const & config ):
    Mu2eG4SensitiveDetector(name,config)
  {}

  G4bool CRVSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if ( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of steps reached in "
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    ProcessCode endCode(_processInfo->
                findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    // VisibleEnergyDepositition suggested by Ralf E

    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetNonIonizingEnergyDeposit(),
                            G4LossTableManager::Instance()->EmSaturation()->
                            VisibleEnergyDeposition(aStep->GetTrack()->GetParticleDefinition(),
                                                    aStep->GetTrack()->GetMaterialCutsCouple(),
                                                    aStep->GetStepLength(),
                                                    aStep->GetTotalEnergyDeposit(),
                                                    aStep->GetNonIonizingEnergyDeposit()),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPostStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetPostStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));
      return true;

  }//ProcessHits

} //namespace mu2e
