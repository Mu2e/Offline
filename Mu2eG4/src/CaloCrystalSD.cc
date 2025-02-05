//
// Define a sensitive detector for CaloCrystal Detectors
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/CaloCrystalSD.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Offline/Mu2eG4/inc/SimParticleHelper.hh"
#include "Offline/Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "Geant4/G4LossTableManager.hh"
#include "Geant4/G4Step.hh"


namespace mu2e {

  CaloCrystalSD::CaloCrystalSD(G4String name, SimpleConfig const & config ):
    Mu2eG4SensitiveDetector(name,config)
  { }

  G4bool CaloCrystalSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
  {
      G4double edep = aStep->GetTotalEnergyDeposit();
      if (edep < 1e-6) return false;


      _currentSize += 1;
      if( _sizeLimit>0 && _currentSize>_sizeLimit && (_currentSize - _sizeLimit)==1)
      {
          mf::LogWarning("G4") << "Maximum number of steps reached in "
                               << SensitiveDetectorName<< ": "<< _currentSize <<std::endl;
          return false;
      }


    ProcessCode endCode(_processInfo->findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    G4int copyNo = touchableHandle->GetCopyNumber(1);  // Make sure to get the right copy number level here
    G4ThreeVector posWorld = aStep->GetPreStepPoint()->GetPosition();

    G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
    G4ThreeVector posLocal           = toLocal.TransformPoint(posWorld);
    // diagnosis purposes only when playing with the geometry, uncomment next two line
    for (int i=0;i<=touchableHandle->GetHistoryDepth();++i) std::cout<<"Cry Transform level "<<i<<"   "<<touchableHandle->GetCopyNumber(i)
    <<"   "<<touchableHandle->GetHistory()->GetTransform(touchableHandle->GetHistoryDepth()-i).TransformPoint(posWorld)
    <<"  "<<touchableHandle->GetVolume(i)->GetName()<<std::endl;

    // VisibleEnergyDeposition following Birks law
    _collection->push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                                       copyNo,
                                       edep,
                                       aStep->GetNonIonizingEnergyDeposit(),
                                       G4LossTableManager::Instance()->EmSaturation()->
                                       VisibleEnergyDeposition(aStep->GetTrack()->GetParticleDefinition(),
                                                               aStep->GetTrack()->GetMaterialCutsCouple(),
                                                               aStep->GetStepLength(),
                                                               edep,
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
  }

}
