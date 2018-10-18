//
// Define a sensitive detector for CaloCrystal Detectors
//
// $Id: CaloCrystalSD.cc,v 1.25 2013/08/28 05:58:17 gandr Exp $
// $Author: gandr $
// $Date: 2013/08/28 05:58:17 $
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/WorldG4.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"


namespace mu2e {

    CaloCrystalSD::CaloCrystalSD(G4String name, SimpleConfig const & config ): 
      Mu2eSensitiveDetector(name,config)
    { }

    G4bool CaloCrystalSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
    {
        G4double edep = aStep->GetTotalEnergyDeposit();
	if (edep < 1e-6) return false;


	_currentSize += 1;
	if( _sizeLimit>0 && _currentSize>_sizeLimit && (_currentSize - _sizeLimit)==1) 
	{
            mf::LogWarning("G4") << "Maximum number of particles reached in " << SensitiveDetectorName<< ": "<< _currentSize <<std::endl;
            return false;
	}


	ProcessCode endCode(_processInfo->findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

	const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

	G4int copyNo = touchableHandle->GetCopyNumber(1);  // Make sure to get the right copy number level here

	G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
	G4ThreeVector posWorld           = aStep->GetPreStepPoint()->GetPosition();
	G4ThreeVector posLocal           = toLocal.TransformPoint(posWorld);


	// diagnosis purposes only when playing with the geometry, uncomment next two line
	//for (int i=0;i<=touchableHandle->GetHistoryDepth();++i) std::cout<<"Cry Transform level "<<i<<"   "<<touchableHandle->GetCopyNumber(i)
        //<<"   "<<touchableHandle->GetHistory()->GetTransform(touchableHandle->GetHistoryDepth()-i).TransformPoint(posWorld)
	//<<"  "<<touchableHandle->GetSolid(i)->GetName()<<"   "<<touchableHandle->GetVolume(i)->GetName()<<std::endl;


	_collection->push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                        		   copyNo,
                        		   edep,
                        		   aStep->GetNonIonizingEnergyDeposit(),
                        		   aStep->GetPreStepPoint()->GetGlobalTime(),
                        		   aStep->GetPreStepPoint()->GetProperTime(),
                        		   aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                        		   aStep->GetPreStepPoint()->GetMomentum(),
                        		   aStep->GetStepLength(),
                        		   endCode
                        		   ));

	return true;

    }

} 
