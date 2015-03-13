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
#include "cetlib/exception.h"

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

	// Calculate energy deposition in the crystal
	G4double edep = aStep->GetTotalEnergyDeposit();
	if( edep<=0 ) return false;

	_currentSize += 1;

	if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
	  if( (_currentSize - _sizeLimit)==1 ) {
            mf::LogWarning("G4") << "Maximum number of particles reached in " 
                        	 << SensitiveDetectorName
                        	 << ": "
                        	 << _currentSize << std::endl;
	  }
	  return false;
	}

	const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

	// Get crystal ID
	G4int copyNo = touchableHandle->GetCopyNumber(0);

	ProcessCode endCode(_processInfo->findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

	// Originally the hit position was saved in local crystal frame.
	// Not it is saved in Mu2e frame, hence the following code is
	// commented out.
	// Calculate enerdy deposition position along the crystal
	 G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
	 G4ThreeVector posWorld = aStep->GetPreStepPoint()->GetPosition();
	 G4ThreeVector posLocal = toLocal.TransformPoint(posWorld);


         //for diagnosis only, uncomment next two lines
	 //if (_currentSize < 2)
	 //  for (int i=0;i<8;++i) std::cout<<"Transform level "<<i<<"   "<<touchableHandle->GetCopyNumber(i)<<"     "<<touchableHandle->GetHistory()->GetTransform(i).TransformPoint(posWorld)<<std::endl;
         

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
