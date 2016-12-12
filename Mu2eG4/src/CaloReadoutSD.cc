//
// Define a sensitive detector for calorimetric readout
//
// $Id: CaloReadoutSD.cc,v 1.22 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"



namespace mu2e {

    CaloReadoutSD::CaloReadoutSD(G4String name, SimpleConfig const & config ):
      Mu2eSensitiveDetector(name,config),_nro(0)
    {
	GeomHandle<Calorimeter> cg;
	_nro  = cg->caloGeomInfo().nROPerCrystal();
    }



    G4bool CaloReadoutSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
    {

	 if( aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0 ) return false;
         if( aStep->GetTotalEnergyDeposit() < 1e-6 ) return false;

	 _currentSize += 1;
	 if( _sizeLimit>0 && _currentSize>_sizeLimit && (_currentSize - _sizeLimit)==1) 
	 {
             mf::LogWarning("G4") << "Maximum number of particles reached in " << SensitiveDetectorName<< ": "<< _currentSize <<std::endl;
             return false;
	 }


	 ProcessCode endCode(_processInfo->findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

	 const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();
	 int idro = touchableHandle->GetCopyNumber(0) + touchableHandle->GetCopyNumber(2)*_nro;  //the idro is always Number(0) + _nro*number(X), make sure X is right


	 //for diagnosis purposes only when playing with the geometry, uncomment next line
	 //for (int i=0;i<8;++i) std::cout<<"RO Transform level "<<i<<"   "<<touchableHandle->GetCopyNumber(i)<<"     "
	 //                               <<touchableHandle->GetHistory()->GetTransform(i).TransformPoint(aStep->GetPreStepPoint()->GetPosition())
	 //				  <<"   "<<idro<<"   "<<touchableHandle->GetSolid(i)->GetName()<<"   "<<touchableHandle->GetVolume(i)->GetName()<<std::endl;

	 _collection->push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                                	    idro,
                                	    aStep->GetTotalEnergyDeposit(),
                                	    aStep->GetNonIonizingEnergyDeposit(),
                                	    aStep->GetPreStepPoint()->GetGlobalTime(),
                                	    aStep->GetPreStepPoint()->GetProperTime(),
                                	    aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                                	    aStep->GetPreStepPoint()->GetMomentum(),
                                	    aStep->GetStepLength(),
                                	    endCode
                                	    ) );

	 return true;
    }


} 
