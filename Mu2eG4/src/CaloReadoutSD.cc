//
// Define a sensitive detector for calorimetric readout
//
// $Id: CaloReadoutSD.cc,v 1.18 2012/11/20 23:02:09 echenard Exp $
// $Author: echenard $
// $Date: 2012/11/20 23:02:09 $
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  CaloReadoutSD::CaloReadoutSD(G4String name, SimpleConfig const & config ):
    Mu2eSensitiveDetector(name,config),
    _nro(0),
    _minE(0.0)
  {

    GeomHandle<Calorimeter> cg;
    _nro  = cg->nROPerCrystal();
    _minE = cg->getElectronEmin();

  }

  G4bool CaloReadoutSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    // Only handle charged events with kinetic energy > 0.1 MeV

    if( aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0 ) return false;
    if( aStep->GetTrack()->GetKineticEnergy() < _minE ) return false;

    // Check that number of steps did not exceed the limit

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in " 
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    // Get readout ID
    int idro = touchableHandle->GetCopyNumber(0);
    // in the previous version of calorimeter geometry the RO id
    // had to be calculated this way:
    // int idro = touchableHandle->GetCopyNumber(0) + touchableHandle->GetCopyNumber(1)*_nro;

    // Which process caused this step to end?
    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    ProcessCode endCode(_processInfo->findAndCount(pname));

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.

      push_back(StepPointMC(art::Ptr<SimParticle>
                            ( *_simID,
                              aStep->GetTrack()->GetTrackID(),
                              _event->productGetter(*_simID) ),
                            idro,
                            aStep->GetTotalEnergyDeposit(),
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


} //namespace mu2e
