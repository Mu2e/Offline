//
// Define a sensitive detector for calorimetric readout
//
// $Id: CaloReadoutSD.cc,v 1.7 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e incldues
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  G4ThreeVector CaloReadoutSD::_mu2eOrigin;

  CaloReadoutSD::CaloReadoutSD(G4String name, const SimpleConfig& config)
    : G4VSensitiveDetector(name),
      _collection(0),
      _nro(0),
      _minE(0.0),
      _debugList(0),
      _sizeLimit(config.get<int>("g4.stepsSizeLimit",0)),
      _currentSize(0)
  {

    // Get list of events for which to make debug printout.
    string key("g4.calorimeterSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

  }


  CaloReadoutSD::~CaloReadoutSD(){ }

  void CaloReadoutSD::Initialize(G4HCofThisEvent* HCE){

    _currentSize=0;

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
	mf::LogWarning("G4") << "Maximum number of particles reached in CaloCrystalSD: "
			      << _currentSize << endl;
      }
      return false;
    }

    // Get readout ID
    int idro = touchableHandle->GetCopyNumber(0);
    // in the previous version of calorimeter geometry the RO id
    // had to be calculated this way:
    // int idro = touchableHandle->GetCopyNumber(0) + touchableHandle->GetCopyNumber(1)*_nro;

    // The points coordinates are saved in the mu2e world

    // we add the hit to the framework collection
    _collection->
      push_back(StepPointMC(aStep->GetTrack()->GetTrackID(),
                            idro,
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetStepLength()
                            ));

    return true;

  }


  void CaloReadoutSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize
			    << " calorimeter RO hits were generated in the event."
			    << endl
			    << "Only " << _sizeLimit << " are saved in output collection."
			    << endl;
    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " RO hits in the calorimeter: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }

  }


  void CaloReadoutSD::beforeG4Event(StepPointMCCollection& outputHits) {
    _collection = &outputHits;
    return;
  } // end of beforeG4Event


} //namespace mu2e
