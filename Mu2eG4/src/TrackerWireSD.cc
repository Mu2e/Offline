//
//
//  $Id: TrackerWireSD.cc,v 1.4 2013/08/28 05:58:17 gandr Exp $
//  $Author: gandr $
//  $Date: 2013/08/28 05:58:17 $
//
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/TrackerWireSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  G4ThreeVector TrackerWireSD::_mu2eDetCenter;

  TrackerWireSD::TrackerWireSD(G4String name, const SimpleConfig& config) :
                  Mu2eSensitiveDetector(name,config) { }


  TrackerWireSD::~TrackerWireSD(){ }

  G4bool TrackerWireSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){

          _currentSize += 1;

          if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
                  if( (_currentSize - _sizeLimit)==1 ) {
                          mf::LogWarning("G4") << "Maximum number of particles reached "
                                          << SensitiveDetectorName
                                          << ": "
                                          << _currentSize << endl;
                  }
                  return false;
          }

          G4double edep  = aStep->GetTotalEnergyDeposit();
          G4double nidep = aStep->GetNonIonizingEnergyDeposit();
          G4double step  = aStep->GetStepLength();
          //G4double idep  = edep-nidep;

          if ( _debugList.inList() )  cout<<"edep "<<edep<<" nidep "<<nidep<<" step "<<step<<endl;
          // I am not sure why we get these cases but we do.  Skip them.
          if ( (edep == 0. /*|| idep == 0.*/)/*&& step == 0.*/ ) {
                  if ( _debugList.inList() )  cout<<"Skipped"<<endl;
                  return false;
          }

          G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
          if ( preStepPoint->GetPhysicalVolume()->GetName().find("vol")!=string::npos ) { return false; }

          G4int motherCopyNo = preStepPoint->GetTouchableHandle()->GetReplicaNumber(1);

          if ( _debugList.inList() )  {
        	  cout<<"Step vol name "<<aStep->GetTrack()->GetVolume()->GetName()<<endl;
        	  cout<<"Step vol copyNo "<<" mother copyNo "<<motherCopyNo<<endl;
          }

          // Which process caused this step to end?
          ProcessCode endCode(_processInfo->
                              findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

          // Add the hit to the framework collection.
          // The point's coordinates are saved in the mu2e coordinate system.
          _collection->
            push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                                  motherCopyNo,
                                  aStep->GetTotalEnergyDeposit(),
                                  aStep->GetNonIonizingEnergyDeposit(),
                                  aStep->GetPreStepPoint()->GetGlobalTime(),
                                  aStep->GetPreStepPoint()->GetProperTime(),
                                  aStep->GetPreStepPoint()->GetPosition() - _mu2eDetCenter,
                                  aStep->GetPreStepPoint()->GetMomentum(),
                                  aStep->GetStepLength(),
                                  endCode
                                  ));

          return true;
   }

} //namespace mu2e
