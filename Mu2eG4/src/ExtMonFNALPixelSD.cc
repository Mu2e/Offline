// $Id: ExtMonFNALPixelSD.cc,v 1.1 2012/08/23 23:36:14 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:36:14 $
//
// Original author Andrei Gaponenko

#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"

#include <sstream>

#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4ios.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

namespace mu2e {

  //================================================================
  ExtMonFNALPixelSD::ExtMonFNALPixelSD(const SimpleConfig& config)
    : G4VSensitiveDetector(SensitiveDetectorName::ExtMonFNAL())
    , _sizeLimit(std::max(0, config.getInt("g4.stepsSizeLimit",0)))
    , _simID()
    , _event()
  {
    verboseLevel = config.getInt("extMonFNAL.sd.verbosityLevel");
  }

  //================================================================
  void ExtMonFNALPixelSD::beforeG4Event(ExtMonFNALSimHitCollection *outputHits,
                                        const art::ProductID& simID,
                                        const art::Event& event )
  {
    _collection  = outputHits;
    _simID       = &simID;
    _event       = &event;
  }


  //================================================================
  G4bool ExtMonFNALPixelSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    bool retval = false;

    G4double totalEDep = aStep->GetTotalEnergyDeposit();

    if(totalEDep > 0.) {
      if(!_sizeLimit||(_collection->size() < _sizeLimit)) {

        retval = true;

        // These are in G4 world coordinates
        const G4ThreeVector globalStartPos = aStep->GetPreStepPoint()->GetPosition();
        const G4ThreeVector globalEndPos = aStep->GetPostStepPoint()->GetPosition();

        // Transform the global coordinates into the sensor frame.
        G4TouchableHistory*  tch = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
        const G4AffineTransform transformation = tch->GetHistory()->GetTopTransform();
        G4ThreeVector localStartPos = transformation.TransformPoint(globalStartPos);
        G4ThreeVector localEndPos = transformation.TransformPoint(globalEndPos);

        art::Ptr<SimParticle> particle( *_simID, aStep->GetTrack()->GetTrackID(), _event->productGetter(*_simID) );
        ExtMonFNALSensorId sid(aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());

        // Add the hit to the framework collection.
        _collection->
          push_back(ExtMonFNALSimHit(sid,
                                     particle,

                                     totalEDep,
                                     aStep->GetNonIonizingEnergyDeposit(),

                                     localStartPos,
                                     aStep->GetPreStepPoint()->GetGlobalTime(),

                                     localEndPos,
                                     aStep->GetPostStepPoint()->GetGlobalTime()
                                ));

      }
    }

    return retval;

  }

  //================================================================
  void ExtMonFNALPixelSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && (_collection->size() >=_sizeLimit) ){
      std::ostringstream os;
      os<< "Total of " << _collection->size()<< " ExtMonFNAL hits were generated in the event.\n"
        << "Only " << _sizeLimit << " are saved in output collection.\n";

      mf::LogWarning("G4")<<os.str();
      G4cout<<os.str()<<G4endl;
    }

    if (verboseLevel>0) {
      G4cout << "\n-------->Hits Collection: in this event there are " << _collection->size()
             << " ExtMonFNAL hits" << G4endl;
    }
  }

  //================================================================

} //namespace mu2e
