#ifndef Mu2eG4_EventAction_hh
#define Mu2eG4_EventAction_hh
//
// G4 begin and end of event actions for Mu2e.
//
// Author: Lisa Goodeough
// Date: 2017/05/04
//

//G4 includes
#include "Geant4/G4UserEventAction.hh"
#include "Geant4/G4Threading.hh"

//Mu2e includes
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "Offline/MCDataProducts/inc/StepInstanceName.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"

//art includes
#include "art/Framework/Principal/Handle.h"

//C++ includes
#include <vector>
#include <memory>

class G4Event;
class G4Timer;

namespace art { class Event; }

namespace mu2e {

  class Mu2eG4TrackingAction;
  class Mu2eG4SteppingAction;
  class SensitiveDetectorHelper;
  class PhysicsProcessInfo;
  class IMu2eG4Cut;

  typedef std::map<art::Ptr<SimParticle>, art::Ptr<SimParticle> >  SimParticleRemapping;


  class Mu2eG4EventAction : public G4UserEventAction
  {
  public:

    Mu2eG4EventAction(const Mu2eG4Config::Top&,
                      Mu2eG4TrackingAction*,
                      Mu2eG4SteppingAction*,
                      SensitiveDetectorHelper*,
                      Mu2eG4PerThreadStorage* pts,
                      PhysicsProcessInfo*,
                      const CLHEP::Hep3Vector&
                      );

    virtual ~Mu2eG4EventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);


  private:

    Mu2eG4PerThreadStorage* perThreadObjects_;

    //these are set using fhicl pset
    SimParticleCollectionPrinter simParticlePrinter_;

    Mu2eG4TrackingAction* _trackingAction;
    Mu2eG4SteppingAction* _steppingAction;

    SensitiveDetectorHelper* _sensitiveDetectorHelper = nullptr;

    const CLHEP::Hep3Vector& _originInWorld;

    // local Mu2e per Geant4 event timer
    std::unique_ptr<G4Timer> _timer;

    PhysicsProcessInfo *_processInfo = nullptr;

    bool _g4InternalFiltering;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_EventAction_hh */
