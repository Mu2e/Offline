#ifndef Mu2eG4_EventAction_hh
#define Mu2eG4_EventAction_hh
//
// G4 begin and end of event actions for Mu2e.
//
// Author: Lisa Goodeough
// Date: 2017/05/04
//

//G4 includes
#include "G4UserEventAction.hh"
#include "G4Threading.hh"

//Mu2e includes
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

//art includes
#include "art/Framework/Principal/Handle.h"

//C++ includes
#include <vector>
#include <memory>

    class G4Event;
    class G4Timer;

namespace fhicl { class ParameterSet; }
namespace art { class Event; }

namespace mu2e {
    
    class TrackingAction;
    class Mu2eG4SteppingAction;
    class SensitiveDetectorHelper;
    class SimParticleHelper;
    class SimParticlePrimaryHelper;
    class PhysicsProcessInfo;
    class IMu2eG4Cut;
    class GenEventBroker;
    class PerEventObjectsManager;
    class EventStash;
    
    typedef std::map<art::Ptr<SimParticle>, art::Ptr<SimParticle> >  SimParticleRemapping;
    

class Mu2eG4EventAction : public G4UserEventAction
{
  public:
    
    Mu2eG4EventAction(const fhicl::ParameterSet&,
                      TrackingAction*,
                      Mu2eG4SteppingAction*,
                      SensitiveDetectorHelper*,
                      IMu2eG4Cut&,
                      IMu2eG4Cut&,
                      IMu2eG4Cut&,
                      GenEventBroker*,
                      PerEventObjectsManager*,
                      PhysicsProcessInfo*,
                      const CLHEP::Hep3Vector&
                      );
    
    virtual ~Mu2eG4EventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    
  private:
    
    //used to set the art::Event from the GenEventBroker
    void setEventData();
    
    //these are set using fhicl pset
    Mu2eG4TrajectoryControl trajectoryControl_;
    SimParticleCollectionPrinter simParticlePrinter_;
    std::vector<double> timeVDtimes_;
    Mu2eG4MultiStageParameters multiStagePars_;

    
    TrackingAction* _trackingAction;
    Mu2eG4SteppingAction* _steppingAction;
    
    SensitiveDetectorHelper* _sensitiveDetectorHelper;
    
    IMu2eG4Cut* _stackingCuts;
    IMu2eG4Cut* _steppingCuts;
    IMu2eG4Cut* _commonCuts;
    const CLHEP::Hep3Vector& _originInWorld;
    
    GenEventBroker *_genEventBroker;
    PerEventObjectsManager* perEvtObjManager;
    
    const StepInstanceName _tvdOutputName;

    // local Mu2e per Geant4 event timer
    std::unique_ptr<G4Timer> _timer;
    
    // Create empty data products.
    std::unique_ptr<SimParticleCollection> simParticles;
    std::unique_ptr<StepPointMCCollection> tvdHits;
    std::unique_ptr<MCTrajectoryCollection> mcTrajectories;
    std::unique_ptr<SimParticleRemapping> simsRemap;
    std::unique_ptr<ExtMonFNALSimHitCollection> extMonFNALHits;
    
    SimParticleHelper *_spHelper;
    SimParticlePrimaryHelper *_parentHelper;
    PhysicsProcessInfo *_processInfo;
    
    //these are set in setEventData
    art::Event *_artEvent;
    EventStash *_stashForEventData;
    int eventNumberInProcess;
    
    bool _g4InternalFiltering;

    // products for the g4study
    //std::unique_ptr<StepPointMCCollection> steppingPoints;    
};

}  // end namespace mu2e
#endif /* Mu2eG4_EventAction_hh */


