#ifndef Mu2eG4_GenEventBroker_hh
#define Mu2eG4_GenEventBroker_hh
//
// GenEventBroker.hh *********************
// ***************************************
//
// Author: Lisa Goodenough
// Date: 2017/05/23
//
//


//Mu2e includes
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollections.hh"


//art includes
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"


//G4 includes
#include "G4Threading.hh"


namespace art { class ProductID; }
namespace art { class Event; }
namespace art { class EDProductGetter; }

namespace mu2e {
    
    typedef std::vector< art::ValidHandle<StepPointMCCollection> > HitHandles;
    class EventStash;
    
class GenEventBroker
{
  public:
    GenEventBroker(const bool using_MT);
        
    // called in Mu2eG4::produce, takes in the GenParticleCollection for the event
    void loadEvent(HitHandles const&, art::ProductID const&,
                   art::Event*, art::InputTag, EventStash*,
                   art::EDProductGetter const*);

    void setEventPtrToZero();
    
    //we need this handle to give to the SimParticlePrimaryHelper in PrimaryGenAction
    art::Handle<GenParticleCollection> const& getGenParticleHandle() const;
    
    inline HitHandles const* getHitHandles() const { return _genInputHits; }
    inline art::ProductID getproductID() const { return _simParticleID; }
    inline art::Event* getartEvent() const { return _artEvent; }
    inline EventStash* getEventStash() const { return _stashForEventData; }
    inline int getStashSize() const { return eventStashSize; }//this can only be called after loadEvent is called in produce
    inline art::EDProductGetter const* getSimProductGetter() const { return _simProductGetter; }
    
  private:
    
    const HitHandles* _genInputHits;    
    art::ProductID _simParticleID;
    art::Event* _artEvent;
    EventStash* _stashForEventData;
    const art::EDProductGetter* _simProductGetter;
    
    art::Handle<GenParticleCollection> gensHandle;
    art::Handle<GenParticleCollections> genCollectionsHandle;
    
    art::InputTag _generatorModuleLabel;
    
    int eventStashSize;
    int eventInstanceToGet;
    
    GenParticleCollections genParticleCollectionStash;
    G4Mutex GenParticleCollectionsMutex = G4MUTEX_INITIALIZER;
    
    const bool usingG4MT;
    

  public:
    
    struct GenParticleCollectionInstance{
        
        explicit GenParticleCollectionInstance(int instance_int):
        instanceNumber(instance_int)//atomic int counter needed here?
        {}
        
        int instanceNumber;
        const GenParticleCollection* genCollection;
        
        };

    GenParticleCollectionInstance getNextGenPartCollectionInstance();
    
    
};

}  // end namespace mu2e
#endif /* Mu2eG4_GenEventBroker_hh */


