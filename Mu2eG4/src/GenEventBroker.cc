//
// GenEventBroker.cc ********************
//
// Author: Lisa Goodenough
// Date: 2017/05/23
//
// Notes:
// 1) In MT mode both the GenParticleCollection and the GenParticleCollections
//    data products must be present.
//    In sequential mode GenParticleCollections should never be present.
//    In sequential mode GenParticleCollection may be absent when multistage mode is enabled.
//

//Mu2e includes
#include "Mu2eG4/inc/GenEventBroker.hh"
#include "Mu2eG4/inc/EventStash.hh"


// framework includes
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

//C++ includes
#include <iostream>

// G4 includes
#include "G4AutoLock.hh"


namespace mu2e {

    GenEventBroker::GenEventBroker(const bool using_MT)
    :
    _genInputHits(0),
    _simParticleID(),
    _artEvent(0),
    _stashForEventData(0),
    _simProductGetter(0),
    eventStashSize(0),
    eventInstanceToGet(0),
    usingG4MT(using_MT)

    {}


void GenEventBroker::loadEvent(HitHandles const& gen_input_hits,
                               art::ProductID const& sim_part_ID,
                               art::Event* art_event,
                               art::InputTag genmoduleLabel,
                               EventStash* stash_for_event_data,
                               art::EDProductGetter const* sim_prod_getter
                               )
    {

        _genInputHits = &gen_input_hits;
        _simParticleID = sim_part_ID;
        _artEvent = art_event;
        _generatorModuleLabel = genmoduleLabel;
        _stashForEventData = stash_for_event_data;
        _simProductGetter = sim_prod_getter;

        //need to reset this value every time loadEvent() is called
        //which corresponds to an empty stash
        eventInstanceToGet = 0;

        //no lock on art::event needed here because GEB exists only in Master Thread
        if(!(_generatorModuleLabel == art::InputTag())) {
            _artEvent->getByLabel(_generatorModuleLabel, gensHandle);
        }

        if(!(_generatorModuleLabel == art::InputTag())) {
            _artEvent->getByLabel(_generatorModuleLabel, genCollectionsHandle);
        }

        // In MT mode both must be present.  See Note 1.
        if ( usingG4MT && !( genCollectionsHandle.isValid() && gensHandle.isValid() ) )
        {
            throw cet::exception("CONFIG")
            << "Error in GenEventBroker::loadEvent. You are trying to run in MT mode and there is no GenParticleCollection(s)!\n";
        }
        
        if (usingG4MT) {//MT mode
            eventStashSize = genCollectionsHandle.product()->size();
        }
        else//sequential mode
        {
            eventStashSize = 1;
        }
     }


void GenEventBroker::setEventPtrToZero()
    {
        _artEvent = 0;
    }


art::Handle<GenParticleCollection> const& GenEventBroker::getGenParticleHandle() const
    {
        //no lock needed here
        return gensHandle;
    }


GenEventBroker::GenParticleCollectionInstance GenEventBroker::getNextGenPartCollectionInstance()
    {
        //need a lock here to protect from multiple threads accessing at the same time
        G4AutoLock GPCs_lock(&GenParticleCollectionsMutex);

        GenParticleCollectionInstance GPCInstance(eventInstanceToGet);
        GPCInstance.genCollection = &(genCollectionsHandle.product()->at(eventInstanceToGet));

        eventInstanceToGet++;

        return GPCInstance;

    }


} // end namespace mu2e
