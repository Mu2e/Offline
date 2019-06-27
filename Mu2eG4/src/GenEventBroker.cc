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


// framework includes
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

//C++ includes
#include <iostream>

namespace mu2e {

    GenEventBroker::GenEventBroker(const bool using_MT)
    :
    _genInputHits(0),
    _simParticleID(),
    _artEvent(0),
    _simProductGetter(0),
    usingG4MT(using_MT)
    {}


void GenEventBroker::loadEvent(HitHandles const& gen_input_hits,
                               art::ProductID const& sim_part_ID,
                               art::Event* art_event,
                               art::InputTag genmoduleLabel,
                               art::EDProductGetter const* sim_prod_getter
                               )
    {
        _genInputHits = &gen_input_hits;
        _simParticleID = sim_part_ID;
        _artEvent = art_event;
        _generatorModuleLabel = genmoduleLabel;
        _simProductGetter = sim_prod_getter;


        //no lock on art::event needed here because GEB exists only in Master Thread
        if(!(_generatorModuleLabel == art::InputTag())) {
            _artEvent->getByLabel(_generatorModuleLabel, gensHandle);
        }

        // In MT mode both must be present.  See Note 1.
        if ( usingG4MT && !gensHandle.isValid() )
        {
            throw cet::exception("CONFIG")
            << "Error in GenEventBroker::loadEvent. You are trying to run in MT mode and there is no GenParticleCollection!\n";
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



} // end namespace mu2e
