//
// PerEventObjectsManager.cc ********************
//
// Author: Lisa Goodenough
// Date: 2017/05/26
//
//

//G4 includes
#include "G4Threading.hh"

//Mu2e includes
#include "Mu2eG4/inc/PerEventObjectsManager.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"

// framework includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

//C++ includes
#include <iostream>
#include <sstream>


namespace mu2e {
    
    PerEventObjectsManager::PerEventObjectsManager(const fhicl::ParameterSet& pset,
                                                   int number_of_threads)
    :
    multiStagePars_(pset.get<fhicl::ParameterSet>("MultiStageParameters")),
    eventInstanceNumber(-1)
    {}
    
    PerEventObjectsManager::~PerEventObjectsManager()
    {}
    

void PerEventObjectsManager::createSimParticleHelpers(const art::ProductID& sim_part_ID,
                                                      const art::Event* evt,
                                                      const art::Handle<GenParticleCollection>* genparthandle,
                                                      const art::EDProductGetter* sim_product_getter)
    {
        //_spHelper = new SimParticleHelper(multiStagePars_.simParticleNumberOffset(), sim_part_ID, evt);
        //_parentHelper = new SimParticlePrimaryHelper(evt, sim_part_ID, *genparthandle);
        
        //_spHelper = std::make_unique<SimParticleHelper> ( SimParticleHelper(multiStagePars_.simParticleNumberOffset(), sim_part_ID, evt) );
        
        _spHelper = std::make_shared<SimParticleHelper> ( multiStagePars_.simParticleNumberOffset(), sim_part_ID, evt, sim_product_getter );
        _parentHelper = std::make_shared<SimParticlePrimaryHelper> ( evt, sim_part_ID, *genparthandle, sim_product_getter );
    }
    

void PerEventObjectsManager::storeEventInstanceNumber(int instance_number)
    {
        eventInstanceNumber = instance_number;
    }


} // end namespace mu2e
