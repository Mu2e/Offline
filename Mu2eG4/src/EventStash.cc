//
// EventStash.cc ********************
//
// Author: Lisa Goodenough
// Date: 2017/10/11
//
//

//Mu2e includes
#include "Mu2eG4/inc/EventStash.hh"

// framework includes
#include "art/Framework/Principal/Event.h"

//C++ includes
#include <iostream>


namespace mu2e {

    
    EventStash::EventStash(const fhicl::ParameterSet& pset)
    :
    initialStashSize(0),
    simParticlePrinter_(pset.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet()))
    {}
    
    
void EventStash::initializeStash(int stash_size)
    {
        initialStashSize = stash_size;
        
        _myVectorOfPerEventData.reserve(initialStashSize);
        
        //this initializes the vector of StashPerEventData to have "stash_size" elements
        for (int i = 0; i < initialStashSize ; i++) {
            _myVectorOfPerEventData.emplace_back();
        }
        
    }
 
    
    /* *********************************************************** */
    /* ******************* put data into Stash ******************* */
    /* *********************************************************** */
    
void EventStash::printInfo(int position)
    
    {
        simParticlePrinter_.print(std::cout, *_myVectorOfPerEventData[position].sims);
    }
    
    
void EventStash::insertData(int position_to_insert,
                            std::unique_ptr<StatusG4> g4_status,
                            std::unique_ptr<SimParticleCollection> sim_collection)
    {
        _myVectorOfPerEventData[position_to_insert].sims = std::move(sim_collection);
        _myVectorOfPerEventData[position_to_insert].status = std::move(g4_status);
        //simParticlePrinter_.print(std::cout, *_myVectorOfPerEventData[position_to_insert].sims);
        //printInfo(position_to_insert);
    }
    
    
void EventStash::insertTVDHits(int position_to_insert,
                               std::unique_ptr<StepPointMCCollection> tvd_hits,
                               std::string tvd_output_name)
    {
        _myVectorOfPerEventData[position_to_insert].tvd.first = tvd_output_name;
        _myVectorOfPerEventData[position_to_insert].tvd.second = std::move(tvd_hits);
    }
    
    
void EventStash::insertMCTrajectoryCollection(int position_to_insert,
                                                       std::unique_ptr<MCTrajectoryCollection> mc_trajectories)
    {
        _myVectorOfPerEventData[position_to_insert].trajectories = std::move(mc_trajectories);
    }

    
void EventStash::insertSimsRemapping(int position_to_insert,
                                              std::unique_ptr<SimParticleRemapping> sims_remap)
    {
        _myVectorOfPerEventData[position_to_insert].simRemapping = std::move(sims_remap);
    }
    
    
void EventStash::insertExtMonFNALSimHits(int position_to_insert,
                                                  std::unique_ptr<ExtMonFNALSimHitCollection> ext_mon_fnal_hits)
    {
        _myVectorOfPerEventData[position_to_insert].extMonFNALHits = std::move(ext_mon_fnal_hits);
    }
    
    

void EventStash::insertSDStepPointMC(int position_to_insert,
                                       std::unique_ptr<StepPointMCCollection> step_point_mc,
                                       std::string instance_name)
    {
        _myVectorOfPerEventData[position_to_insert].sensitiveDetectorSteps[instance_name] = std::move(step_point_mc);
    }

    
    
void EventStash::insertCutsStepPointMC(int position_to_insert,
                                         std::unique_ptr<StepPointMCCollection> step_point_mc,
                                         std::string instance_name)
    {
        _myVectorOfPerEventData[position_to_insert].cutsSteps[instance_name] = std::move(step_point_mc);
    }


    
/* *********************************************************** */
/* ******************* get data from Stash ******************* */
/* *********************************************************** */
 

std::unique_ptr<SimParticleCollection> EventStash::getSimPartCollection(int position_to_get)
    {
        return std::move(_myVectorOfPerEventData[position_to_get].sims);
    }

    
std::unique_ptr<StatusG4> EventStash::getG4Status(int position_to_get)
    {
        return std::move(_myVectorOfPerEventData[position_to_get].status);
    }
    
    
std::string EventStash::getTVDName(int position_to_get)
    {
        return _myVectorOfPerEventData[position_to_get].tvd.first;
    }

std::unique_ptr<StepPointMCCollection> EventStash::getTVDHits(int position_to_get)
    {
        return std::move(_myVectorOfPerEventData[position_to_get].tvd.second);
    }
    
    
std::unique_ptr<MCTrajectoryCollection> EventStash::getMCTrajCollection(int position_to_get)
    {
        return std::move(_myVectorOfPerEventData[position_to_get].trajectories);
    }
    
    
std::unique_ptr<SimParticleRemapping> EventStash::getSimParticleRemap(int position_to_get)
    {
        return std::move(_myVectorOfPerEventData[position_to_get].simRemapping);
    }
    
    
std::unique_ptr<ExtMonFNALSimHitCollection> EventStash::getExtMonFNALSimHitCollection(int position_to_get)
    {
        return std::move(_myVectorOfPerEventData[position_to_get].extMonFNALHits);
    }
  
    
// Put all of the SD data into the event.
// To be called at the end of each event.
void EventStash::putSensitiveDetectorData(int position_to_put, art::Event& event, art::EDProductGetter const* sim_product_getter)
    {
    
        std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > steps_map =
            std::move(_myVectorOfPerEventData[position_to_put].sensitiveDetectorSteps);

        for (std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = steps_map.begin();
             i != steps_map.end(); ++i) {
            
            for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){
                
                StepPointMC& step = *j;
                
                if ( step.simParticle().isNonnull() ){
                    step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                               step.simParticle().key(),
                                                               sim_product_getter );
                }//if
            }//for StepPointMCCollection::iterator
            
            event.put(std::move(i->second), i->first );
        }//for (std::unordered_map...
    }
    
    
// Put all of the cuts data into the event.
// To be called at the end of each event.
void EventStash::putCutsData(int position_to_put, art::Event& event, art::EDProductGetter const* sim_product_getter)
    {

        std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > cuts_map =
            std::move(_myVectorOfPerEventData[position_to_put].cutsSteps);
        
        for (std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = cuts_map.begin();
             i != cuts_map.end(); ++i) {
            
            for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){
                StepPointMC& step = *j;
                
                if ( step.simParticle().isNonnull() ){
                    step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                               step.simParticle().key(),
                                                               sim_product_getter );
                }//if
            }//for StepPointMCCollection::iterator
                        
            event.put(std::move(i->second), i->first );
        }
        
    }
    
void EventStash::clearStash()
    {
        _myVectorOfPerEventData.clear();
    }

    
} // end namespace mu2e
