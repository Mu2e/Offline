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


//G4 includes
#include "G4AutoLock.hh"


namespace mu2e {

    
    EventStash::EventStash(const fhicl::ParameterSet& pset)
    :
    initialStashSize(0),
    simParticlePrinter_(pset.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet()))

    {
        std::cout << "AT new EventStash c'tor" << std::endl;
    }
    
    
void EventStash::initializeStash(int stash_size)
    {
        initialStashSize = stash_size;
        
        _myVectorOfPerEventData.reserve(initialStashSize);
        
        //this initializes the vector of StashPerEventData to have "stash_size" elements
        for (int i = 0; i < initialStashSize ; i++) {
            //std::cout << "calling _myVectorOfPerEventData.emplace_back for the " << i+1 <<  "th time" << std::endl;
            _myVectorOfPerEventData.emplace_back(0);
        }
        
    }
 
    
    /* *********************************************************** */
    /* ******************* put data into Stash ******************* */
    /* *********************************************************** */
    
    
//NOTE: these functions are protected in MT mode by the mutex
//which is locked inside EventAction::EndofEventAction
//where these functions are called
void EventStash::printInfo(int position)
    
    {
        simParticlePrinter_.print(std::cout, *_myVectorOfPerEventData.at(position).sims);
    }
    
    
void EventStash::insertData(int position_to_insert, int instance_num,
                                     std::unique_ptr<StatusG4> g4_status,
                                     std::unique_ptr<SimParticleCollection> sim_collection)
    {
        //_myVectorOfPerEventData.at(position_to_insert) = StashPerEventData(integer_to_insert);
        _myVectorOfPerEventData.at(position_to_insert).instanceNumber = instance_num;
        
        if (position_to_insert%10 == 0) {
        std::cout << "calling EventStash::insertData() for the --" << position_to_insert << "-- element of the stash"  << std::endl;
        }
        
        _myVectorOfPerEventData.at(position_to_insert).sims = std::move(sim_collection);
        _myVectorOfPerEventData.at(position_to_insert).status = std::move(g4_status);
        
        //simParticlePrinter_.print(std::cout, *_myVectorOfPerEventData.at(position_to_insert).sims);
        
        //printInfo(position_to_insert);
        
    }
    
    
void EventStash::insertTVDHits(int position_to_insert,
                               std::unique_ptr<StepPointMCCollection> tvd_hits,
                               std::string tvd_output_name)
    {
        _myVectorOfPerEventData.at(position_to_insert).tvd.first = tvd_output_name;
        _myVectorOfPerEventData.at(position_to_insert).tvd.second = std::move(tvd_hits);
    }
    
    
void EventStash::insertMCTrajectoryCollection(int position_to_insert,
                                                       std::unique_ptr<MCTrajectoryCollection> mc_trajectories)
    {
        _myVectorOfPerEventData.at(position_to_insert).trajectories = std::move(mc_trajectories);
    }

    
void EventStash::insertSimsRemapping(int position_to_insert,
                                              std::unique_ptr<SimParticleRemapping> sims_remap)
    {
        _myVectorOfPerEventData.at(position_to_insert).simRemapping = std::move(sims_remap);
    }
    
    
void EventStash::insertExtMonFNALSimHits(int position_to_insert,
                                                  std::unique_ptr<ExtMonFNALSimHitCollection> ext_mon_fnal_hits)
    {
        _myVectorOfPerEventData.at(position_to_insert).extMonFNALHits = std::move(ext_mon_fnal_hits);
    }
    
    

void EventStash::insertSDStepPointMC(int position_to_insert,
                                       std::unique_ptr<StepPointMCCollection> step_point_mc,
                                       std::string instance_name)
    {
        _myVectorOfPerEventData.at(position_to_insert).sensitiveDetectorSteps[instance_name] = std::move(step_point_mc);
    }

    
    
void EventStash::insertCutsStepPointMC(int position_to_insert,
                                         std::unique_ptr<StepPointMCCollection> step_point_mc,
                                         std::string instance_name)
    {
        _myVectorOfPerEventData.at(position_to_insert).cutsSteps[instance_name] = std::move(step_point_mc);
    }


    
    /* *********************************************************** */
    /* ******************* get data from Stash ******************* */
    /* *********************************************************** */
 
    
int EventStash::getInstanceNumber(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return _myVectorOfPerEventData.at(position_to_get).instanceNumber;
    }
    
    
std::unique_ptr<SimParticleCollection> EventStash::getSimPartCollection(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return std::move(_myVectorOfPerEventData.at(position_to_get).sims);
    }

    
std::unique_ptr<StatusG4> EventStash::getG4Status(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return std::move(_myVectorOfPerEventData.at(position_to_get).status);
    }
    
    
    //COMBINE THESE TWO?
//*****************************
std::string EventStash::getTVDName(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return _myVectorOfPerEventData.at(position_to_get).tvd.first;
    }

std::unique_ptr<StepPointMCCollection> EventStash::getTVDHits(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return std::move(_myVectorOfPerEventData.at(position_to_get).tvd.second);
    }
//*****************************
    
    
std::unique_ptr<MCTrajectoryCollection> EventStash::getMCTrajCollection(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return std::move(_myVectorOfPerEventData.at(position_to_get).trajectories);
    }
    
    
std::unique_ptr<SimParticleRemapping> EventStash::getSimParticleRemap(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return std::move(_myVectorOfPerEventData.at(position_to_get).simRemapping);
    }
    
    
std::unique_ptr<ExtMonFNALSimHitCollection> EventStash::getExtMonFNALSimHitCollection(int position_to_get)
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        return std::move(_myVectorOfPerEventData.at(position_to_get).extMonFNALHits);
    }
  
 
/*
std::unique_ptr<StepPointMCCollection> EventStash::getSDStepPointMC(int position_to_get,
                                       std::string instance_name)
    {
        return std::move(_myVectorOfPerEventData.at(position_to_get).sensitiveDetectorSteps[instance_name]);
    }
*/


    
// Put all of the SD data into the event.
// To be called at the end of each event.
void EventStash::putSensitiveDetectorData(int position_to_put, art::Event& event)
    {
    
        std::map< std::string, std::unique_ptr<StepPointMCCollection> > steps_map =
            std::move(_myVectorOfPerEventData.at(position_to_put).sensitiveDetectorSteps);

        for (std::map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = steps_map.begin();
             i != steps_map.end(); i++) {
            
            
            //***NEW STUFF
            for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){
                
                StepPointMC& step = *j;
                
                if ( step.simParticle().isNonnull() ){
                    step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                               step.simParticle().key(),
                                                               event.productGetter( step.simParticle().id() ) );
                }
            }//***END NEW STUFF
            
            event.put(std::move(i->second), i->first );
        }//for (std::map...
    }
    
    
// Put all of the cuts data into the event.
// To be called at the end of each event.
void EventStash::putCutsData(int position_to_put, art::Event& event)
    {
        //this is idenitcal in format to the put for the SD data
        //does this work correctly for the different kinds of cuts?

        //G4AutoLock stash_lock(&EventStashMutex);
        std::map< std::string, std::unique_ptr<StepPointMCCollection> > steps_map =
            std::move(_myVectorOfPerEventData.at(position_to_put).cutsSteps);
        
        for (std::map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = steps_map.begin();
             i != steps_map.end(); i++) {
            
            event.put(std::move(i->second), i->first );
        }
        
    }
    
/*
void EventStash::popOffLastElement()
    {
        //G4AutoLock stash_lock(&EventStashMutex);
        _myVectorOfPerEventData.pop_back();
        
    }
    

void EventStash::eraseStoredData()
    {
        _myVectorOfPerEventData.erase(_myVectorOfPerEventData.begin());
        
    }
*/
void EventStash::clearStash()
    {
        
        _myVectorOfPerEventData.clear();
        
    }

    
} // end namespace mu2e
