#ifndef Mu2eG4_EventStash_hh
#define Mu2eG4_EventStash_hh
//
// EventStash.hh *********************
// ***************************************
//
// Author: Lisa Goodenough
// Date: 2017/10/11
//
//


//Mu2e includes
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"


//G4 includes
#include "G4Threading.hh"


//C++ includes
#include <vector>
#include <memory>
#include <utility>


namespace art   { class Event; }
namespace art { class EDProductGetter; }

namespace mu2e {
        
class EventStash
{
  public:
    //EventStash(int stash_size);
    EventStash(const fhicl::ParameterSet&);
    
    inline int getStashSize() { return _myVectorOfPerEventData.size(); }

    void initializeStash(int stash_size);
    
    
    //methods to insert event data into the stash
    void insertData(int position_to_insert,
                    std::unique_ptr<StatusG4> g4_status,
                    std::unique_ptr<SimParticleCollection> sim_collection);
    
    void insertTVDHits(int position_to_insert,
                       std::unique_ptr<StepPointMCCollection> tvd_hits,
                       std::string tvd_output_name);
    
    void insertMCTrajectoryCollection(int position_to_insert,
                                      std::unique_ptr<MCTrajectoryCollection> mc_trajectories);
                                               
    void insertSimsRemapping(int position_to_insert,
                             std::unique_ptr<SimParticleRemapping> sims_remap);

    void insertExtMonFNALSimHits(int position_to_insert,
                                 std::unique_ptr<ExtMonFNALSimHitCollection> ext_mon_fnal_hits);
    
    void insertSDStepPointMC(int position_to_insert,
                             std::unique_ptr<StepPointMCCollection> step_point_mc,
                             std::string instance_name);
    
    void insertCutsStepPointMC(int position_to_insert,
                               std::unique_ptr<StepPointMCCollection> step_point_mc,
                               std::string instance_name);
    
    
    //a helper function for testing
    void printInfo(int position);


    //methods to get event data from the stash
    
    std::unique_ptr<SimParticleCollection> getSimPartCollection(int position_to_get);
    std::unique_ptr<StatusG4> getG4Status(int position_to_get);
    std::string getTVDName(int position_to_get);
    std::unique_ptr<StepPointMCCollection> getTVDHits(int position_to_get);
    std::unique_ptr<MCTrajectoryCollection> getMCTrajCollection(int position_to_get);
    std::unique_ptr<SimParticleRemapping> getSimParticleRemap(int position_to_get);
    std::unique_ptr<ExtMonFNALSimHitCollection> getExtMonFNALSimHitCollection(int position_to_get);
    
    void putSensitiveDetectorData(int position_to_put, art::Event& event, art::EDProductGetter const* sim_product_getter);
    
    void putCutsData(int position_to_put, art::Event& event, art::EDProductGetter const* sim_product_getter);
    
    std::string getStepInstanceName(int position_to_get);
    
    void clearStash();
    
    
   private:
                            
    //this is the size of the input GenParticleCollection, and is used to set the size of the stash
    int initialStashSize;
    
    // A helper class to hold the ptrs to the event data
    struct StashPerEventData {
        
        StashPerEventData():
        sims(nullptr),
        status(nullptr),
        trajectories(nullptr),
        simRemapping(nullptr),
        extMonFNALHits(nullptr)
        {
            tvd.first = "";
            tvd.second = nullptr;
        }
        
        std::unique_ptr<SimParticleCollection>                                  sims;
        std::unique_ptr<StatusG4>                                               status;
        std::pair< std::string, std::unique_ptr<StepPointMCCollection> >        tvd;
        std::unique_ptr<MCTrajectoryCollection>                                 trajectories;
        std::unique_ptr<SimParticleRemapping>                                   simRemapping;
        std::unique_ptr<ExtMonFNALSimHitCollection>                             extMonFNALHits;
        //should we use a different data structure here?
        std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >         sensitiveDetectorSteps;
        std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >         cutsSteps;        
    };
    
    std::vector<StashPerEventData> _myVectorOfPerEventData;
    SimParticleCollectionPrinter simParticlePrinter_;

    
};

}  // end namespace mu2e
#endif /* Mu2eG4_EventStash_hh */


