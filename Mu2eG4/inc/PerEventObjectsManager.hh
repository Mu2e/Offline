#ifndef Mu2eG4_PerEventObjectsManager_hh
#define Mu2eG4_PerEventObjectsManager_hh
//
// PerEventObjectsManager.hh *********************
// ***************************************
//
// Author: Lisa Goodeough
// Date: 2017/05/26
//
//


//Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"


//art includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"

#include <memory>


namespace art { class ProductID; }
namespace art { class Event; }

namespace mu2e {
    
    class PhysicsProcessInfo;
    class SimParticleHelper;
    class SimParticlePrimaryHelper;
    
    
class PerEventObjectsManager
{
  public:
    PerEventObjectsManager(const fhicl::ParameterSet& pset, int number_of_threads);
    
    ~PerEventObjectsManager();
    
    
    //called by PrimaryGeneratorAction to create
    //SimParticleHelper and SimParticlePrimaryHelper
    void createSimParticleHelpers(const art::ProductID& sim_part_ID,
                                  const art::Event* evt,
                                  const art::Handle<GenParticleCollection>* genparthandle,
                                  const art::EDProductGetter* sim_product_getter);
    
    //stores the number of the event instance, corresponding to the instance
    //in the stash of the generated particles, that is currently being processed
    void storeEventInstanceNumber(int instance_number);
    
    
    //accessors
    inline SimParticleHelper* getSimParticleHelper() const { return _spHelper.get(); }
    inline SimParticlePrimaryHelper* getSimParticlePrimaryHelper() const { return _parentHelper.get(); }
    
    inline int getEventInstanceNumber() const {return eventInstanceNumber; }

    
   private:
    
    Mu2eG4MultiStageParameters multiStagePars_;
    //SimParticleHelper *_spHelper;
    //SimParticlePrimaryHelper *_parentHelper;
    
    //::unique_ptr<SimParticleHelper> _spHelper;//causes compile error in AI

    std::shared_ptr<SimParticleHelper> _spHelper;
    std::shared_ptr<SimParticlePrimaryHelper> _parentHelper;
    
    //for managing the stashes
    int eventInstanceNumber;
    
};

}  // end namespace mu2e
#endif /* Mu2eG4_PerEventObjectsManager_hh */


