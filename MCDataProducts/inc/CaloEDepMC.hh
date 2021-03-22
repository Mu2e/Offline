#ifndef MCDataProducts_CaloEDepMC_hh
#define MCDataProducts_CaloEDepMC_hh

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCRelationship.hh"

namespace mu2e {
 
   class CaloEDepMC
   {
       public:                              
          CaloEDepMC() : 
             sim_(),eDep_(0.0f),eDepG4_(1e6),time_(0.0f),momentumIn_(0.0f),rel_() 
          {};
          
          CaloEDepMC(const art::Ptr<SimParticle>& sim, float eDep, float eDepG4, float time, float mom, const MCRelationship& rel) : 
             sim_(sim),eDep_(eDep),eDepG4_(eDepG4),time_(time),momentumIn_(mom),rel_(rel)
          {};
          
          const art::Ptr<SimParticle>&   sim        ()  const  {return sim_;}
          const MCRelationship&          rel        ()  const  {return rel_;}
          float                          energyDep  ()  const  {return eDep_;}
          float                          energyDepG4()  const  {return eDepG4_;}
          float                          time       ()  const  {return time_;}
          float                          momentumIn ()  const  {return momentumIn_;}
          
          void  addEDep  (float val)                           {eDep_+=val;}
          void  addEDepG4(float val)                           {eDepG4_+=val;}
          void  addTime  (float val)                           {time_ = std::min(time_,val);}
          void  addMom   (float val)                           {momentumIn_ = std::max(momentumIn_,val);}
          void  resetSim (const art::Ptr<SimParticle>& sim)    {sim_ = sim;}
                 
                   
       private:
          art::Ptr<SimParticle> sim_;
          float                 eDep_;
          float                 eDepG4_;
          float                 time_;
          float                 momentumIn_; //the momentum of the SimParticle when entering in the disk
          MCRelationship        rel_; 
   };   
} 

#endif
