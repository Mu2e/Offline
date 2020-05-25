#ifndef MCDataProducts_CaloEDepMC_hh
#define MCDataProducts_CaloEDepMC_hh

#include "MCDataProducts/inc/SimParticle.hh"


namespace mu2e {
 
   class CaloEDepMC
   {
       public:          
                    
          CaloEDepMC() : 
            sim_(),eDep_(0.0f),eDepG4_(0.0f),time_(0.0f),momentumIn_(0.0f) 
          {};
          
          CaloEDepMC(const art::Ptr<SimParticle>& sim, float eDep, float eDepG4, float time, float mom) : 
            sim_(sim),eDep_(eDep),eDepG4_(eDepG4),time_(time),momentumIn_(mom)
          {};
                    
          const art::Ptr<SimParticle>& sim       () const {return sim_;}
          float                        eDep      () const {return eDep_;}
          float                        eDepG4    () const {return eDepG4_;}
          float                        time      () const {return time_;}
          float                        momentumIn() const {return momentumIn_;}
          
          void  set(float eDep, float eDepG4, float time, float mom) {eDep_=eDep;eDepG4_=eDepG4;time_=time;momentumIn_=mom;} 
                 
                   
       private:
          art::Ptr<SimParticle> sim_;
          float                 eDep_;
          float                 eDepG4_;
          float                 time_;
          float                 momentumIn_; //the momentum of the SimParticle when entering in the disk
   };   
} 

#endif
