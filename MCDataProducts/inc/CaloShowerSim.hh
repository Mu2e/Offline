#ifndef MCDataProducts_CaloShowerSim_hh
#define MCDataProducts_CaloShowerSim_hh

#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include <iostream>
#include <numeric>


namespace 
{
   float sumEnergyDepG4(float x, const art::Ptr<mu2e::CaloShowerStep>& step)                              {return x+step->energyDepG4();}
   bool  momentumOrder (const art::Ptr<mu2e::CaloShowerStep>& a, const art::Ptr<mu2e::CaloShowerStep>& b) {return a->momentumIn()<b->momentumIn();}
   bool  timeOrder     (const art::Ptr<mu2e::CaloShowerStep>& a, const art::Ptr<mu2e::CaloShowerStep>& b) {return a->time()<b->time();}
}


namespace mu2e {
      
   class CaloShowerSim
   {
       public:          
          using SimPtr   = art::Ptr<SimParticle>;
          using StepPtr  = art::Ptr<CaloShowerStep>;
          using StepPtrs = std::vector<StepPtr>;
                    
          CaloShowerSim() : 
             steps_(),energyCorr_(0),timeCorr_(0)
          {}          
          
          CaloShowerSim(const StepPtrs& steps, float energyCorr, float timeCorr) : 
             steps_(steps),energyCorr_(energyCorr),timeCorr_(timeCorr) 
          {}

          
          const SimPtr&    sim()             const {return steps_[0]->simParticle();}
          const StepPtrs&  caloShowerSteps() const {return steps_;}
                int        crystalId()       const {return steps_[0]->volumeId(); }
                float      time()            const {return timeCorr_;}
                float      energyDep()       const {return energyCorr_;}
                float      energyDepG4()     const {return std::accumulate(steps_.begin(),steps_.end(),0.0f,sumEnergyDepG4);}
                float      timeOrig()        const {return (*std::min_element(steps_.begin(),steps_.end(),timeOrder))->time();}
                float      momentumIn()      const {return (*std::max_element(steps_.begin(),steps_.end(),momentumOrder))->momentumIn();}
 
          
          void setCaloShowerSteps(const StepPtrs& steps) {steps_ = steps;}

       private:
            StepPtrs steps_;          
            float    energyCorr_;
            float    timeCorr_;
   };
   
   
   typedef std::vector<mu2e::CaloShowerSim> CaloShowerSimCollection;

} 

#endif

