#ifndef MCDataProducts_CaloShowerSim_hh
#define MCDataProducts_CaloShowerSim_hh

//
// Summary information of energy deposited after all correction inside each crystal for a given SimParticle
//
// Original author Bertrand Echenard
//

#include <iostream>
#include <numeric>
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/SimParticle.hh"


namespace mu2e {
   
  
   namespace 
   {
      float sumEnergyDepG4(float x, const art::Ptr<CaloShowerStep>& step)                        {return x+step->energyDepG4();}
      bool  momentumOrder (const art::Ptr<CaloShowerStep>& a, const art::Ptr<CaloShowerStep>& b) {return a->momentumIn()<b->momentumIn();}
      bool  timeOrder     (const art::Ptr<CaloShowerStep>& a, const art::Ptr<CaloShowerStep>& b) {return a->time()<b->time();}
   }
   

   class CaloShowerSim
   {
       public:          
          using SimPtr   = art::Ptr<SimParticle>;
          using StepPtr  = art::Ptr<CaloShowerStep>;
          using StepPtrs = std::vector<StepPtr>;
                    
          CaloShowerSim()                                         : steps_(),energyCorr_(0)               {}          
          CaloShowerSim(const StepPtrs& steps, double energyCorr) : steps_(steps),energyCorr_(energyCorr) {}

          int             crystalId()       const {return steps_[0]->volumeId(); }
          const SimPtr&   sim()             const {return steps_[0]->simParticle();}
          const StepPtrs& caloShowerSteps() const {return steps_;}
          float           energyDep()       const {return energyCorr_;}
          float           energyDepG4()     const {return std::accumulate(steps_.begin(),steps_.end(),0.0f,sumEnergyDepG4);}
          float           time()            const {return (*std::min_element(steps_.begin(),steps_.end(),timeOrder))->time();}
          float           momentumIn()      const {return (*std::max_element(steps_.begin(),steps_.end(),momentumOrder))->momentumIn();}
 
          void setCaloShowerSteps(const StepPtrs& steps) {steps_ = steps;}

       private:
            StepPtrs  steps_;          
            double    energyCorr_;
   };
} 

#endif

