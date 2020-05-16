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
      float sumEnergyDepG4(float x, const art::Ptr<CaloShowerStep>& step) {return x+step->energyDepG4();}
      bool  momentumOrder(const art::Ptr<CaloShowerStep>& a, const art::Ptr<CaloShowerStep>& b) {return a->momentumIn()<b->momentumIn();}
   }
   

   class CaloShowerSim {

       public:          
          using SimPtr   = art::Ptr<SimParticle>;
          using StepPtr  = art::Ptr<CaloShowerStep>;
          using StepPtrs = std::vector<StepPtr>;
                    
          CaloShowerSim()                                                  : steps_(),time_(0),energy_(0) {}          
          CaloShowerSim(const StepPtrs& steps, double time, double energy) : steps_(steps),time_(time),energy_(energy) {}

          int             crystalId()       const {return steps_[0]->volumeId(); }
          const SimPtr&   sim()             const {return steps_[0]->simParticle();}
          const StepPtrs& caloShowerSteps() const {return steps_;}
          float           time()            const {return time_;}
          float           energyDep()       const {return energy_;}
          float           energyDepG4()     const {return std::accumulate(steps_.begin(),steps_.end(),0.0f,sumEnergyDepG4);}
          float           momentumIn()      const {return (*std::max_element(steps_.begin(),steps_.end(),momentumOrder))->momentumIn();}
 
          void setCaloShowerSteps(const StepPtrs& steps) {steps_ = steps;}

       private:
            StepPtrs  steps_;          
            double    time_;          
            double    energy_;
   };
} 

#endif

