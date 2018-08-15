#ifndef MCDataProducts_CaloShowerROStep_hh
#define MCDataProducts_CaloShowerROStep_hh

//
// Propagation of each CaloShowerStep to each readout
//
// Original author Bertrand Echenard
//

#include "MCDataProducts/inc/CaloShowerStep.hh"
#include <iostream>



namespace mu2e {

   class CaloShowerStepRO 
   {
       public:
	  
	  typedef art::Ptr<CaloShowerStep> caloStepPtr;
	  
	  CaloShowerStepRO(): ROID_(-1),step_(),time_(0), energy_(0) {}
	  
	  CaloShowerStepRO(int ROID, const caloStepPtr& step, double time, double energy) : 
	    ROID_(ROID),step_(step),time_(time),energy_(energy)
	  {}


	  int                        ROID()            const {return ROID_;}
	  const caloStepPtr&         caloShowerStep()  const {return step_;}
	  double                     time()            const {return time_;}
	  double                     energy()          const {return energy_;}

          void setCaloShowerStep(const caloStepPtr& step) {step_ = step;}
       private:
	    
	    int                 ROID_;      
	    caloStepPtr         step_;
	    double              time_;          
	    double              energy_;
   };

} 

#endif
