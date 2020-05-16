#ifndef MCDataProducts_CaloShowerROStep_hh
#define MCDataProducts_CaloShowerROStep_hh
//
// PE generated in the photo-sensors from energy deposit in the crystal
//
// Original author Bertrand Echenard
//
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include <vector>

namespace mu2e {

   class CaloShowerStepRO 
   {
       public:
	 	  
	  CaloShowerStepRO(): ROID_(-1),step_() {}	  
	  CaloShowerStepRO(int ROID, const art::Ptr<CaloShowerStep>& step, const std::vector<float> PETime) : 
	    ROID_(ROID),step_(step),PETime_(PETime) {}

	  int                               ROID()            const {return ROID_;}
	  const art::Ptr<CaloShowerStep>&   caloShowerStep()  const {return step_;}
	  const std::vector<float>          PETime()          const {return PETime_;}
          unsigned                          NPE()             const {return PETime_.size();}
          
          void setCaloShowerStep(const art::Ptr<CaloShowerStep>& step) {step_ = step;}
       
       private:	    
	    int                       ROID_;      
	    art::Ptr<CaloShowerStep>  step_;
	    std::vector<float>        PETime_;          
   };
} 

#endif
