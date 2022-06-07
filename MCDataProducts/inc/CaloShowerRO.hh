#ifndef MCDataProducts_CaloShowerROStep_hh
#define MCDataProducts_CaloShowerROStep_hh

#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include <vector>

namespace mu2e {


   class CaloShowerRO
   {
       public:
          CaloShowerRO(): SiPMID_(-1),step_() {}

          CaloShowerRO(int SiPMID, const art::Ptr<CaloShowerStep>& step, const std::vector<float>& PETime) :
             SiPMID_(SiPMID),step_(step),PETime_(PETime)
          {}

          const art::Ptr<CaloShowerStep>&   caloShowerStep()  const {return step_;}
          const std::vector<float>          PETime()          const {return PETime_;}
          int                               SiPMID()          const {return SiPMID_;}
          unsigned                          NPE()             const {return PETime_.size();}

          void setCaloShowerStep(const art::Ptr<CaloShowerStep>& step) {step_ = step;}

       private:
          int                       SiPMID_;
          art::Ptr<CaloShowerStep>  step_;
          std::vector<float>        PETime_;
   };

   using CaloShowerROCollection = std::vector<mu2e::CaloShowerRO>;
}

#endif
