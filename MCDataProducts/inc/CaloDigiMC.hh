#ifndef MCDataProducts_CaloDigiMC_hh
#define MCDataProducts_CaloDigiMC_hh

#include "MCDataProducts/inc/SimParticle.hh"
#include <iostream>
#include <vector>
#include <numeric>

namespace mu2e
{   
   class CaloDigiMC
   {
       public:          
          using SimPtr = art::Ptr<SimParticle>;
                    
          CaloDigiMC() : 
            simParticle_(),eDep_(),eDepG4_(),time_(),momentumIn_(),isConversion_(false) 
          {};
          
          CaloDigiMC(const std::vector<SimPtr>& sims,         const std::vector<float>& eDeps, 
                     const std::vector<float>& eDepG4s,       const std::vector<float>& times, 
                     const std::vector<float>& momentumIns,   bool isConversion) : 
            simParticle_(sims),eDep_(eDeps),eDepG4_(eDepG4s),time_(times),momentumIn_(momentumIns), 
            isConversion_(isConversion) 
         {}
          
          
          unsigned nParticles      ()           const {return simParticle_.size();}
          SimPtr   simParticle     (unsigned i) const {return simParticle_.at(i); }
          float    time            (unsigned i) const {return time_.at(i);}
          float    energyDep       (unsigned i) const {return eDep_.at(i);}
          float    energyDepG4     (unsigned i) const {return eDepG4_.at(i);}
          float    momentumInDisk  (unsigned i) const {return momentumIn_.at(i);}
          float    totalEnergyDep  ()           const {return std::accumulate(eDep_.begin(),eDep_.end(),0.0f);}
          float    totalEnergyDepG4()           const {return std::accumulate(eDepG4_.begin(),eDepG4_.end(),0.0f);}
          bool     isConversion    ()           const {return isConversion_;}
                    

       private:
            std::vector<SimPtr> simParticle_;
            std::vector<float>  eDep_;
            std::vector<float>  eDepG4_;
            std::vector<float>  time_;
            std::vector<float>  momentumIn_; //the momentum of the SimParticle when entering in the disk
            bool                isConversion_;
   };
} 

#endif
