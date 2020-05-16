#ifndef CaloSimSummary_HH_
#define CaloSimSummary_HH_

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include <map>
#include <vector>


namespace mu2e {
     
     struct CaloSimSummaryEntry
     {
         CaloSimSummaryEntry(float edep, float edepG4, float time, float mom) : 
           edep_(edep), edepG4_(edepG4), time_(time),mom_(mom) 
         {};	       
         
         void update(float edep, float edepG4, float time, float hitTime, float mom)
         {
            edep_ += edep; edepG4_ += edepG4;mom_ = std::max(mom, mom_);
            if (std::abs(hitTime-time)<std::abs(hitTime-time_)) time_ = time;
         }
         
         float edep_,edepG4_,time_,mom_;
     };


     class CaloSimSummary
     {
         public:
             CaloSimSummary(): simMap_() {};
         
             void add(art::Ptr<SimParticle> sim, float edep, float edepG4, float time, float hitTime, float pIn);
             
             const std::vector<art::Ptr<SimParticle>> sims()           const;
             const std::vector<float>                 energyDep()      const;
             const std::vector<float>                 energyDepG4()    const;
             const std::vector<float>                 time()           const;
             const std::vector<float>                 momentumIn()     const;
             bool                                     isConversion()   const;

         private:
             std::map<art::Ptr<SimParticle>, CaloSimSummaryEntry> simMap_;         
     };
} 

#endif
