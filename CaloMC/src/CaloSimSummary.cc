#include "CaloMC/inc/CaloSimSummary.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {

    void CaloSimSummary::add(art::Ptr<SimParticle> sim, float edep, float edepG4, float time, float hitTime, float pIn)
    {
        auto mfind = simMap_.find(sim);                     
        if (mfind != simMap_.end()) mfind->second.update(edep,edepG4,time,hitTime,pIn);
        else simMap_.insert(std::pair<art::Ptr<SimParticle>,CaloSimSummaryEntry>(sim,CaloSimSummaryEntry(edep,edepG4,time,pIn)));         
    }

    const std::vector<art::Ptr<SimParticle>> CaloSimSummary::sims() const 
    {
        std::vector<art::Ptr<SimParticle>> content;
        for (auto& kv : simMap_) content.push_back(kv.first);
        return content;         
    }

    const std::vector<float> CaloSimSummary::energyDep() const 
    {
        std::vector<float> content;
        for (auto& kv : simMap_) content.push_back(kv.second.edep_);
        return content;         
    }

    const std::vector<float> CaloSimSummary::energyDepG4() const 
    {
        std::vector<float> content;
        for (auto& kv : simMap_) content.push_back(kv.second.edepG4_);
        return content;         
    }

    const std::vector<float> CaloSimSummary::time() const 
    {
        std::vector<float> content;
        for (auto& kv : simMap_) content.push_back(kv.second.time_);
        return content;         
    }

    const std::vector<float> CaloSimSummary::momentumIn() const 
    {
        std::vector<float> content;
        for (auto& kv : simMap_) content.push_back(kv.second.mom_);
        return content;         
    }

    bool CaloSimSummary::isConversion() const
    {     
       for (auto& kv : simMap_)
       {
          auto parent(kv.first);
          while (parent->hasParent()) parent = parent->parent();                     
	  if (parent->genParticle() && parent->genParticle()->generatorId().isConversion() ) return true;
       }    		
       return false;
    }
}
