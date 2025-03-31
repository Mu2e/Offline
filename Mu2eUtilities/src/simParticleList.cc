#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

namespace mu2e {

  // Convert the SimParticleCollection into a list without a selection
  std::vector<art::Ptr<SimParticle> > simParticleList(art::ValidHandle<SimParticleCollection> simh)
  {
    std::vector<art::Ptr<SimParticle> > res;

    const auto& sims = *simh;
    for(auto i = sims.begin(); i != sims.end(); ++i) {
      res.emplace_back(simh, i->first.asInt());
    }

    return res;
  }

  std::vector<art::Ptr<SimParticle> > simParticleList(art::ValidHandle<SimParticleCollection> simh,
                                                      PDGCode::type pdgId)
  {
    std::vector<art::Ptr<SimParticle> > res;

    const auto& sims = *simh;
    for(auto i = sims.begin(); i != sims.end(); ++i) {
      const auto& inpart = i->second;
      if(inpart.pdgId() == pdgId)
        {
          res.emplace_back(simh, i->first.asInt());
        }
    }

    return res;
  }

  std::vector<art::Ptr<SimParticle> > simParticleList(art::ValidHandle<SimParticleCollection> simh,
                                                      PDGCode::type pdgId,
                                                      ProcessCode stoppingCode)
  {
    std::vector<art::Ptr<SimParticle> > res;

    const auto& sims = *simh;
    for(auto i = sims.begin(); i != sims.end(); ++i) {
      const auto& inpart = i->second;
      if((inpart.pdgId() == pdgId) &&
         (inpart.stoppingCode() == stoppingCode))
        {
          res.emplace_back(simh, i->first.asInt());
        }
    }

    return res;
  }

}
