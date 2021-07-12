#include "Mu2eUtilities/inc/simParticleList.hh"

namespace mu2e {

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
