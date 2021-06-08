#include "Mu2eUtilities/inc/stoppedMuonList.hh"

namespace mu2e {

  //================================================================
  namespace {

    template<class H>
    std::vector<art::Ptr<SimParticle> > stoppedMuonListImpl(H simh) {

      std::vector<art::Ptr<SimParticle> > res;

      const auto& sims = *simh;

      for(auto i = sims.begin(); i != sims.end(); ++i) {
        const auto& inpart = i->second;
        if((inpart.pdgId() == PDGCode::mu_minus) &&
           (inpart.stoppingCode() == ProcessCode::muMinusCaptureAtRest)) // G4 sets this code for both decay and capture cases
          {
            res.emplace_back(simh, i->first.asInt());
          }
      }

      return res;
    }
  }

  //================================================================

  std::vector<art::Ptr<SimParticle> > stoppedMuonList(art::ValidHandle<SimParticleCollection> simh) {
    return stoppedMuonListImpl(simh);
  }

  std::vector<art::Ptr<SimParticle> > stoppedMuonList(art::Handle<SimParticleCollection> simh) {
    return stoppedMuonListImpl(simh);
  }

}
