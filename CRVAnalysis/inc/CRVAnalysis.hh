#include <string>

#include "CRVAnalysis/inc/CrvHitInfoReco.hh"
#include "CRVAnalysis/inc/CrvHitInfoMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e
{
  class CRVAnalysis
  {
    public:

    CRVAnalysis() {}

    static void FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterModuleLabel,
                                          const std::string &crvCoincidenceClusterMCModuleLabel,
                                          const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo);

    private:
    static const art::Ptr<SimParticle> &FindPrimaryParticle(const art::Ptr<SimParticle> &simParticle) 
    {
      return simParticle->hasParent() ? FindPrimaryParticle(simParticle->parent()) : simParticle;
    }
  };

}


