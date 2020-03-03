#include <string>

#include "CRVAnalysis/inc/CrvHitInfoReco.hh"
#include "CRVAnalysis/inc/CrvHitInfoMC.hh"
#include "CRVAnalysis/inc/CrvSummaryReco.hh"
#include "CRVAnalysis/inc/CrvSummaryMC.hh"
#include "CRVAnalysis/inc/CrvPlaneInfoMC.hh"
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
                                          const std::string &crvRecoPulseLabel,
                                          const std::string &crvStepPointMCLabel,
                                          const std::string &mcTrajectoryLabel,
                                          const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo,
                                          CrvSummaryReco &recoSummary, CrvSummaryMC &MCSummary,
                                          CrvPlaneInfoMCCollection &MCInfoPlane, double crvPlaneY);

    private:
    static const art::Ptr<SimParticle> &FindPrimaryParticle(const art::Ptr<SimParticle> &simParticle) 
    {
      return simParticle->hasParent() ? FindPrimaryParticle(simParticle->parent()) : simParticle;
    }
  };

}


