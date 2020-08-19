#include <string>

#include "CRVAnalysis/inc/CrvHitInfoReco.hh"
#include "CRVAnalysis/inc/CrvHitInfoMC.hh"
#include "CRVAnalysis/inc/CrvWaveformInfo.hh"
#include "CRVAnalysis/inc/CrvSummaryReco.hh"
#include "CRVAnalysis/inc/CrvSummaryMC.hh"
#include "CRVAnalysis/inc/CrvPlaneInfoMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "CRVAnalysis/inc/CrvPulseInfoReco.hh"
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
                                          const std::string &simParticleLabel,
                                          const std::string &mcTrajectoryLabel,
                                          const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo,
                                          CrvSummaryReco &recoSummary, CrvSummaryMC &MCSummary,
                                          CrvPlaneInfoMCCollection &MCInfoPlane, double crvPlaneY);

    static void FillCrvPulseInfoCollections(const std::string &crvRecoPulseCollection,
                                            const std::string &crvWaveformsModuleLabel,
                                            const std::string &crvDigiModuleLabel,
                                            const SimParticleTimeOffset &timeOffsets,
                                            const art::Event& event, CrvPulseInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo, CrvWaveformInfoCollection &waveformInfo);

    private:
    static const art::Ptr<SimParticle> &FindPrimaryParticle(const art::Ptr<SimParticle> &simParticle)
    {
      return simParticle->hasParent() ? FindPrimaryParticle(simParticle->parent()) : simParticle;
    }
    static const SimParticle &FindPrimaryParticle(const SimParticle &simParticle)
    {
      return simParticle.hasParent() ? *FindPrimaryParticle(simParticle.parent()) : simParticle;
    }
    static const art::Ptr<SimParticle> &FindParentParticle(const art::Ptr<SimParticle> &simParticle)
    {
      return simParticle->hasParent() ? simParticle->parent() : simParticle;
    }
    static const art::Ptr<SimParticle> &FindGParentParticle(const art::Ptr<SimParticle> &simParticle)
    {
      const art::Ptr<SimParticle> &parentParticle = simParticle->hasParent() ? simParticle->parent() : simParticle;
      return parentParticle->hasParent() ? parentParticle->parent() : parentParticle;
    }

    static const int _trajectorySimParticleId = 300001;  //only temporarily here for some tests
  };

}
