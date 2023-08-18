#include <string>

#include "Offline/CRVAnalysis/inc/CrvHitInfoReco.hh"
#include "Offline/CRVAnalysis/inc/CrvHitInfoMC.hh"
#include "Offline/CRVAnalysis/inc/CrvWaveformInfo.hh"
#include "Offline/CRVAnalysis/inc/CrvSummaryReco.hh"
#include "Offline/CRVAnalysis/inc/CrvSummaryMC.hh"
#include "Offline/CRVAnalysis/inc/CrvPlaneInfoMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/CRVAnalysis/inc/CrvPulseInfoReco.hh"
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
                                          const std::string &crvStepLabel,
                                          const std::string &simParticleLabel,
                                          const std::string &mcTrajectoryLabel,
                                          const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo,
                                          CrvSummaryReco &recoSummary, CrvSummaryMC &MCSummary,
                                          CrvPlaneInfoMCCollection &MCInfoPlane, double crvPlaneY);

    static void FillCrvPulseInfoCollections(const std::string &crvRecoPulseCollection,
                                            const std::string &crvWaveformsModuleLabel,
                                            const std::string &crvDigiModuleLabel,
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
