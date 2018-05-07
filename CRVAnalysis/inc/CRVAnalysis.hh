#include <string>

#include "CRVAnalysis/inc/CrvHitInfoReco.hh"
#include "CRVAnalysis/inc/CrvHitInfoMC.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e
{
  class CRVAnalysis
  {
    public:

    CRVAnalysis() {}

    static void FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterSummarizerModuleLabel,
                                          const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo);

  };

}


