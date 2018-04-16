#include <iostream>
#include <map>
#include <string>

#include "CRVAnalysis/inc/CrvHitInfoReco.hh"
#include "CRVAnalysis/inc/CrvHitInfoMC.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusters.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e
{
  class CRVAnalysis
  {
    public:

    CRVAnalysis() {}

//FIXME
/*
    static void FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterFinderModuleLabel,
                                          const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo);

    private:

    static constexpr double _overlapTime = 20.0; //ns
    static void FillCrvHitInfoMCCollection(const std::vector<CrvCoincidenceClusters::Hit> &hits,
                                           const std::vector<art::Handle<mu2e::StepPointMCCollection> > &CrvStepsVector,
                                           CrvHitInfoMCCollection &MCCInfo);
*/
  };

}


