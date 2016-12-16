#include <iostream>
#include <map>
#include <string>

#include "CRVAnalysis/inc/CRVAnalysisInfo.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e
{
  class CRVAnalysis
  {
    public:

    CRVAnalysis(const std::string &crvCoincidenceModuleLabel, double leadingTime, double trailingTime);

    void FillCRVHitsCollections(const art::Event& event, CRVHitsRecoCollection &recoInfo, CRVHitsMCCollection &MCInfo);

    private:

    std::string _crvCoincidenceModuleLabel;
    double      _leadingTime;
    double      _trailingTime;

    void FillCRVHitsMCCollection(const std::vector<const CrvCoincidenceCheckResult::CoincidenceHit*> &cluster,
                                 const std::vector<art::Handle<mu2e::StepPointMCCollection> > &CRVStepsVector,
                                 CRVHitsMCCollection &MCInfo);

  };

}


