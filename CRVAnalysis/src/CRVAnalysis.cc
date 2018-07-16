#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterSummaryCollection.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  void CRVAnalysis::FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterSummarizerModuleLabel,
                                              const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo)
  {
    art::Handle<CrvCoincidenceClusterSummaryCollection> crvCoincidenceClusterSummaryCollection;
    event.getByLabel(crvCoincidenceClusterSummarizerModuleLabel,"",crvCoincidenceClusterSummaryCollection);

    if(crvCoincidenceClusterSummaryCollection.product()==NULL) return;

    std::vector<CrvCoincidenceClusterSummary>::const_iterator iter;
    for(iter=crvCoincidenceClusterSummaryCollection->begin(); iter!=crvCoincidenceClusterSummaryCollection->end(); iter++)
    {
      const CrvCoincidenceClusterSummary &cluster = *iter;

      //fill the Reco collection
      recoInfo.emplace_back(cluster.GetCrvSectorType(), cluster.GetAvgCounterPos(), cluster.GetStartTime(), cluster.GetEndTime(), cluster.GetPEs(), 
                            cluster.GetPulses().size());

      //fill the MC collection
      if(cluster.HasMCInfo())
      {
        const art::Ptr<SimParticle> &simParticle = cluster.GetMostLikelySimParticle();
        MCInfo.emplace_back(true, 
                            simParticle->pdgId(), 
0, //                            simParticle->genParticle()->pdgId(),
0, //                            simParticle->genParticle()->generatorId().id(),
                            cluster.GetEarliestHitPos(),
                            cluster.GetEarliestHitTime(),
                            cluster.GetTotalEnergyDeposited());
      }
      else MCInfo.emplace_back();
    }//loop through all clusters
  }//FillCrvInfoStructure

}
