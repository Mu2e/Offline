#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  void CRVAnalysis::FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterSummarizerModuleLabel,
                                              const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo)
  {
    art::Handle<CrvCoincidenceClustersSummary> crvCoincidenceClustersSummary;
    event.getByLabel(crvCoincidenceClusterSummarizerModuleLabel,"",crvCoincidenceClustersSummary);

    if(crvCoincidenceClustersSummary.product()==NULL) return;

    const std::vector<CrvCoincidenceClustersSummary::Cluster> &clusters = crvCoincidenceClustersSummary->GetClusters();
    std::vector<CrvCoincidenceClustersSummary::Cluster>::const_iterator iter;
    for(iter=clusters.begin(); iter!=clusters.end(); iter++)
    {
      const CrvCoincidenceClustersSummary::Cluster &cluster = *iter;

      //fill the Reco collection
      recoInfo.emplace_back(cluster._crvSectorType, cluster._avgCounterPos, cluster._startTime, cluster._endTime, cluster._PEs, 
                            cluster._pulses.size());

      //fill the MC collection
      if(cluster._hasMCInfo)
      {
        const art::Ptr<SimParticle> &simParticle = cluster._mostLikelySimParticle;
        MCInfo.emplace_back(true, 
                            simParticle->pdgId(), 
                            simParticle->genParticle()->pdgId(),
                            simParticle->genParticle()->generatorId().id(),
                            cluster._earliestHitPos,
                            cluster._earliestHitTime,
                            cluster._totalEnergyDeposited);
      }
      else MCInfo.emplace_back();
    }//loop through all clusters
  }//FillCrvInfoStructure

}
