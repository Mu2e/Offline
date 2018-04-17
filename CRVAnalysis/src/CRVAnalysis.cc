#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  void CRVAnalysis::FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterFinderModuleLabel,
                                              const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo)
  {
    art::Handle<CrvCoincidenceClusters> crvCoincidenceClusters;
    event.getByLabel(crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusters);

    if(crvCoincidenceClusters.product()==NULL) return;

    const std::vector<CrvCoincidenceClusters::Cluster> &clusters = crvCoincidenceClusters->GetClusters();
    std::vector<CrvCoincidenceClusters::Cluster>::const_iterator iter;
    for(iter=clusters.begin(); iter!=clusters.end(); iter++)
    {
      const CrvCoincidenceClusters::Cluster &cluster = *iter;

      //fill the Reco collection
      recoInfo.emplace_back(cluster._crvSectorType, cluster._avgCounterPos, cluster._startTime, cluster._endTime, cluster._PEs, 
                            cluster._crvRecoPulses.size());

      //fill the MC collection
      if(cluster._simParticles.size()>0)
      {
        const art::Ptr<SimParticle> simParticle = cluster._simParticles[0]; //FIXME: checking only the most frequent particle
        MCInfo.emplace_back(true, 
                            simParticle->pdgId(), 
                            simParticle->genParticle()->pdgId(),
                            simParticle->genParticle()->generatorId().id(),
                            cluster._earliestHitPos,
                            cluster._earliestHitTime,
                            cluster._energyDeposited);
      }
    }//loop through all clusters
  }//FillCrvInfoStructure

}
