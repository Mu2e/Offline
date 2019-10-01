#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvCoincidenceClusterMCCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  void CRVAnalysis::FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterModuleLabel,
                                              const std::string &crvCoincidenceClusterMCModuleLabel,
                                              const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo)
  {
    art::Handle<CrvCoincidenceClusterCollection>   crvCoincidenceClusterCollection;
    art::Handle<CrvCoincidenceClusterMCCollection> crvCoincidenceClusterMCCollection;

    event.getByLabel(crvCoincidenceClusterModuleLabel,"",crvCoincidenceClusterCollection);
    event.getByLabel(crvCoincidenceClusterMCModuleLabel,"",crvCoincidenceClusterMCCollection);

    if(!crvCoincidenceClusterCollection.isValid()) return;
    size_t nClusters=crvCoincidenceClusterCollection->size();

    for(size_t i=0; i<nClusters; i++)
    {
      const CrvCoincidenceCluster &cluster = crvCoincidenceClusterCollection->at(i);

      //fill the Reco collection
      recoInfo.emplace_back(cluster.GetCrvSectorType(), cluster.GetAvgCounterPos(), 
                            cluster.GetStartTime(), cluster.GetEndTime(), 
                            cluster.GetPEs(), cluster.GetCrvRecoPulses().size());
// test Ptrs
//      for(auto const& pulse : cluster.GetCrvRecoPulses())
//	if(pulse.isNull())
//	  std::cout << "Found invalid CrvRecoPulsePtr" << std::endl;
//	else
//	  std::cout << "CrvRecoPulse has " << pulse->GetWaveformIndices().size()  << " Digis" << std::endl;
    }

    //fill the MC collection
    if(crvCoincidenceClusterMCCollection.isValid())
    {
      size_t nClustersMC=crvCoincidenceClusterMCCollection->size();
      if(nClusters!=nClustersMC) std::cout<<"The number of MC and reco CRV coincidence clusters does not match!"<<std::endl;
      for(size_t i=0; i<nClustersMC; i++)
      {
        const CrvCoincidenceClusterMC &clusterMC = crvCoincidenceClusterMCCollection->at(i);
        if(clusterMC.HasMCInfo())
        {
          const art::Ptr<SimParticle> &simParticle = clusterMC.GetMostLikelySimParticle();
          const art::Ptr<SimParticle> &primaryParticle = FindPrimaryParticle(simParticle);
          MCInfo.emplace_back(true, 
                              simParticle->pdgId(), 
                              primaryParticle->pdgId(),
                              primaryParticle->startMomentum().e(),
                              primaryParticle->startPosition(),
                              clusterMC.GetEarliestHitPos(),
                              clusterMC.GetEarliestHitTime(),
                              clusterMC.GetTotalEnergyDeposited());
        }
        else MCInfo.emplace_back();
      }
    }//loop through all clusters
  }//FillCrvInfoStructure

}
