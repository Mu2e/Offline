#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvCoincidenceClusterMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  void CRVAnalysis::FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterModuleLabel,
                                              const std::string &crvCoincidenceClusterMCModuleLabel,
                                              const std::string &crvRecoPulseLabel,
                                              const std::string &crvStepPointMCLabel,
                                              const std::string &mcTrajectoryLabel,
                                              const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo,
                                              CrvSummaryReco &recoSummary, CrvSummaryMC &MCSummary,
                                              CrvPlaneInfoMCCollection &MCInfoPlane, double crvPlaneY)
  {
    art::Handle<CrvCoincidenceClusterCollection>   crvCoincidenceClusterCollection;
    art::Handle<CrvCoincidenceClusterMCCollection> crvCoincidenceClusterMCCollection;
    art::Handle<CrvRecoPulseCollection>            crvRecoPulseCollection;
    art::Handle<StepPointMCCollection>             crvStepPointMCCollection;
    art::Handle<MCTrajectoryCollection>            mcTrajectoryCollection;

    event.getByLabel(crvCoincidenceClusterModuleLabel,"",crvCoincidenceClusterCollection);
    event.getByLabel(crvCoincidenceClusterMCModuleLabel,"",crvCoincidenceClusterMCCollection);
    event.getByLabel(crvRecoPulseLabel,"",crvRecoPulseCollection);
    event.getByLabel(crvStepPointMCLabel,"CRV",crvStepPointMCCollection);
    event.getByLabel(mcTrajectoryLabel,"",mcTrajectoryCollection);

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

    if(!crvRecoPulseCollection.isValid()) return;
    size_t nRecoPulses=crvRecoPulseCollection->size();
    recoSummary._totalPEs=0;
    std::set<CRSScintillatorBarIndex> counters;
    for(size_t i=0; i<nRecoPulses; i++)
    {
      recoSummary._totalPEs+=crvRecoPulseCollection->at(i).GetPEs();
      counters.insert(crvRecoPulseCollection->at(i).GetScintillatorBarIndex());
    }
    recoSummary._nHitCounters=counters.size();


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
    }

    if(crvStepPointMCCollection.isValid())
    {
      size_t nStepPoints=crvStepPointMCCollection->size();
      MCSummary._totalEnergyDeposited=0;
      std::set<CRSScintillatorBarIndex> counters;
      for(size_t i=0; i<nStepPoints; i++)
      {
        MCSummary._totalEnergyDeposited+=crvStepPointMCCollection->at(i).totalEDep();
        counters.insert(crvStepPointMCCollection->at(i).barIndex());
      }
      MCSummary._nHitCounters=counters.size();
    }


    //locate points where the cosmic MC trajectories cross the xz plane of CRV-T
    if(mcTrajectoryCollection.isValid())
    {
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
      for(trajectoryIter=mcTrajectoryCollection->begin(); trajectoryIter!=mcTrajectoryCollection->end(); trajectoryIter++)
      {
        const art::Ptr<SimParticle> &trajectorySimParticle = trajectoryIter->first;
        const art::Ptr<SimParticle> &trajectoryPrimaryParticle = FindPrimaryParticle(trajectorySimParticle);
        GenId genId = trajectoryPrimaryParticle->genParticle()->generatorId();
        if(genId==GenId::cosmicToy || genId==GenId::cosmicDYB || genId==GenId::cosmic || genId==GenId::cosmicCRY)
        {
          const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
          if(points.size()<1) continue;
          CLHEP::Hep3Vector previousPos=points[0].pos();
          for(size_t i=1; i<points.size(); i++)
          {
            CLHEP::Hep3Vector pos=points[i].pos();
            if((previousPos.y()>crvPlaneY && pos.y()<=crvPlaneY) || (previousPos.y()<crvPlaneY && pos.y()>=crvPlaneY))
            {
              double fraction=(crvPlaneY-pos.y())/(previousPos.y()-pos.y());
              CLHEP::Hep3Vector planePos=fraction*(previousPos-pos)+pos;
              double planeTime=fraction*(points[i-1].t()-points[i].t())+points[i].t();
              double planeKineticEnergy=fraction*(points[i-1].kineticEnergy()-points[i].kineticEnergy())+points[i].kineticEnergy();
              MCInfoPlane.emplace_back(trajectorySimParticle->pdgId(), 
                                       trajectoryPrimaryParticle->pdgId(),
                                       trajectoryPrimaryParticle->startMomentum().e(),
                                       trajectoryPrimaryParticle->startPosition(),
                                       planePos,
                                       planeTime,
                                       planeKineticEnergy);
            }
            previousPos=pos;
          }
        }
      }
    }
  }//FillCrvInfoStructure

}
