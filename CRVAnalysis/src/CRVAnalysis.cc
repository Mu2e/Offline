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
    GeomHandle<CosmicRayShield> CRS;

    std::vector<art::Handle<mu2e::StepPointMCCollection> > CrvStepsVector;
    art::Selector selector(art::ProductInstanceNameSelector("CRV") && 
                           art::ProcessNameSelector("*"));
    event.getMany(selector, CrvStepsVector);

    std::string crvCoincidenceInstanceName="";
    art::Handle<CrvCoincidenceClusters> crvCoincidenceClusters;
    event.getByLabel(crvCoincidenceClusterFinderModuleLabel,crvCoincidenceInstanceName,crvCoincidenceClusters);

    if(crvCoincidenceClusters.product()==NULL) return;

    const std::vector<CrvCoincidenceClusters::Cluster> &clusters = crvCoincidenceClusters->GetClusters();
    std::vector<CrvCoincidenceClusters::Cluster>::const_iterator iter = clusters.begin();
    for( ; iter!=clusters.end(); iter++)
    {
      //fill the Reco collection
      recoInfo.emplace_back(iter->_crvSectorType, iter->_avgPos, iter->_startTime, iter->_endTime, iter->_PEs, iter->_hits.size());
      //fill the MC collection
      FillCrvHitInfoMCCollection(iter->_hits, CrvStepsVector, MCInfo);
    }//loop through all clusters
  }//FillCrvInfoStructure

  //fill the MC collection
  void CRVAnalysis::FillCrvHitInfoMCCollection(const std::vector<CrvCoincidenceClusters::Hit> &hits,
                                               const std::vector<art::Handle<mu2e::StepPointMCCollection> > &CrvStepsVector,
                                               CrvHitInfoMCCollection &MCInfo) 
  {

    struct trackInfo
    {
      double totalEnergy;
      double earliestTime;
      CLHEP::Hep3Vector pos;
      CLHEP::Hep3Vector momentum;
      int pdgId, primaryPdgId;
      int generator;
    };
    std::map<cet::map_vector_key, trackInfo> tracks; 
    std::map<cet::map_vector_key, trackInfo>::iterator iterTrack, theTrack;

    //find all step points for the cluster
    for(size_t i=0; i<CrvStepsVector.size(); i++)  //vector will be empty for non-MC events
    {
      const art::Handle<StepPointMCCollection> &CrvStepsCollection = CrvStepsVector[i];
      StepPointMCCollection::const_iterator iterStepPoint;
      for(iterStepPoint=CrvStepsCollection->begin(); iterStepPoint!=CrvStepsCollection->end(); iterStepPoint++)
      {
        const StepPointMC &step = *iterStepPoint;
        for(size_t j=0; j<hits.size(); j++)
        {
          double timeDiff = hits[j]._time-step.time();
          if(timeDiff<100 && timeDiff>0 && step.barIndex()==hits[j]._counter)
          {  
            //found a relevant step point for this cluster
            cet::map_vector_key trackID=step.trackId();

            iterTrack = tracks.find(trackID);
            if(iterTrack==tracks.end())  //this is a new track
            {
              trackInfo t;
              t.totalEnergy=step.totalEDep();
              t.earliestTime=step.time();
              t.pos=step.position();
              t.momentum=step.momentum();

              art::Ptr<SimParticle> simparticle = step.simParticle();
              t.pdgId=step.simParticle()->pdgId();

              //trying to find the primary
              while(simparticle->isSecondary()) simparticle=simparticle->parent();
              if(simparticle->isPrimary())
              {
                t.primaryPdgId=simparticle->genParticle()->pdgId();
                t.generator=simparticle->genParticle()->generatorId().id();
              }
              else
              {
                t.primaryPdgId=0;
                t.generator=0;
              }
              tracks[trackID]=t;
            }
            else  //this track has already been collected; just update the deposited energy, and the earliest hit time
            {
              iterTrack->second.totalEnergy+=step.totalEDep();
              if(iterTrack->second.earliestTime>step.time())
              {
                iterTrack->second.earliestTime=step.time();
                iterTrack->second.pos=step.position();
                iterTrack->second.momentum=step.momentum().unit();
              }
            }
          }
        }
      }
    }

    //loop through all tracks of this cluster to find the one which deposited the most energy
    double maxEnergy=0;
    double depositedEnergy=0;
    for(iterTrack=tracks.begin(); iterTrack!=tracks.end(); iterTrack++)
    {
      depositedEnergy+=iterTrack->second.totalEnergy;
      if(iterTrack->second.totalEnergy>maxEnergy)
      {
        maxEnergy=iterTrack->second.totalEnergy;
        theTrack=iterTrack;
      }
    }

    if(maxEnergy>0) 
      MCInfo.emplace_back(true, 
                          theTrack->second.pdgId, 
                          theTrack->second.primaryPdgId, 
                          theTrack->second.generator, 
                          theTrack->second.pos, 
                          theTrack->second.momentum, 
                          theTrack->second.earliestTime, 
                          depositedEnergy);
    else MCInfo.emplace_back(CrvHitInfoMC());

  }//FillCrvHitInfoMCCollection 

}


