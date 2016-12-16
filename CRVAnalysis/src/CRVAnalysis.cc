#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  CRVAnalysis::CRVAnalysis(const std::string &crvCoincidenceModuleLabel, double leadingTime, double trailingTime) : 
                           _crvCoincidenceModuleLabel(crvCoincidenceModuleLabel),
                           _leadingTime(leadingTime),
                           _trailingTime(trailingTime)
  {
  }

  void CRVAnalysis::FillCRVHitsCollections(const art::Event& event, CRVHitsRecoCollection &recoInfo, CRVHitsMCCollection &MCInfo)
  {
    GeomHandle<CosmicRayShield> CRS;

    std::vector<art::Handle<mu2e::StepPointMCCollection> > CRVStepsVector;
    art::Selector selector(art::ProductInstanceNameSelector("CRV") && 
                           art::ProcessNameSelector("*"));
    event.getMany(selector, CRVStepsVector);

    std::string crvCoincidenceInstanceName="";
    art::Handle<CrvCoincidenceCheckResult> crvCoincidenceCheckResult;
    event.getByLabel(_crvCoincidenceModuleLabel,crvCoincidenceInstanceName,crvCoincidenceCheckResult);

    if(crvCoincidenceCheckResult.product()==NULL) return;

    //loop through all coincidence combinations
    //extract all the hits of all coincidence combinations
    //distribute them into the crv sector types
    //and put them into a (time ordered) set
    std::map<int, std::set<CrvCoincidenceCheckResult::CoincidenceHit> > sectorTypeMap;

    const std::vector<CrvCoincidenceCheckResult::CoincidenceCombination> &coincidenceCombinationsAll = crvCoincidenceCheckResult->GetCoincidenceCombinations();
    std::vector<CrvCoincidenceCheckResult::CoincidenceCombination>::const_iterator iter;
    for(iter=coincidenceCombinationsAll.begin(); iter!=coincidenceCombinationsAll.end(); iter++)
    {
      for(int i=0; i<3; i++)
      {
        const mu2e::CRSScintillatorBarIndex &crvBarIndex = iter->_counters[i]; 
        const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
        int crvSectorNumber = crvCounter.id().getShieldNumber();
        int crvSectorType = CRS->getCRSScintillatorShield(crvSectorNumber).getSectorType();
        // 0: R
        // 1: L
        // 2: T
        // 3: D
        // 4: U
        // 5,6,7: C1,C2,C3

        sectorTypeMap[crvSectorType].emplace(iter->_time[i], iter->_PEs[i], iter->_counters[i], iter->_SiPMs[i]);
      }
    }

    //loop through all crv sectors types
    std::map<int, std::set<CrvCoincidenceCheckResult::CoincidenceHit> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); sectorTypeMapIter++)
    {
      int crvSectorType = sectorTypeMapIter->first;
      const std::set<CrvCoincidenceCheckResult::CoincidenceHit> &crvHits = sectorTypeMapIter->second;

      //loop through the set of crv hits for this particular crv sector type
      //remember: the hits in this set are time ordered
      std::set<CrvCoincidenceCheckResult::CoincidenceHit>::const_iterator h=crvHits.begin();
      while(h!=crvHits.end())
      {
        int nCoincidenceHits=1;
        int PEs=h->_PEs;
        CLHEP::Hep3Vector pos=CRS->getBar(h->_counter).getPosition();
        double startTime = h->_time - _leadingTime;
        double endTime = h->_time + _trailingTime;

        std::vector<const CrvCoincidenceCheckResult::CoincidenceHit*> cluster;
        cluster.push_back(&(*h));

        while(++h != crvHits.end())  //go to next (time ordered) hit
        {
          if(endTime >= h->_time - _leadingTime)  //the veto time window of this next hit overlaps the time window of the previous hit
          {                                       //so it can be added to the cluster
            endTime = h->_time + _trailingTime;
            nCoincidenceHits++;
            PEs+=h->_PEs;
            pos+=CRS->getBar(h->_counter).getPosition();

            cluster.push_back(&(*h));
          }
          else break;
        }

        pos/=nCoincidenceHits; //find the average position of the crv counters

        //insert the cluster information into the vector of the reco hits
        recoInfo.emplace_back(crvSectorType, pos, startTime, endTime, PEs, nCoincidenceHits);

        //fill the MC collection
        FillCRVHitsMCCollection(cluster, CRVStepsVector, MCInfo);

      }//loop through all hits
    }//loop through all sector types
  }//FillCrvInfoStructure

  //fill the MC collection
  void CRVAnalysis::FillCRVHitsMCCollection(const std::vector<const CrvCoincidenceCheckResult::CoincidenceHit*> &cluster,
                                            const std::vector<art::Handle<mu2e::StepPointMCCollection> > &CRVStepsVector,
                                            CRVHitsMCCollection &MCInfo) 
  {

    struct trackInfo
    {
      double totalEnergy;
      double earliestTime;
      CLHEP::Hep3Vector pos;
      CLHEP::Hep3Vector momentum;
      int pdgId, primaryPdgId;
      std::string generator;
    };
    std::map<cet::map_vector_key, trackInfo> tracks; 
    std::map<cet::map_vector_key, trackInfo>::iterator iterTrack, theTrack;

    //find all step points for the cluster
    for(size_t i=0; i<CRVStepsVector.size(); i++)  //vector will be empty for non-MC events
    {
      const art::Handle<StepPointMCCollection> &CRVStepsCollection = CRVStepsVector[i];
      StepPointMCCollection::const_iterator iterStepPoint;
      for(iterStepPoint=CRVStepsCollection->begin(); iterStepPoint!=CRVStepsCollection->end(); iterStepPoint++)
      {
        const StepPointMC &step = *iterStepPoint;
        for(size_t j=0; j<cluster.size(); j++)
        {
          double timeDiff = cluster[j]->_time-step.time();
          if(timeDiff<100 && timeDiff>0 && step.barIndex()==cluster[j]->_counter)
          {  
            //found a relevant step point for this cluster
            cet::map_vector_key trackID=step.trackId();

            iterTrack = tracks.find(trackID);
            if(iterTrack==tracks.end())
            {
              trackInfo t;
              t.totalEnergy=step.totalEDep();
              t.earliestTime=step.time();
              t.pos=step.position();
              t.momentum=step.momentum();
              t.pdgId=step.simParticle()->pdgId();
              t.primaryPdgId=step.simParticle()->genParticle()->pdgId();
              t.generator=step.simParticle()->genParticle()->generatorId().name();
              tracks[trackID]=t;
            }
            else
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
    else MCInfo.emplace_back(false, 0, 0, "", CLHEP::Hep3Vector(), CLHEP::Hep3Vector(), 0, 0);

  }//FillCRVHitsMCCollection 

}


