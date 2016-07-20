#include "CRVAnalysis/inc/CRVAnalysis.hh"

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

namespace mu2e
{
  CRVAnalysis::CRVAnalysis(const std::string &g4ModuleLabel, const std::string &crvCoincidenceModuleLabel) : 
                           _g4ModuleLabel(g4ModuleLabel), _crvCoincidenceModuleLabel(crvCoincidenceModuleLabel)
  {
    GeomHandle<DetectorSystem> det;
    _detSysOrigin = det->getOrigin();

    GeomHandle<CosmicRayShield> CRS;
    for(size_t i=0; i<CRS->getCRSScintillatorShields().size(); i++)
    {
      const CRSScintillatorShield &sector = CRS->getCRSScintillatorShield(i);
      int sectorType = sector.getSectorType();
      int coordinate = abs(sector.getCRSScintillatorBarDetail().getThicknessDirection());

      if(sector.nModules()>0)
      {
        const CRSScintillatorModule &module = sector.getModule(0);
        double x0 = module.getAluminumSheet(0).getPosition()[coordinate];
        double x1 = module.getAluminumSheet(1).getPosition()[coordinate];
        double position=(x0+x1)/2.0;
        CrvPlane plane(coordinate,position);
        _crvPlanes[sectorType]=plane;
      }
    }
  }

  bool CRVAnalysis::FindCrvPlaneCrossings(const art::Event& event, const cet::map_vector_key& particleKey, const int &sectorType,
                                          CLHEP::Hep3Vector &point, CLHEP::Hep3Vector &direction)
  {
    const CrvPlane &plane = _crvPlanes[sectorType];
    int coordinate=plane._coordinate;
    double position=plane._position;

    std::string productInstanceName="";
    art::Handle<mu2e::MCTrajectoryCollection> _mcTrajectories;
    if(event.getByLabel(_g4ModuleLabel,productInstanceName,_mcTrajectories))
    {
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator traj_iter;
      for(traj_iter=_mcTrajectories->begin(); traj_iter!=_mcTrajectories->end(); traj_iter++)
      {
        if(traj_iter->first->id()==particleKey) 
        {
          const std::vector<CLHEP::HepLorentzVector> &trajectoryPoints = traj_iter->second.points();
          for(unsigned int i=1; i<trajectoryPoints.size(); i++)
          {
            CLHEP::Hep3Vector point1=trajectoryPoints[i-1];
            CLHEP::Hep3Vector point2=trajectoryPoints[i];
            CLHEP::Hep3Vector diffVector=point2-point1;
            if(diffVector[coordinate]==0) continue;  //these two points are both on the same plane, try to find another pair
            if((point1[coordinate]>=position && point2[coordinate]<=position)
            || (point1[coordinate]<=position && point2[coordinate]>=position))
            {
              double ratio=(position - point1[coordinate])/diffVector[coordinate];
              point=ratio*diffVector+point1-_detSysOrigin;
              direction=diffVector.unit();
              return(true); //don't look for another trajectory
            }
          } //loop over trajectory points
          return(false); //don't look for another trajectory
        }
      }//loop over trajectories
    }
    return(false);
  }

  void CRVAnalysis::FillCRVInfoMCStructure(const art::Event& event, CRVAnalysisInfo &info)
  {
    GeomHandle<CosmicRayShield> CRS;

    std::vector<art::Handle<mu2e::StepPointMCCollection> > CRVStepsVector;
    art::Selector selector(art::ProductInstanceNameSelector("CRV") && 
                           art::ProcessNameSelector("*"));
    event.getMany(selector, CRVStepsVector);

    std::string simParticleProductInstanceName="";
    art::Handle<mu2e::SimParticleCollection> simParticleCollection;
    event.getByLabel(_g4ModuleLabel,simParticleProductInstanceName,simParticleCollection);


    std::string crvCoincidenceInstanceName="";
    art::Handle<CrvCoincidenceCheckResult> crvCoincidenceCheckResult;
    event.getByLabel(_crvCoincidenceModuleLabel,crvCoincidenceInstanceName,crvCoincidenceCheckResult);

    if(crvCoincidenceCheckResult.product()==NULL) return;

    //loop through all coincidence combinations
    std::vector<CoincidenceCluster> coincidenceClusters;

    const std::vector<CrvCoincidenceCheckResult::CoincidenceCombination> &coincidenceCombinations = crvCoincidenceCheckResult->GetCoincidenceCombinations();
    std::vector<CrvCoincidenceCheckResult::CoincidenceCombination>::const_iterator iter;
    for(iter=coincidenceCombinations.begin(); iter!=coincidenceCombinations.end(); iter++)
    {
      std::vector<int> matchingClusters;
      for(int i=0; i<3; i++)
      {
        for(size_t j=0; j<coincidenceClusters.size(); j++)
        for(size_t k=0; k<coincidenceClusters[j].size(); k++)
        {
          const CoincidenceHit &hit=coincidenceClusters[j][k];
          if(hit._counter==iter->_counters[i] && fabs(hit._time-iter->_time[i])<50) matchingClusters.push_back(j);
        }
      }

      if(matchingClusters.size()==0)
      {
        //found no matching cluster: this coincidence triplet is its own cluster
        CoincidenceCluster newCluster;
        for(int i=0; i<3; i++) newCluster.emplace_back(iter->_time[i], iter->_PEs[i], iter->_counters[i]);
        coincidenceClusters.push_back(newCluster);
      }
      if(matchingClusters.size()>=1)
      {
        //found one matching cluster: this coincidence triplet is added to this cluster
        CoincidenceCluster &cluster0=coincidenceClusters[matchingClusters[0]];
        for(int i=0; i<3; i++) cluster0.emplace_back(iter->_time[i], iter->_PEs[i], iter->_counters[i]);
      }
      if(matchingClusters.size()>=2)
      {
        //found two or three matching cluster: additionally add the content of the 2nd (and 3rd cluster), and remove the 2nd (and 3rd cluster)
        if(matchingClusters[0]!=matchingClusters[1])
        {
          CoincidenceCluster &cluster0=coincidenceClusters[matchingClusters[0]];
          CoincidenceCluster &cluster1=coincidenceClusters[matchingClusters[1]];
          cluster0.insert(cluster0.end(),cluster1.begin(),cluster1.end());
          cluster1.clear();  //shrink to size 0 instead of removing in order to preserve the indices
        }
      }
      if(matchingClusters.size()==3)
      {
        if(matchingClusters[0]!=matchingClusters[2])
        {
          CoincidenceCluster &cluster0=coincidenceClusters[matchingClusters[0]];
          CoincidenceCluster &cluster2=coincidenceClusters[matchingClusters[2]];
          cluster0.insert(cluster0.end(),cluster2.begin(),cluster2.end());
          cluster2.clear();  //shrink to size 0 instead of removing in order to preserve the indices
        }
      };
    }//loop through coincidence combinations

    //loop through all coincidence clusters
    std::vector<CoincidenceCluster>::const_iterator iterCluster;
    for(iterCluster=coincidenceClusters.begin(); iterCluster!=coincidenceClusters.end(); iterCluster++)
    {
      const CoincidenceCluster &clusterOriginal=*iterCluster;
      if(clusterOriginal.empty()) continue;

      //remove duplicate hits
      CoincidenceCluster::const_iterator iterHit, iterHitOriginal;
      CoincidenceCluster cluster;
      for(iterHitOriginal=clusterOriginal.begin(); iterHitOriginal!=clusterOriginal.end(); iterHitOriginal++)
      {
        bool duplicateHit=false;
        for(iterHit=cluster.begin(); iterHit!=cluster.end(); iterHit++)
        {
          if(*iterHit==*iterHitOriginal) {duplicateHit=true; continue;} //found a duplicate hit
        } 
        if(!duplicateHit) cluster.push_back(*iterHitOriginal); //TODO use a vector of references or pointers to reduce time for copying
      }

      //loop through all hits of the cluster to find 
      //-the CRV sector
      //-the average position
      //-the first and last hit time
      //-the total number of PEs
      //-the most likely track responsible
      CRVInfoReco infoReco;
      infoReco._PEs=0;  //used also as sum of all weights
      infoReco._timeWindowStart=NAN;
      infoReco._timeWindowEnd=NAN;
      std::map<cet::map_vector_key, int> tracks; 
      for(iterHit=cluster.begin(); iterHit!=cluster.end(); iterHit++)
      {
        const CoincidenceHit &hit = *iterHit;

        //get CRV counter type
        const CRSScintillatorBar &crvCounter = CRS->getBar(hit._counter);
        int sectorNumber = crvCounter.id().getShieldNumber();
        infoReco._crvSectorType = CRS->getCRSScintillatorShield(sectorNumber).getSectorType();
        // 0: R
        // 1: L
        // 2: T
        // 3: D
        // 4: U
        // 5,6,7: C1,C2,C3

        //get CRV average counter position (weighted by PEs)
        infoReco._pos+=crvCounter.getPosition()*hit._PEs;

        //get total number of PEs
        infoReco._PEs+=hit._PEs;

        //get first and last hit
        if(isnan(infoReco._timeWindowStart) || infoReco._timeWindowStart>hit._time) infoReco._timeWindowStart=hit._time;
        if(isnan(infoReco._timeWindowEnd) || infoReco._timeWindowEnd<hit._time) infoReco._timeWindowEnd=hit._time;

        //collect possible StepPoint MC to find potential track IDs and weight them by PEs
        for(size_t i=0; i<CRVStepsVector.size(); i++)  //vector will be empty for non-MC events
        {
          const art::Handle<StepPointMCCollection> &CRVSteps = CRVStepsVector[i];
          StepPointMCCollection::const_iterator iterStepPoint;
          for(iterStepPoint=CRVSteps->begin(); iterStepPoint!=CRVSteps->end(); iterStepPoint++)
          {
            StepPointMC const& step(*iterStepPoint);
            double timeDiff = hit._time-step.time();
            if(timeDiff<0 || timeDiff>80) continue;
            cet::map_vector_key id = step.trackId();
            tracks[id]+=hit._PEs;
          }
        }
      }//loop over all hits in this coincidence cluster

      infoReco._pos/=infoReco._PEs;  //find average position

      info._infoReco.push_back(infoReco);

      if(simParticleCollection.product()!=NULL)
      {
        //find track ID with largest weight
        cet::map_vector_key trackID;  //FIXME uninitialized
        int largestWeight=0;
        std::map<cet::map_vector_key,int>::const_iterator iterTrack;
        for(iterTrack=tracks.begin(); iterTrack!=tracks.end(); iterTrack++)
        {
          if(iterTrack->second>largestWeight) 
          {
            largestWeight=iterTrack->second;
            trackID=iterTrack->first;
          }
        }

        //find out particle type, whether it is of cosmic origin, where it crosses the CRV, and with which direction
        CRVInfoMC infoMC;
        infoMC._cosmicOrigin=false;
        infoMC._validPdgId=false;
        infoMC._validPosDirTime=false;
        cet::map_vector<mu2e::SimParticle>::const_iterator iterSimParticle=simParticleCollection->find(trackID);
        if(iterSimParticle!=simParticleCollection->end())
        {
          const mu2e::SimParticle& particle = iterSimParticle->second;
          infoMC._validPdgId=true;
          infoMC._pdgId=particle.pdgId();
          int pdgIdParent=infoMC._pdgId;
          art::Ptr<SimParticle> parent = particle.parent();
          while(parent.isNonnull())
          {
            pdgIdParent=parent->pdgId();
            parent=parent->parent();
          }
          if(abs(pdgIdParent)==13) infoMC._cosmicOrigin=true;
          if(infoMC._cosmicOrigin)
          {
            infoMC._validPosDirTime=FindCrvPlaneCrossings(event,trackID, infoReco._crvSectorType, infoMC._pos, infoMC._direction);
          }
        }

        info._infoMC.push_back(infoMC);
      }

    }//loop over all coincidence clusters


  }//fillCRVInfoMCStructure
}


