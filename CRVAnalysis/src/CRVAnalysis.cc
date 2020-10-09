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
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "CRVResponse/inc/CrvHelper.hh"
#include "ConditionsService/inc/CrvParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

namespace mu2e
{
  void CRVAnalysis::FillCrvHitInfoCollections(const std::string &crvCoincidenceClusterModuleLabel,
                                              const std::string &crvCoincidenceClusterMCModuleLabel,
                                              const std::string &crvRecoPulseLabel,
                                              const std::string &crvStepPointMCLabel,
                                              const std::string &simParticleLabel,
                                              const std::string &mcTrajectoryLabel,
                                              const art::Event& event, CrvHitInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo,
                                              CrvSummaryReco &recoSummary, CrvSummaryMC &MCSummary,
                                              CrvPlaneInfoMCCollection &MCInfoPlane, double crvPlaneY)
  {

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvCoincidenceClusterCollection>   crvCoincidenceClusterCollection;
    art::Handle<CrvCoincidenceClusterMCCollection> crvCoincidenceClusterMCCollection;
    art::Handle<CrvRecoPulseCollection>            crvRecoPulseCollection;
    art::Handle<StepPointMCCollection>             crvStepPointMCCollection;
    art::Handle<SimParticleCollection>             simParticleCollection;
    art::Handle<MCTrajectoryCollection>            mcTrajectoryCollection;

    event.getByLabel(crvCoincidenceClusterModuleLabel,"",crvCoincidenceClusterCollection);
    event.getByLabel(crvCoincidenceClusterMCModuleLabel,"",crvCoincidenceClusterMCCollection);
    event.getByLabel(crvRecoPulseLabel,"",crvRecoPulseCollection);
    event.getByLabel(crvStepPointMCLabel,"CRV",crvStepPointMCCollection);
    event.getByLabel(simParticleLabel,"",simParticleCollection);
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
          const art::Ptr<SimParticle> &parentParticle  = FindParentParticle(simParticle);
          const art::Ptr<SimParticle> &gparentParticle = FindGParentParticle(simParticle);
          MCInfo.emplace_back(true,
                              simParticle->pdgId(),
                              primaryParticle->pdgId(), primaryParticle->startMomentum().e() - primaryParticle->startMomentum().m(), primaryParticle->startPosition(),
                              parentParticle->pdgId(),  parentParticle->startMomentum().e()  - parentParticle->startMomentum().m(),  parentParticle->startPosition(),
                              gparentParticle->pdgId(), gparentParticle->startMomentum().e() - gparentParticle->startMomentum().m(), gparentParticle->startPosition(),
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
      double totalStep[] = {0, 0, 0, 0};
      for(size_t i=0; i<nStepPoints; i++)
      {
        MCSummary._totalEnergyDeposited+=crvStepPointMCCollection->at(i).totalEDep();
        counters.insert(crvStepPointMCCollection->at(i).barIndex());
        const CRSScintillatorBarId &CRVCounterId = CRS->getBar(crvStepPointMCCollection->at(i).barIndex()).id();
        int layer = CRVCounterId.getLayerNumber();
        int pdgId = crvStepPointMCCollection->at(i).simParticle()->pdgId();
        if(abs(pdgId)==13)
          totalStep[layer] = totalStep[layer] + crvStepPointMCCollection->at(i).stepLength();
      }
      MCSummary._nHitCounters=counters.size();
      MCSummary._minPathLayer=*std::min_element(totalStep,totalStep+4);
      MCSummary._maxPathLayer=*std::max_element(totalStep,totalStep+4);
    }

#if 1
    //locate points where the cosmic MC trajectories cross the xz plane of CRV-T
    if(mcTrajectoryCollection.isValid())
    {
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
      for(trajectoryIter=mcTrajectoryCollection->begin(); trajectoryIter!=mcTrajectoryCollection->end(); trajectoryIter++)
      {
        const art::Ptr<SimParticle> &trajectorySimParticle = trajectoryIter->first;
        if(abs(trajectorySimParticle->pdgId())!=13) continue;
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
              CLHEP::Hep3Vector planeDir=(pos-previousPos).unit();
              double planeTime=fraction*(points[i-1].t()-points[i].t())+points[i].t();
              double planeKineticEnergy=fraction*(points[i-1].kineticEnergy()-points[i].kineticEnergy())+points[i].kineticEnergy();
              MCInfoPlane.emplace_back(trajectorySimParticle->pdgId(),
                                       trajectoryPrimaryParticle->pdgId(),
                                       trajectoryPrimaryParticle->startMomentum().e(),
                                       trajectoryPrimaryParticle->startPosition(),
                                       planePos,
                                       planeDir,
                                       planeTime,
                                       planeKineticEnergy,
                                       0);  //unused
            }
            previousPos=pos;
          }
        }
      }
    }
#else   //only a temporary section due to the missing trajectories at the CRV. needs to be deleted soon
    {
      GlobalConstantsHandle<ParticleDataTable> pdt;
      const HepPDT::ParticleData& mu_data = pdt->particle(PDGCode::mu_minus).ref();
      const double _mMu = mu_data.mass().value();

      int pdgId=0;
      CLHEP::Hep3Vector primaryPos, planePos, planeDir;
      double primaryEnergy=NAN;
      double planeTime=NAN;
      double earliestTime=NAN;
      double planeKineticEnergy=0;
      int dataSource=0;
      if(crvStepPointMCCollection.isValid())
      {
        size_t nStepPoints=crvStepPointMCCollection->size();
        for(size_t i=0; i<nStepPoints; i++)
        {
          if(crvStepPointMCCollection->at(i).simParticle()->id().asInt()==0) //only sim particles with ID 0
          {
            if(crvStepPointMCCollection->at(i).time()<earliestTime || isnan(earliestTime))
            {
              pdgId=crvStepPointMCCollection->at(i).simParticle()->pdgId();
              primaryPos=crvStepPointMCCollection->at(i).simParticle()->startPosition();
              primaryEnergy=crvStepPointMCCollection->at(i).simParticle()->startMomentum().e();

              planePos=crvStepPointMCCollection->at(i).position();
              planeDir=crvStepPointMCCollection->at(i).momentum().unit();
              double fraction=(crvPlaneY-planePos.y())/planeDir.y();
              planePos=fraction*planeDir+planePos;

              planeTime=crvStepPointMCCollection->at(i).time();
              earliestTime=planeTime;

              double momentum = crvStepPointMCCollection->at(i).momentum().mag();
              planeKineticEnergy=sqrt(momentum*momentum+_mMu*_mMu)-_mMu;

              dataSource=1;
              if(crvStepPointMCCollection->at(i).position().y()>=crvPlaneY) dataSource=2;
            }
          }
        }
      }
      if(dataSource==0 && mcTrajectoryCollection.isValid())  //StepPoints not found for this event
      {
        std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
        for(trajectoryIter=mcTrajectoryCollection->begin(); trajectoryIter!=mcTrajectoryCollection->end(); trajectoryIter++)
        {
          if(trajectoryIter->second.simid()==_trajectorySimParticleId) //seems to point to the primary
          {
            const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
            if(points.size()>=2 && abs(trajectoryIter->first->pdgId())==13)
            {
              pdgId=trajectoryIter->first->pdgId();
              primaryPos=trajectoryIter->first->startPosition();
              primaryEnergy=trajectoryIter->first->startMomentum().e();
              double fraction=(crvPlaneY-points[1].pos().y())/(points[0].pos().y()-points[1].pos().y());
              planePos=fraction*(points[0].pos()-points[1].pos())+points[1].pos();
              planeDir=(points[1].pos()-points[0].pos()).unit();
              planeTime=fraction*(points[0].t()-points[1].t())+points[1].t();
              planeKineticEnergy=points[0].kineticEnergy();  //use Ekin of first point
              dataSource=3;
            }
            break;
          }
        }
      }
      MCInfoPlane.emplace_back(pdgId,
                               pdgId,
                               primaryEnergy,
                               primaryPos,
                               planePos,
                               planeDir,
                               planeTime,
                               planeKineticEnergy,
                               dataSource);
    }
#endif

  }//FillCrvInfoStructure


   void CRVAnalysis::FillCrvPulseInfoCollections (const std::string &crvRecoPulseCollectionModuleLabel,
                                                  const std::string &crvWaveformsModuleLabel,
                                                  const std::string &crvDigiModuleLabel,
                                                  const SimParticleTimeOffset &timeOffsets,
                                                  const art::Event& event, CrvPulseInfoRecoCollection &recoInfo, CrvHitInfoMCCollection &MCInfo, CrvWaveformInfoCollection &waveformInfo){


     art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
     event.getByLabel(crvRecoPulseCollectionModuleLabel,"",crvRecoPulseCollection);
     if(!crvRecoPulseCollection.isValid()) return;

     mu2e::ConditionsHandle<mu2e::CrvParams> crvPar("ignored");
     double _digitizationPeriod  = crvPar->digitizationPeriod;
     double _recoPulsePedestal  = crvPar->pedestal;
     GeomHandle<CosmicRayShield> CRS;

     // Create SiPM map to extract sequantial SiPM IDs
     const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
     std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter;
     int iSiPM = 0;
     std::map<int,int> sipm_map;
     for(iter=counters.begin(); iter!=counters.end(); iter++)
       {
       const CRSScintillatorBarIndex &barIndex = (*iter)->index();
       for(int SiPM=0; SiPM<4; SiPM++)
         {
         if(!(*iter)->getBarDetail().hasCMB(SiPM%2)) continue;
         sipm_map[barIndex.asInt()*4 + SiPM] = iSiPM;
         iSiPM++;
         }
       }

     // Loop through all reco pulses
     for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
       {
         const art::Ptr<CrvRecoPulse> crvRecoPulse(crvRecoPulseCollection, recoPulseIndex);
         //get information about the counter
         const CRSScintillatorBarIndex &barIndex = crvRecoPulse->GetScintillatorBarIndex();
         int sectorNumber  = -1;
         int moduleNumber  = -1;
         int layerNumber   = -1;
         int counterNumber = -1;
         CrvHelper::GetCrvCounterInfo(CRS, barIndex, sectorNumber, moduleNumber, layerNumber, counterNumber);

         //Reco pulses information
         int SiPM = crvRecoPulse->GetSiPMNumber();
         int SiPMId = sipm_map.find(barIndex.asInt()*4 + SiPM)->second;
         CLHEP::Hep3Vector HitPos = CrvHelper::GetCrvCounterPos(CRS, barIndex);
         recoInfo.emplace_back(HitPos, barIndex.asInt(), sectorNumber, SiPMId,
                               crvRecoPulse->GetPEs(), crvRecoPulse->GetPEsPulseHeight(), crvRecoPulse->GetPulseHeight()+_recoPulsePedestal,
                               crvRecoPulse->GetPulseBeta(), crvRecoPulse->GetPulseFitChi2(), crvRecoPulse->GetPulseTime());

         //MCtruth pulses information
         art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
         if(crvWaveformsModuleLabel!="") event.getByLabel(crvWaveformsModuleLabel,"",crvDigiMCCollection); //this is an optional part for MC information
         double totalEnergyDeposited    = 0;
         double ionizingEnergyDeposited = 0;
         double earliestHitTime         = NAN;
         CLHEP::Hep3Vector earliestHitPos;
         art::Ptr<SimParticle> mostLikelySimParticle;
         //for this reco pulse
         CrvHelper::GetInfoFromCrvRecoPulse(crvRecoPulse, crvDigiMCCollection, timeOffsets,
                                            totalEnergyDeposited, ionizingEnergyDeposited,
                                            earliestHitTime, earliestHitPos, mostLikelySimParticle);

         bool hasMCInfo = (mostLikelySimParticle.isNonnull()?true:false); //MC
         if(hasMCInfo)
           {
             const art::Ptr<SimParticle> &primaryParticle = FindPrimaryParticle(mostLikelySimParticle);
             const art::Ptr<SimParticle> &parentParticle = FindParentParticle(mostLikelySimParticle);
             const art::Ptr<SimParticle> &gparentParticle = FindGParentParticle(mostLikelySimParticle);
             MCInfo.emplace_back(true, mostLikelySimParticle->pdgId(),
                                 primaryParticle->pdgId(), primaryParticle->startMomentum().e() - primaryParticle->startMomentum().m(), primaryParticle->startPosition(),
                                 parentParticle->pdgId(),  parentParticle->startMomentum().e()  - parentParticle->startMomentum().m(),  parentParticle->startPosition(),
                                 gparentParticle->pdgId(), gparentParticle->startMomentum().e() - gparentParticle->startMomentum().m(), gparentParticle->startPosition(),
                                 earliestHitPos, earliestHitTime, ionizingEnergyDeposited);
           }
         else
           MCInfo.emplace_back();
       }

     //    Fill waveforms struct
     art::Handle<mu2e::CrvDigiCollection> crvDigis;
     if(crvDigiModuleLabel!="") event.getByLabel(crvDigiModuleLabel,"",crvDigis);
     for(size_t j=0; j<crvDigis->size(); j++)
       {
         mu2e::CrvDigi const& digi(crvDigis->at(j));
         int _SiPMId = sipm_map.find(digi.GetScintillatorBarIndex().asInt()*4 + digi.GetSiPMNumber())->second;
         for(size_t k=0; k<mu2e::CrvDigi::NSamples; k++)
           waveformInfo.emplace_back(digi.GetADCs()[k], (digi.GetStartTDC()+k)*_digitizationPeriod, _SiPMId);
       }
   } //FillCrvPulseInfoCollections

}
