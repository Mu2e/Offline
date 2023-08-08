#include "Offline/CRVResponse/inc/CrvMCHelper.hh"

namespace mu2e
{
  void CrvMCHelper::GetStepPointsFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                                const art::Handle<CrvDigiMCCollection> &digis,
                                                std::set<art::Ptr<CrvStep> > &steps)
  {
    if(!digis.isValid()) return;

    const std::vector<size_t> &waveformIndices = crvRecoPulse->GetWaveformIndices();
    for(size_t i=0; i<waveformIndices.size(); i++)
    {
      size_t waveformIndex = waveformIndices[i];
      const CrvDigiMC &digi = digis->at(waveformIndex);
      const std::vector<art::Ptr<CrvStep> > &stepPoints = digi.GetCrvSteps();
      for(size_t j=0; j<stepPoints.size(); j++)
      {
        if(stepPoints[j].isNonnull()) steps.insert(stepPoints[j]);
      }
    }
  }

  void CrvMCHelper::GetInfoFromStepPoints(const std::set<art::Ptr<CrvStep> > &steps,
                                        double &visibleEnergyDeposited,
                                        double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                        art::Ptr<SimParticle> &mostLikelySimParticle)
  {
    visibleEnergyDeposited=0;
    std::map<art::Ptr<SimParticle>,double> simParticleMap;
    std::set<art::Ptr<CrvStep> >::const_iterator stepPointIter;
    for(stepPointIter=steps.begin(); stepPointIter!=steps.end(); stepPointIter++)
    {
      const CrvStep &step = **stepPointIter;
      visibleEnergyDeposited+=step.visibleEDep();
      simParticleMap[step.simParticle()]+=step.visibleEDep();
    }

    std::map<art::Ptr<SimParticle>,double>::const_iterator simParticleIter;
    double simParticleDepEnergy=0;
    for(simParticleIter=simParticleMap.begin(); simParticleIter!=simParticleMap.end(); simParticleIter++)
    {
      if(simParticleIter->second>simParticleDepEnergy)
      {
        simParticleDepEnergy=simParticleIter->second;
        mostLikelySimParticle=simParticleIter->first;
      }
    }

    //TODO: Is this still necessary? Removing it for now.
    //time folding is not applied here, but was used to create the digis, ...
    //so we need to avoid that some step points from a different micro bunch
    //could be accidentally found to be the step point with the earliest hit time.
    //therefore, only step points of the most likely sim particle will be considered.
    bool firstLoop=true;
    for(stepPointIter=steps.begin(); stepPointIter!=steps.end(); stepPointIter++)
    {
      const CrvStep &step = **stepPointIter;
      double t = step.startTime();
      if(firstLoop || earliestHitTime>t)
      {
        firstLoop=false;
        earliestHitTime=t;
        earliestHitPos=step.startPosition();
      }
    }
  }

  void CrvMCHelper::GetInfoFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                          const art::Handle<CrvDigiMCCollection> &digis,
                                          double &visibleEnergyDeposited,
                                          double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                          art::Ptr<SimParticle> &mostLikelySimParticle)
  {
    std::set<art::Ptr<CrvStep> > steps;

    CrvMCHelper::GetStepPointsFromCrvRecoPulse(crvRecoPulse, digis, steps);
    CrvMCHelper::GetInfoFromStepPoints(steps, visibleEnergyDeposited,
                                     earliestHitTime, earliestHitPos, mostLikelySimParticle);
  }


} // end namespace mu2e
