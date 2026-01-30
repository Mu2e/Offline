#include "Offline/CRVResponse/inc/CrvMCHelper.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

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
                                        double &avgHitTime, CLHEP::Hep3Vector &avgHitPos,
                                        art::Ptr<SimParticle> &mostLikelySimParticle,
                                          const art::Handle<mu2e::EventWindowMarker>& ewmh)
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

    //get the event window length from the global constants
    if (!ewmh.isValid()) {
      throw cet::exception("CrvMCHelper") << "EventWindowMarker handle is not valid" << std::endl;
    }
    const bool onspill = ewmh->spillType() == EventWindowMarker::onspill;
    const double ewm_window_length = ewmh->eventLength(); // get the offspill window length from the EWM, ~100 us
    const double event_window_length = (onspill) ? GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod() : ewm_window_length;

    //TODO: Is this still necessary? Removing it for now.
    //time folding is not applied here, but was used to create the digis, ...
    //so we need to avoid that some step points from a different micro bunch
    //could be accidentally found to be the step point with the earliest hit time.
    //therefore, only step points of the most likely sim particle will be considered.
    bool firstLoop=true;
    avgHitTime=0;
    avgHitPos.set(0,0,0);
    for(stepPointIter=steps.begin(); stepPointIter!=steps.end(); stepPointIter++)
    {
      const CrvStep &step = **stepPointIter;
      double t_s = step.startTime();
      double t_e = step.endTime();
      if(firstLoop || earliestHitTime>t_s)
      {
        firstLoop=false;
        earliestHitTime=t_s;
        earliestHitPos=step.startPosition();
      }
      // Need to account for time wrapping when calculating average time
      avgHitTime+=std::fmod(0.5*(t_s+t_e), event_window_length)*step.visibleEDep();
      avgHitPos+=0.5*(step.startPosition()+step.endPosition())*step.visibleEDep();
    }
    if(visibleEnergyDeposited>0.0)
    {
      avgHitTime/=visibleEnergyDeposited;
      avgHitPos/=visibleEnergyDeposited;
    }
  }

  void CrvMCHelper::GetInfoFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                          const art::Handle<CrvDigiMCCollection> &digis,
                                          double &visibleEnergyDeposited,
                                          double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                          double &avgHitTime, CLHEP::Hep3Vector &avgHitPos,
                                          art::Ptr<SimParticle> &mostLikelySimParticle,
                                          const art::Handle<mu2e::EventWindowMarker>& ewmh)
  {
    std::set<art::Ptr<CrvStep> > steps;

    CrvMCHelper::GetStepPointsFromCrvRecoPulse(crvRecoPulse, digis, steps);
    CrvMCHelper::GetInfoFromStepPoints(steps, visibleEnergyDeposited,
                                       earliestHitTime, earliestHitPos, avgHitTime, avgHitPos, mostLikelySimParticle, ewmh);
  }


} // end namespace mu2e
