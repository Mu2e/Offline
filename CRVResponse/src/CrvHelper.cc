#include "CRVResponse/inc/CrvHelper.hh"

namespace mu2e 
{
  void CrvHelper::GetStepPointsFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse,
                                                const art::Handle<CrvDigiMCCollection> &digis,
                                                std::set<art::Ptr<StepPointMC> > &steps)
  {
    if(!digis.isValid()) return;

    const std::vector<size_t> &waveformIndices = crvRecoPulse->GetWaveformIndices();
    for(size_t i=0; i<waveformIndices.size(); i++) 
    {
      size_t waveformIndex = waveformIndices[i];
      const CrvDigiMC &digi = digis->at(waveformIndex);
      const std::vector<art::Ptr<StepPointMC> > &stepPoints = digi.GetStepPoints();
      for(size_t j=0; j<stepPoints.size(); j++)
      {
        if(stepPoints[j].isNonnull()) steps.insert(stepPoints[j]);
      }
    }
  }

  void CrvHelper::GetInfoFromStepPoints(const std::set<art::Ptr<StepPointMC> > &steps, 
                                        const SimParticleTimeOffset &timeOffsets,
                                        double &totalEnergyDeposited, double &ionizingEnergyDeposited,
                                        double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                        art::Ptr<SimParticle> &mostLikelySimParticle)
  {
    totalEnergyDeposited=0;
    ionizingEnergyDeposited=0;
    std::map<art::Ptr<SimParticle>,double> simParticleMap;
    std::set<art::Ptr<StepPointMC> >::const_iterator stepPointIter;
    for(stepPointIter=steps.begin(); stepPointIter!=steps.end(); stepPointIter++)
    {
      const StepPointMC &step = **stepPointIter;
      totalEnergyDeposited+=step.totalEDep();
      ionizingEnergyDeposited+=step.ionizingEdep();
      simParticleMap[step.simParticle()]+=step.totalEDep();
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
    earliestHitTime=NAN;
    for(stepPointIter=steps.begin(); stepPointIter!=steps.end(); stepPointIter++)
    {
      const StepPointMC &step = **stepPointIter;
      double t = timeOffsets.timeWithOffsetsApplied(step);
      if(isnan(earliestHitTime) || earliestHitTime>t)
      {
        earliestHitTime=t;
        earliestHitPos=step.position();
      }
    }
  }

  void CrvHelper::GetInfoFromCrvRecoPulse(const art::Ptr<CrvRecoPulse> &crvRecoPulse, 
                                          const art::Handle<CrvDigiMCCollection> &digis,
                                          const SimParticleTimeOffset &timeOffsets,
                                          double &totalEnergyDeposited, double &ionizingEnergyDeposited,
                                          double &earliestHitTime, CLHEP::Hep3Vector &earliestHitPos,
                                          art::Ptr<SimParticle> &mostLikelySimParticle)
  {
    std::set<art::Ptr<StepPointMC> > steps;

    CrvHelper::GetStepPointsFromCrvRecoPulse(crvRecoPulse, digis, steps);
    CrvHelper::GetInfoFromStepPoints(steps, timeOffsets, totalEnergyDeposited, ionizingEnergyDeposited, 
                                     earliestHitTime, earliestHitPos, mostLikelySimParticle);
  }

  void CrvHelper::GetCrvCounterInfo(const GeomHandle<CosmicRayShield> &CRS, 
                                    mu2e::CRSScintillatorBarIndex crvBarIndex,
                                    int &sectorNumber, int &moduleNumber, int &layerNumber, int &counterNumber)
  {
    const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
    const CRSScintillatorBarId &crvCounterId = crvCounter.id();

    counterNumber=crvCounterId.getBarNumber();
    layerNumber  =crvCounterId.getLayerNumber();
    moduleNumber =crvCounterId.getModuleNumber();
    sectorNumber =crvCounterId.getShieldNumber();
  }

  std::string CrvHelper::GetSectorName(const GeomHandle<CosmicRayShield> &CRS, int sectorNumber)
  {
    const CRSScintillatorShield &sector = CRS->getCRSScintillatorShields().at(sectorNumber);
    return sector.getName();
  }

  int CrvHelper::GetSectorType(const GeomHandle<CosmicRayShield> &CRS, int sectorNumber)
  {
    const CRSScintillatorShield &sector = CRS->getCRSScintillatorShields().at(sectorNumber);
    return sector.getSectorType();
  }

  CLHEP::Hep3Vector CrvHelper::GetCrvCounterPos(const GeomHandle<CosmicRayShield> &CRS, 
                                     mu2e::CRSScintillatorBarIndex crvBarIndex)
  {
    const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
    return crvCounter.getPosition();
  }

} // end namespace mu2e
