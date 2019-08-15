//
// A module to find the MC information for the clusters of coincidences of CRV pulses
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "MCDataProducts/inc/CrvCoincidenceClusterMCCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>
#include <TH2D.h>

namespace mu2e 
{
  class CrvCoincidenceClusterMatchMC : public art::EDProducer 
  {
    public:
    explicit CrvCoincidenceClusterMatchMC(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    std::string _crvCoincidenceClusterFinderModuleLabel;  //module label of the CrvCoincidenceClusterFinder module
                                                          //it is possible to have more than one instance of the CrvCoincidenceClusterFinder module
    std::string _crvWaveformsModuleLabel;  //module label of the CrvWaveform module. 
                                           //this is optional. only needed, if MC information is required
    SimParticleTimeOffset _timeOffsets;
    void CollectStepPoints(std::set<art::Ptr<StepPointMC> > &steps, 
                           const art::Handle<CrvDigiMCCollection> &digis, 
                           size_t waveformIndex);
    void ScanStepPoints(const std::set<art::Ptr<StepPointMC> > &steps, 
                        double &energyDeposited, double &earliestHitTime,
                        CLHEP::Hep3Vector &earliestHitPos,
                        art::Ptr<SimParticle> &mostLikelySimParticle);
  };

  CrvCoincidenceClusterMatchMC::CrvCoincidenceClusterMatchMC(fhicl::ParameterSet const& pset) :
   art::EDProducer{pset},
    _crvCoincidenceClusterFinderModuleLabel(pset.get<std::string>("crvCoincidenceClusterFinderModuleLabel")),
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel","")),
    _timeOffsets(pset.get<fhicl::ParameterSet>("timeOffsets"))
  {
    produces<CrvCoincidenceClusterMCCollection>();
  }

  void CrvCoincidenceClusterMatchMC::beginJob()
  {
  }

  void CrvCoincidenceClusterMatchMC::endJob()
  {
  }

  void CrvCoincidenceClusterMatchMC::beginRun(art::Run &run)
  {
  }

  void CrvCoincidenceClusterMatchMC::produce(art::Event& event) 
  {
    _timeOffsets.updateMap(event);

    std::unique_ptr<CrvCoincidenceClusterMCCollection> crvCoincidenceClusterMCCollection(new CrvCoincidenceClusterMCCollection);

    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;
    event.getByLabel(_crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusterCollection);

    if(crvCoincidenceClusterCollection.product()==NULL) return;

    art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
    if(_crvWaveformsModuleLabel!="") event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection); //this is an optional part for MC information

    CrvCoincidenceClusterCollection::const_iterator iter;
    for(iter=crvCoincidenceClusterCollection->begin(); iter!=crvCoincidenceClusterCollection->end(); iter++)
    {
      bool   hasMCInfo                = (crvDigiMCCollection.isValid()?true:false); //MC
      double totalEnergyDeposited     = 0;         //MC
      double earliestHitTime          = NAN;       //MC
      art::Ptr<SimParticle> mostLikelySimParticle; //MC
      CLHEP::Hep3Vector     earliestHitPos;        //MC

      //loop through all reco pulses and try to find the MC information
      std::vector<CrvCoincidenceClusterMC::PulseInfo> pulses; //collection all pulses (sim particles, energy dep.)
      std::map<art::Ptr<SimParticle>, int> simParticleMap;  //counting all simParticles
      std::set<art::Ptr<StepPointMC> > stepsAllPulses;  //collecting all step points for total energy calculation
                                            //use a set to avoid double counts, e.g. for two pulses coming from same digi
      const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses = iter->GetCrvRecoPulses();
      for(size_t i=0; i<crvRecoPulses.size(); i++)
      {
        const art::Ptr<CrvRecoPulse> crvRecoPulse = crvRecoPulses[i];
        art::Ptr<SimParticle> simParticle;
        double energyDeposited = 0;

        //get MC information, if available
        if(hasMCInfo)
        {
          std::set<art::Ptr<StepPointMC> > stepsThisPulse;
          const std::vector<size_t> &waveformIndices = crvRecoPulse->GetWaveformIndices();
          for(size_t j=0; j<waveformIndices.size(); j++) 
          {
            size_t waveformIndex = waveformIndices[j];
            CollectStepPoints(stepsThisPulse, crvDigiMCCollection, waveformIndex);
            CollectStepPoints(stepsAllPulses, crvDigiMCCollection, waveformIndex);
          }

          double earliestHitTimeThisPulse; //not used here
          CLHEP::Hep3Vector earliestHitPosThisPulse; //not used here
          ScanStepPoints(stepsThisPulse, energyDeposited, earliestHitTimeThisPulse, earliestHitPosThisPulse, simParticle);
        }

        pulses.emplace_back(simParticle,energyDeposited);
      }//loop over reco pulses

      ScanStepPoints(stepsAllPulses, totalEnergyDeposited, earliestHitTime, earliestHitPos, mostLikelySimParticle);

      //insert the cluster information into the vector of the crv coincidence clusters
      crvCoincidenceClusterMCCollection->emplace_back(hasMCInfo, pulses, mostLikelySimParticle, totalEnergyDeposited, earliestHitTime, earliestHitPos);
    }//loop over all clusters

    event.put(std::move(crvCoincidenceClusterMCCollection));
  } // end produce

  void CrvCoincidenceClusterMatchMC::CollectStepPoints(std::set<art::Ptr<StepPointMC> > &steps, 
                                                  const art::Handle<CrvDigiMCCollection> &digis, 
                                                  size_t waveformIndex)
  {
    const CrvDigiMC &digi = digis->at(waveformIndex);
    const std::vector<art::Ptr<StepPointMC> > &stepPoints = digi.GetStepPoints();
    for(size_t k=0; k<stepPoints.size(); k++)
    {
      if(stepPoints[k].isNonnull()) steps.insert(stepPoints[k]);
    }
  }

  void CrvCoincidenceClusterMatchMC::ScanStepPoints(const std::set<art::Ptr<StepPointMC> > &steps, 
                                               double &energyDeposited, double &earliestHitTime,
                                               CLHEP::Hep3Vector &earliestHitPos,
                                               art::Ptr<SimParticle> &mostLikelySimParticle)
  {
    energyDeposited=0;
    std::map<art::Ptr<SimParticle>,double> simParticleMap;
    std::set<art::Ptr<StepPointMC> >::const_iterator i;
    for(i=steps.begin(); i!=steps.end(); i++)
    {
      const StepPointMC &step = **i;
      energyDeposited+=step.totalEDep();
      simParticleMap[step.simParticle()]+=step.totalEDep();
    }

    std::map<art::Ptr<SimParticle>,double>::iterator simParticleIter;
    double simParticleDepEnergy=0;
    for(simParticleIter=simParticleMap.begin(); simParticleIter!=simParticleMap.end(); simParticleIter++)
    {
      if(simParticleIter->second>simParticleDepEnergy)
      {
        simParticleDepEnergy=simParticleIter->second;
        mostLikelySimParticle=simParticleIter->first;
      }
    }

    //time folding is not applied here, but was used to create the digis, ...
    //so we need to avoid that some step points from a different micro bunch 
    //could be accidentally found to be the step point with the earliest hit time.
    //therefore, only step points of the most likely sim particle will be considered.
    earliestHitTime=NAN;
    for(i=steps.begin(); i!=steps.end(); i++)
    {
      const StepPointMC &step = **i;
      if(step.simParticle()==mostLikelySimParticle)
      {
        double t = _timeOffsets.timeWithOffsetsApplied(step);
        if(isnan(earliestHitTime) || earliestHitTime>t)
        {
          earliestHitTime=t;
          earliestHitPos=step.position();
        }
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvCoincidenceClusterMatchMC;
DEFINE_ART_MODULE(CrvCoincidenceClusterMatchMC)
