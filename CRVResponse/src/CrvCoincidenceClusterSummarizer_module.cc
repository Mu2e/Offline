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
#include "RecoDataProducts/inc/CrvRecoPulse.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterSummaryCollection.hh"

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
  class CrvCoincidenceClusterSummarizer : public art::EDProducer 
  {
    public:
    explicit CrvCoincidenceClusterSummarizer(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    std::string _crvCoincidenceClusterFinderModuleLabel;  //module label of the CrvCoincidenceClusterFinder module
                                                          //it is possible to have more than one instance of the CrvCoincidenceClusterFinder module
    std::string _crvWaveformsModuleLabel;  //module label of the CrvWaveform module. 
                                           //this is optional. only needed, if MC information is required
  };

  CrvCoincidenceClusterSummarizer::CrvCoincidenceClusterSummarizer(fhicl::ParameterSet const& pset) :
    _crvCoincidenceClusterFinderModuleLabel(pset.get<std::string>("crvCoincidenceClusterFinderModuleLabel")),
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel",""))
  {
    produces<CrvCoincidenceClusterSummaryCollection>();
  }

  void CrvCoincidenceClusterSummarizer::beginJob()
  {
  }

  void CrvCoincidenceClusterSummarizer::endJob()
  {
  }

  void CrvCoincidenceClusterSummarizer::beginRun(art::Run &run)
  {
  }

  void CrvCoincidenceClusterSummarizer::produce(art::Event& event) 
  {
    std::unique_ptr<CrvCoincidenceClusterSummaryCollection> crvCoincidenceClusterSummaryCollection(new CrvCoincidenceClusterSummaryCollection);

    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;
    event.getByLabel(_crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusterCollection);

    if(crvCoincidenceClusterCollection.product()==NULL) return;

    art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
    if(_crvWaveformsModuleLabel!="") event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection); //this is an optional part for MC information

    CrvCoincidenceClusterCollection::const_iterator iter;
    for(iter=crvCoincidenceClusterCollection->begin(); iter!=crvCoincidenceClusterCollection->end(); iter++)
    {
      int    crvSectorType            = iter->GetCrvSectorType();
      CLHEP::Hep3Vector avgCounterPos = iter->GetAvgCounterPos();
      double startTime                = iter->GetStartTime();
      double endTime                  = iter->GetEndTime();
      int    PEs                      = iter->GetPEs();
      bool   hasMCInfo                = (crvDigiMCCollection.isValid()?true:false); //MC
      double totalEnergyDeposited     = 0;         //MC
      double earliestHitTime          = NAN;       //MC
      art::Ptr<SimParticle> mostLikelySimParticle; //MC
      CLHEP::Hep3Vector     earliestHitPos;        //MC

      //loop through all reco pulses and try to find the MC information
      std::vector<CrvCoincidenceClusterSummary::PulseInfo> pulses; //collection all pulses (incl. reco and MC info)
      std::map<art::Ptr<SimParticle>, int> simParticleMap;  //counting all simParticles
      std::set<size_t> allWaveformIndices;  //collecting all waveform indices for total energy calculation
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
          const std::vector<size_t> &waveformIndices = crvRecoPulse->GetWaveformIndices();
          double maxEnergyDeposited = 0;
          for(size_t j=0; j<waveformIndices.size(); j++)
          {
            size_t waveformIndex = waveformIndices[j];
            allWaveformIndices.insert(waveformIndex);

            const CrvDigiMC &digi = crvDigiMCCollection->at(waveformIndex);

            //get the SimParticle responsible for the current pulse from this CrvDigiMC,
            //if it has makes the biggest contribution with respect to the deposited energy
            if(digi.GetEnergyDeposited()>maxEnergyDeposited)
            {
              maxEnergyDeposited = digi.GetEnergyDeposited();
              simParticle = digi.GetSimParticle();
            }

            //add the deposited energy for the current pulse from this CrvDigiMC
            energyDeposited += digi.GetEnergyDeposited();
          }
        }

        pulses.emplace_back(crvRecoPulse,simParticle,energyDeposited);
      }//loop over reco pulses

      //looping through all crvDigiMCs (if available)
      for(std::set<size_t>::const_iterator i=allWaveformIndices.begin(); i!=allWaveformIndices.end(); i++)
      {
        size_t waveformIndex = *i;
        const CrvDigiMC &digi = crvDigiMCCollection->at(waveformIndex);

        //add the deposited energy for the current CrvDigiMC
        totalEnergyDeposited += digi.GetEnergyDeposited();

        //put the SimParticle of this CrvDigiMC into the collection of SimParticles
        simParticleMap[digi.GetSimParticle()]++;
          
        //get the StepPointMCs responsible for this CrvDigiMC
        const std::vector<art::Ptr<StepPointMC> > &stepPoints = digi.GetStepPoints();
        for(size_t k=0; k<stepPoints.size(); k++)
        {
          if(stepPoints[k].isNull()) continue;
          double t=stepPoints[k]->time();
          if(t<earliestHitTime || isnan(earliestHitTime))
          {
            earliestHitTime=t;
            earliestHitPos=stepPoints[k]->position();
          }
        }
      }//loop over crvDigiMCs

      //finding the most likely simParticles
      std::map<art::Ptr<SimParticle>,int >::iterator simParticleIter;
      int simParticleCount=0;
      for(simParticleIter=simParticleMap.begin(); simParticleIter!=simParticleMap.end(); simParticleIter++)
      {
        if(simParticleIter->second>simParticleCount)
        {
          simParticleCount=simParticleIter->second;
          mostLikelySimParticle=simParticleIter->first;
        }
      }//loop over sim particle map

      //insert the cluster information into the vector of the crv coincidence clusters
      crvCoincidenceClusterSummaryCollection->emplace_back(crvSectorType, avgCounterPos, startTime, endTime, PEs, hasMCInfo, 
                                                           pulses, mostLikelySimParticle, totalEnergyDeposited, earliestHitTime, earliestHitPos);
    }//loop over all clusters

    event.put(std::move(crvCoincidenceClusterSummaryCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvCoincidenceClusterSummarizer;
DEFINE_ART_MODULE(CrvCoincidenceClusterSummarizer)
