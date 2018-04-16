//
// A module to find clusters of coincidences of CRV pulses
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
#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusters.hh"

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
  class CrvCoincidenceClusterFinder : public art::EDProducer 
  {
    public:
    explicit CrvCoincidenceClusterFinder(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    int         _verboseLevel;
    double      _maxDistance;
    double      _maxTimeDifference;
    std::string _crvCoincidenceCheckModuleLabel;  //module label of the CrvCoincidenceCheck module
                                                  //it is possible to have more than one instance of the CrvCoincidenceCheck module
    std::string _crvWaveformsModuleLabel;  //module label of the CrvWaveform module. 
                                           //this is optional. only needed, if MC information is required
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _deadTimeWindowStartMargin;
    double      _deadTimeWindowEndMargin;
    double      _totalDeadTime;
    double      _totalTime;

    struct ClusterHit
    {
      art::Ptr<CrvRecoPulse> _crvRecoPulse; 
      double _pos;
      double _time;
      CLHEP::Hep3Vector _counterPosition;
      ClusterHit(const art::Ptr<CrvRecoPulse> crvRecoPulse, double pos, double time, const CLHEP::Hep3Vector &counterPosition) 
            : _crvRecoPulse(crvRecoPulse), _pos(pos), _time(time), _counterPosition(counterPosition) {}
    };

    struct PosCompare
    {
      bool operator() (const ClusterHit &lhs, const ClusterHit &rhs) const 
      {                                                //identical hits will be removed, 
        if(lhs._pos!=rhs._pos) return lhs._pos < rhs._pos; 
        return lhs._crvRecoPulse < rhs._crvRecoPulse;
      }
    };

    struct TimeCompare
    {
      bool operator() (const ClusterHit &lhs, const ClusterHit &rhs) const 
      {                                                //identical hits will be removed, 
        if(lhs._time!=rhs._time) return lhs._time < rhs._time; 
        return lhs._crvRecoPulse < rhs._crvRecoPulse;
      }
    };

  };

  CrvCoincidenceClusterFinder::CrvCoincidenceClusterFinder(fhicl::ParameterSet const& pset) :
    _verboseLevel(pset.get<int>("verboseLevel")),
    _maxDistance(pset.get<double>("maxDistance")),
    _maxTimeDifference(pset.get<double>("maxTimeDifference")),
    _crvCoincidenceCheckModuleLabel(pset.get<std::string>("crvCoincidenceCheckModuleLabel")),
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel","")),
    _timeWindowStart(pset.get<double>("timeWindowStart")),
    _timeWindowEnd(pset.get<double>("timeWindowEnd")),
    _deadTimeWindowStartMargin(pset.get<double>("deadTimeWindowStartMargin")),
    _deadTimeWindowEndMargin(pset.get<double>("deadTimeWindowEndMargin")),
    _totalDeadTime(0),
    _totalTime(0)
  {
    produces<CrvCoincidenceClusters>();
  }

  void CrvCoincidenceClusterFinder::beginJob()
  {
  }

  void CrvCoincidenceClusterFinder::endJob()
  {
    double fractionDeadTime = _totalDeadTime / _totalTime;
    std::cout << "SUMMARY " << _crvCoincidenceCheckModuleLabel << "    Total dead time: " << _totalDeadTime << " / " << _totalTime << " = " << fractionDeadTime*100 << "%   using time window " << _timeWindowStart << " ns ... " << _timeWindowEnd << " ns" << std::endl;
  }

  void CrvCoincidenceClusterFinder::beginRun(art::Run &run)
  {
  }

  void CrvCoincidenceClusterFinder::produce(art::Event& event) 
  {
    _totalTime += _timeWindowEnd - _timeWindowStart;

    std::unique_ptr<CrvCoincidenceClusters> crvCoincidenceClusters(new CrvCoincidenceClusters);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvCoincidenceCheckResult> crvCoincidenceCheckResult;
    event.getByLabel(_crvCoincidenceCheckModuleLabel,"",crvCoincidenceCheckResult);

    if(crvCoincidenceCheckResult.product()==NULL) return;

    art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
    if(_crvWaveformsModuleLabel!="") event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection); //this is an optional part for MC information

    //loop through all coincidence combinations
    //extract all the hits of all coincidence combinations
    //distribute them into the crv sector types
    //and put them into a (pos ordered) set
    std::map<int, std::set<ClusterHit,PosCompare> > sectorTypeMap;

    const std::vector<CrvCoincidenceCheckResult::CoincidenceCombination> &coincidenceCombinationsAll = crvCoincidenceCheckResult->GetCoincidenceCombinations();
    std::vector<CrvCoincidenceCheckResult::CoincidenceCombination>::const_iterator iter;
    for(iter=coincidenceCombinationsAll.begin(); iter!=coincidenceCombinationsAll.end(); iter++)
    {
      int crvSectorType = iter->_sectorType;
      const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses = iter->_crvRecoPulses;
      for(size_t i=0; i<crvRecoPulses.size(); i++)
      {
        const art::Ptr<CrvRecoPulse> crvRecoPulse = crvRecoPulses[i];
 
        const mu2e::CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse->GetScintillatorBarIndex(); 
        const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
        int widthDirection=crvCounter.getBarDetail().getWidthDirection();
        double pos=crvCounter.getPosition()[widthDirection];

        sectorTypeMap[crvSectorType].emplace(crvRecoPulse, pos, crvRecoPulse->GetPulseTime(), crvCounter.getPosition());
      }
    }

    //loop through all crv sectors types
    std::map<int, std::set<ClusterHit,PosCompare> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); sectorTypeMapIter++)
    {
      int crvSectorType = sectorTypeMapIter->first;
      const std::set<ClusterHit,PosCompare> &posOrderedHits = sectorTypeMapIter->second; 

      //loop through the vector of (position ordered) crv hits for this particular crv sector type
      std::set<ClusterHit,PosCompare>::const_iterator hp=posOrderedHits.begin();
      while(hp!=posOrderedHits.end())
      {
        double startPos=hp->_pos;
        std::set<ClusterHit,TimeCompare> posCluster;  //time ordered hits of this position cluster
        posCluster.insert(*hp);
        while(++hp!=posOrderedHits.end())
        {
          if(hp->_pos-startPos>_maxDistance) break;
          posCluster.insert(*hp);
        }
        //filled one full position cluster

        //loop through the vector of (time ordered) crv hits for this position cluster
        std::set<ClusterHit,TimeCompare>::const_iterator ht=posCluster.begin();
        while(ht!=posCluster.end())
        {
          double startTime=ht->_time;
          std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses;
          crvRecoPulses.push_back(ht->_crvRecoPulse);

          double endTime=ht->_time;
          int PEs=ht->_crvRecoPulse->GetPEs();
          CLHEP::Hep3Vector avgCounterPos=ht->_counterPosition;

          while(++ht!=posCluster.end())
          {
            if(ht->_time-startTime>_maxTimeDifference) break;
            crvRecoPulses.push_back(ht->_crvRecoPulse);

            endTime=ht->_time;
            PEs+=ht->_crvRecoPulse->GetPEs();
            avgCounterPos+=ht->_counterPosition;
          }
          //filled one full position&time cluster

          //average counter position
          avgCounterPos/=crvRecoPulses.size();

          //get MC information, if available
          std::vector<art::Ptr<SimParticle> > simParticles;
          double energyDeposited = 0;
          double earliestHitTime = NAN;
          CLHEP::Hep3Vector earliestHitPos;
          if(crvDigiMCCollection.isValid())
          {
            //loop through all CrvRecoPulses to find all all CrvDigiMCs belongin to this cluster 
            std::map<art::Ptr<SimParticle>, int> simParticleMap;  //counting all simParticles
            for(size_t i=0; i<crvRecoPulses.size(); i++)
            {
              const std::vector<size_t> &waveformIndices = crvRecoPulses[i]->GetWaveformIndices();
              for(size_t j=0; j<waveformIndices.size(); j++)
              {
                size_t waveformIndex = waveformIndices[j];
                const CrvDigiMC &digi = crvDigiMCCollection->at(waveformIndex);

                //get the SimParticle responsible for this CrvDigiMC
                const art::Ptr<SimParticle> simParticle = digi.GetSimParticle();
                //there is a special case in which there is no simParticle:
                //a reco pulse based on two consecutive single waveforms, 
                //where the second waveform is either caused by noise, or just the tail end of the first waveform
                if(simParticle.isNonnull()) simParticleMap[simParticle]++;

                //get the StepPointMCs responsible for this CrvDigiMC
                const std::vector<art::Ptr<StepPointMC> > &stepPoints = digi.GetStepPoints();
                for(size_t k=0; k<stepPoints.size(); k++)
                {
                  if(!stepPoints[k].isNonnull()) continue;
                  energyDeposited = stepPoints[k]->totalEDep();
                  double t=stepPoints[k]->time();
                  if(t<earliestHitTime || isnan(earliestHitTime))
                  {
                    earliestHitTime=t;
                    earliestHitPos=stepPoints[k]->position();
                  }
                }
              }
            }
            //finding the three most frequent simParticles
            for(int i=0; i<3; i++)
            {
              art::Ptr<SimParticle> simParticle;
              std::map<art::Ptr<SimParticle>,int >::iterator simParticleIter;
              int simParticleCount=0;
              for(simParticleIter=simParticleMap.begin(); simParticleIter!=simParticleMap.end(); simParticleIter++)
              {
                if(std::find(simParticles.begin(),simParticles.end(),simParticleIter->first)!=simParticles.end()) continue; //ignore particles which have been found already
                if(simParticleIter->second>simParticleCount)
                {
                  simParticleCount=simParticleIter->second;
                  simParticle=simParticleIter->first;
                }
              }
              if(simParticle.isNonnull()) simParticles.push_back(simParticle);  //need to check whether it actually found e.g. a 2nd or 3rd most frequent particle
            }


          }//MC information

          //insert the cluster information into the vector of the crv coincidence clusters
          crvCoincidenceClusters->GetClusters().emplace_back(crvSectorType, avgCounterPos, startTime, endTime, PEs, 
                                                             crvRecoPulses, simParticles, energyDeposited, earliestHitTime, earliestHitPos);

          if(_verboseLevel>0)
          {
            double deadTimeWindowStart = startTime-_deadTimeWindowStartMargin;
            double deadTimeWindowEnd   = endTime+_deadTimeWindowEndMargin;
            if(deadTimeWindowStart<_timeWindowStart) deadTimeWindowStart=_timeWindowStart;
            if(deadTimeWindowEnd<_timeWindowStart) continue;
            if(deadTimeWindowStart>_timeWindowEnd) continue;
            if(deadTimeWindowEnd>_timeWindowEnd) deadTimeWindowEnd=_timeWindowEnd;
   
            double deadTime = deadTimeWindowEnd-deadTimeWindowStart;
            std::cout << "   Found dead time window: " << deadTimeWindowStart << " ns ... " << deadTimeWindowEnd << " ns   (dead time incl. start/end margins: "<<deadTime<<" ns)";
            std::cout << "   in CRV region " << crvSectorType << std::endl;
            _totalDeadTime += deadTime;
            double fractionDeadTime = _totalDeadTime / _totalTime;
            std::cout << "Dead time so far: " << _totalDeadTime << " ns / " << _totalTime << " ns = " << fractionDeadTime*100 << "%    using time window " << _timeWindowStart << " ns ... " << _timeWindowEnd << " ns" << std::endl;
          }


        }//loop over all hits in position cluster
      }//loop over all hits in sector
    }//loop over all sector types

    event.put(std::move(crvCoincidenceClusters));
  } // end produce

} // end namespace mu2e

using mu2e::CrvCoincidenceClusterFinder;
DEFINE_ART_MODULE(CrvCoincidenceClusterFinder)
