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
    std::string _crvRecoPulseModuleLabel; 
    std::string _crvCoincidenceCheckModuleLabel;  //module label of the CrvCoincidenceCheck module
                                                  //it is possible to have more than one instance of the CrvCoincidenceCheck module
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _deadTimeWindowStartMargin;
    double      _deadTimeWindowEndMargin;
    double      _totalDeadTime;
    double      _totalTime;

    struct PosOrderedHit : public CrvCoincidenceClusters::Hit
    {
      PosOrderedHit(double time, int PEs, mu2e::CRSScintillatorBarIndex counter, int SiPM, double pos): Hit(time, PEs, counter, SiPM, pos) {}
      bool operator<(const PosOrderedHit &rhs) const 
      {                                                //identical hits will be removed, 
        if(_pos!=rhs._pos) return _pos < rhs._pos;     //but allow for hits with the same pos but on different counters/SiPMs/times
        if(_counter!=rhs._counter) return _counter < rhs._counter;
        if(_SiPM!=rhs._SiPM) return _SiPM < rhs._SiPM;
        return _time < rhs._time;
      }
    };

    struct TimeOrderedHit : public CrvCoincidenceClusters::Hit
    {
      TimeOrderedHit(double time, int PEs, mu2e::CRSScintillatorBarIndex counter, int SiPM, double pos): Hit(time, PEs, counter, SiPM, pos) {}
      bool operator<(const TimeOrderedHit &rhs) const  
      {                                                //identical hits will be removed, 
        if(_time!=rhs._time) return _time < rhs._time; //but allow for hits with the same time but on different counters/SiPMs
        if(_counter!=rhs._counter) return _counter < rhs._counter;
        return _SiPM < rhs._SiPM;
      }
    };

  };

  CrvCoincidenceClusterFinder::CrvCoincidenceClusterFinder(fhicl::ParameterSet const& pset) :
    _verboseLevel(pset.get<int>("verboseLevel")),
    _maxDistance(pset.get<double>("maxDistance")),
    _maxTimeDifference(pset.get<double>("maxTimeDifference")),
    _crvCoincidenceCheckModuleLabel(pset.get<std::string>("crvCoincidenceCheckModuleLabel")),
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

    std::string crvCoincidenceInstanceName="";
    art::Handle<CrvCoincidenceCheckResult> crvCoincidenceCheckResult;
    event.getByLabel(_crvCoincidenceCheckModuleLabel,crvCoincidenceInstanceName,crvCoincidenceCheckResult);

    if(crvCoincidenceCheckResult.product()==NULL) return;

    //loop through all coincidence combinations
    //extract all the hits of all coincidence combinations
    //distribute them into the crv sector types
    //and put them into a (pos ordered) set
    std::map<int, std::set<PosOrderedHit> > sectorTypeMap;

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

        int widthDirection=crvCounter.getBarDetail().getWidthDirection();
        double pos=crvCounter.getPosition()[widthDirection];

        sectorTypeMap[crvSectorType].emplace(iter->_time[i], iter->_PEs[i], iter->_counters[i], iter->_SiPMs[i], pos);
      }
    }

    //loop through all crv sectors types
    std::map<int, std::set<PosOrderedHit> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); sectorTypeMapIter++)
    {
      int crvSectorType = sectorTypeMapIter->first;
      const std::set<PosOrderedHit> &crvHits = sectorTypeMapIter->second; 

      //loop through the vector of (position ordered) crv hits for this particular crv sector type
      std::set<PosOrderedHit>::const_iterator h=crvHits.begin();
      while(h!=crvHits.end())
      {
        double endPos=h->_pos;
        std::set<TimeOrderedHit> posCluster;
        posCluster.emplace(h->_time,h->_PEs,h->_counter,h->_SiPM,h->_pos);
        while(++h!=crvHits.end())
        {
          if(h->_pos-endPos>_maxDistance) break;
          endPos=h->_pos;
          posCluster.emplace(h->_time,h->_PEs,h->_counter,h->_SiPM,h->_pos);
        }
        //filled one full position cluster

        //loop through the vector of (time ordered) crv hits for this cluster
        std::set<TimeOrderedHit>::const_iterator hh=posCluster.begin();
        while(hh!=posCluster.end())
        {
          double startTime=hh->_time;
          double endTime=hh->_time;
          std::vector<const CrvCoincidenceClusters::Hit*> posTimeCluster;
          posTimeCluster.push_back(&(*hh));

          int PEs=hh->_PEs;
          CLHEP::Hep3Vector avgPos=CRS->getBar(hh->_counter).getPosition();

          while(++hh!=posCluster.end())
          {
            if(hh->_time-endTime>_maxTimeDifference) break;
            endTime=hh->_time;
            posTimeCluster.push_back(&(*hh));

            PEs+=hh->_PEs;
            avgPos+=CRS->getBar(hh->_counter).getPosition();
          }
          //filled one full position&time cluster

          //average position
          avgPos/=posTimeCluster.size();

          //insert the cluster information into the vector of the crv coincidence clusters
          crvCoincidenceClusters->GetClusters().emplace_back(crvSectorType, avgPos, startTime, endTime, PEs, posTimeCluster);

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
