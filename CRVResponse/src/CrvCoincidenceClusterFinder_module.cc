//
// A module to find clusters of coincidences of CRV pulses
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidence.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
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
    bool        _usingPEsPulseHeight;
    double      _maxDistance;
    double      _maxTimeDifference;
    std::string _crvCoincidenceCheckModuleLabel;  //module label of the CrvCoincidenceCheck module
                                                  //it is possible to have more than one instance of the CrvCoincidenceCheck module
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _deadTimeWindowStartMargin;
    double      _deadTimeWindowEndMargin;
    double      _totalDeadTime;
    double      _totalTime;

    struct ClusterHit
    {
      art::Ptr<CrvRecoPulse> _crvRecoPulse;
      float  _x, _y;
      double _time;
      CLHEP::Hep3Vector _counterPosition;
      int    _layer;
      ClusterHit(const art::Ptr<CrvRecoPulse> crvRecoPulse, float x, float y, double time, const CLHEP::Hep3Vector &counterPosition, int layer)
            : _crvRecoPulse(crvRecoPulse), _x(x), _y(y), _time(time), _counterPosition(counterPosition), _layer(layer) {}
    };

    struct PosCompare
    {
      bool operator() (const ClusterHit &lhs, const ClusterHit &rhs) const
      {                                                //identical hits will be removed,
        if(lhs._x!=rhs._x) return lhs._x < rhs._x;
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
    art::EDProducer{pset},
    _verboseLevel(pset.get<int>("verboseLevel")),
    _usingPEsPulseHeight(pset.get<bool>("usingPEsPulseHeight")),
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
    produces<CrvCoincidenceClusterCollection>();
  }

  void CrvCoincidenceClusterFinder::beginJob()
  {
  }

  void CrvCoincidenceClusterFinder::endJob()
  {
    double fractionDeadTime = _totalDeadTime / _totalTime;
    if(_verboseLevel>0)
      {
        std::cout << "SUMMARY " << _crvCoincidenceCheckModuleLabel << "    Total dead time: " << _totalDeadTime << " / " << _totalTime << " = " << fractionDeadTime*100 << "%   using time window " << _timeWindowStart << " ns ... " << _timeWindowEnd << " ns" << std::endl;
      }
  }

  void CrvCoincidenceClusterFinder::beginRun(art::Run &run)
  {
  }

  void CrvCoincidenceClusterFinder::produce(art::Event& event)
  {
    _totalTime += _timeWindowEnd - _timeWindowStart;

    std::unique_ptr<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection(new CrvCoincidenceClusterCollection);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvCoincidenceCollection> crvCoincidenceCollection;
    event.getByLabel(_crvCoincidenceCheckModuleLabel,"",crvCoincidenceCollection);

    if(crvCoincidenceCollection.product()==NULL) return;

    //loop through all coincidence combinations
    //extract all the hits of all coincidence combinations
    //distribute them into the crv sector types
    //and put them into a (x-pos ordered) set
    std::map<int, std::set<ClusterHit,PosCompare> > sectorTypeMap;

    CrvCoincidenceCollection::const_iterator iter;
    for(iter=crvCoincidenceCollection->begin(); iter!=crvCoincidenceCollection->end(); iter++)
    {
      int crvSectorType = std::abs(iter->GetCrvSectorType());
      const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses = iter->GetCrvRecoPulses();
      for(size_t i=0; i<crvRecoPulses.size(); i++)
      {
        const art::Ptr<CrvRecoPulse> crvRecoPulse = crvRecoPulses[i];

        const mu2e::CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse->GetScintillatorBarIndex();
        const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
        int widthDirection=crvCounter.getBarDetail().getWidthDirection();
        int thicknessDirection=crvCounter.getBarDetail().getThicknessDirection();
        float x=crvCounter.getPosition()[widthDirection];
        float y=crvCounter.getPosition()[thicknessDirection];
        int layer=crvCounter.id().getLayerNumber();
        sectorTypeMap[crvSectorType].emplace(crvRecoPulse, x, y, crvRecoPulse->GetPulseTime(), crvCounter.getPosition(), layer);
      }
    }

    //loop through all crv sectors types
    std::map<int, std::set<ClusterHit,PosCompare> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); sectorTypeMapIter++)
    {
      int crvSectorType = sectorTypeMapIter->first;
      const std::set<ClusterHit,PosCompare> &posOrderedHits = sectorTypeMapIter->second;

      //loop through the vector of (x-position ordered) crv hits for this particular crv sector type
      std::set<ClusterHit,PosCompare>::const_iterator hp=posOrderedHits.begin();
      while(hp!=posOrderedHits.end())
      {
        double startPos=hp->_x;
        std::set<ClusterHit,TimeCompare> posCluster;  //time ordered hits of this x-position cluster
        posCluster.insert(*hp);
        while(++hp!=posOrderedHits.end())
        {
          if(hp->_x-startPos>_maxDistance) break;
          posCluster.insert(*hp);
        }
        //filled one full x-position cluster

        //loop through the vector of (time ordered) crv hits for this x-position cluster
        std::set<ClusterHit,TimeCompare>::const_iterator ht=posCluster.begin();
        while(ht!=posCluster.end())
        {
          double startTime=ht->_time;
          std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses;
          crvRecoPulses.push_back(ht->_crvRecoPulse);

          double endTime=ht->_time;
          float PEs=ht->_crvRecoPulse->GetPEs();
          if(_usingPEsPulseHeight) PEs=ht->_crvRecoPulse->GetPEsPulseHeight();
          CLHEP::Hep3Vector avgCounterPos=ht->_counterPosition*PEs;  //PE-weighted average position

          std::set<int> layerSet;
          layerSet.insert(ht->_layer);
          int nPoints=1;
          float sumX =ht->_x;
          float sumY =ht->_y;
          float sumYY=ht->_y*ht->_y;
          float sumXY=ht->_x*ht->_y;
          while(++ht!=posCluster.end())
          {
            if(ht->_time-startTime>_maxTimeDifference) break;
            crvRecoPulses.push_back(ht->_crvRecoPulse);

            endTime=ht->_time;
            float PEsTmp=ht->_crvRecoPulse->GetPEs();
            if(_usingPEsPulseHeight) PEsTmp=ht->_crvRecoPulse->GetPEsPulseHeight();
            PEs+=PEsTmp;
            avgCounterPos+=ht->_counterPosition*PEsTmp;

            layerSet.insert(ht->_layer);
            ++nPoints;
            sumX +=ht->_x;
            sumY +=ht->_y;
            sumYY+=ht->_y*ht->_y;
            sumXY+=ht->_x*ht->_y;
          }
          //filled one full position&time cluster
          if(layerSet.size()<2) continue;

          //average counter position (PE weighted), slope, layers
          avgCounterPos/=PEs;
          float slope=(nPoints*sumXY-sumX*sumY)/(nPoints*sumYY-sumY*sumY);
          std::vector<int> layers(layerSet.begin(), layerSet.end());

          //insert the cluster information into the vector of the crv coincidence clusters
          crvCoincidenceClusterCollection->emplace_back(crvSectorType, avgCounterPos, startTime, endTime, PEs, crvRecoPulses, slope, layers);

          //calculate dead time
          double deadTimeWindowStart = startTime-_deadTimeWindowStartMargin;
          double deadTimeWindowEnd   = endTime+_deadTimeWindowEndMargin;
          if(deadTimeWindowStart<_timeWindowStart) deadTimeWindowStart=_timeWindowStart;
          if(deadTimeWindowEnd<_timeWindowStart) continue;
          if(deadTimeWindowStart>_timeWindowEnd) continue;
          if(deadTimeWindowEnd>_timeWindowEnd) deadTimeWindowEnd=_timeWindowEnd;

          double deadTime = deadTimeWindowEnd-deadTimeWindowStart;
          _totalDeadTime += deadTime;

          if(_verboseLevel>0)
          {
            std::cout << "   Found dead time window: " << deadTimeWindowStart << " ns ... " << deadTimeWindowEnd << " ns   (dead time incl. start/end margins: "<<deadTime<<" ns)";
            std::cout << "   in CRV region " << crvSectorType << std::endl;
            double fractionDeadTime = _totalDeadTime / _totalTime;
            std::cout << "Dead time so far: " << _totalDeadTime << " ns / " << _totalTime << " ns = " << fractionDeadTime*100 << "%    using time window " << _timeWindowStart << " ns ... " << _timeWindowEnd << " ns" << std::endl;
          }


        }//loop over all hits in position cluster
      }//loop over all hits in sector
    }//loop over all sector types

    event.put(std::move(crvCoincidenceClusterCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvCoincidenceClusterFinder;
DEFINE_ART_MODULE(CrvCoincidenceClusterFinder)
