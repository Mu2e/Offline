//
// A module to find clusters of coincidences of CRV pulses
//
// Original Author: Ralf Ehrlich

#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVReco/inc/CrvHelper.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

#include <string>
#include <array>

namespace mu2e
{
  class CrvCoincidenceFinder : public art::EDProducer
  {
    public:
    struct SectorConfig
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> CRVSector{Name("CRVSector"), Comment("CRV sector")};
      fhicl::Atom<int> PEthreshold{Name("PEthreshold"), Comment("PE threshold required for a coincidence")};
      fhicl::Atom<double> maxTimeDifferenceAdjacentPulses{Name("maxTimeDifferenceAdjacentPulses"), Comment("maximum time difference of pulses of adjacent channels considered for coincidences")};
      fhicl::Atom<double> maxTimeDifference{Name("maxTimeDifference"), Comment("maximum time difference of a coincidence hit combination")};
      fhicl::Atom<double> minSlope{Name("minSlope"), Comment("minimum slope allowed for a coincidence")};
      fhicl::Atom<double> maxSlope{Name("maxSlope"), Comment("maximum slope allowed for a coincidence")};
      fhicl::Atom<double> maxSlopeDifference{Name("maxSlopeDifference"), Comment("maximum slope difference between layers allowed for a coincidence")};
      fhicl::Atom<int> coincidenceLayers{Name("coincidenceLayers"), Comment("number of layers required for a coincidence")};
      fhicl::Atom<double> minClusterPEs{Name("minClusterPEs"), Comment("minimum number of PEs in a cluster required for storage")};
      fhicl::Atom<double> initialClusterMaxDistance{Name("initialClusterMaxDistance"), Comment("maximum distances between hits to be considered for a hit cluster at initial clustering process")};
    };
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> verboseLevel{Name("verboseLevel"), Comment("verbose level")};
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("module label of the input CrvRecoPulses")};
      //sector-specific coincidence settings
      fhicl::Sequence<fhicl::Table<SectorConfig> > sectorConfig{Name("sectorConfig"), Comment("sector-specific settings")};
      //other settings
      fhicl::Atom<int> bigClusterThreshold{Name("bigClusterThreshold"), Comment("no coincidence check for clusters with a number of hits above this threshold")};
      fhicl::Atom<double> fiberSignalSpeed{Name("fiberSignalSpeed"), Comment("effective speed of signals inside the CRV fibers in mm/ns")};
      fhicl::Atom<double> maxTimeInterval{Name("maxTimeInterval"), Comment("maximum time interval of hits on one readout side considered for hit position calculation")};
      fhicl::Atom<double> timeOffset{Name("timeOffset"), Comment("additional time delay caused by electronics response and physical processes in ns")};
      fhicl::Sequence<int> compensateChannelStatus{Name("compensateChannelStatus"), Comment("compensate missinge pulses for channels with the following channel statuses")};
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit CrvCoincidenceFinder(const Parameters& config);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    int                      _verboseLevel;
    std::string              _crvRecoPulsesModuleLabel;

    std::vector<SectorConfig> _sectorConfig;

    size_t      _bigClusterThreshold;
    double      _fiberSignalSpeed;
    double      _maxTimeInterval;
    double      _timeOffset;

    mu2e::ProditionsHandle<mu2e::CRVStatus> _sipmStatus;
    std::vector<int>                        _compensateChannelStatus;

    int         _totalEvents;
    int         _totalEventsCoincidence;

    struct sectorCoincidenceProperties
    {
      int  precedingCounters;
      int  nCountersPerModule;
      int  sectorType;
      bool sipmsAtSide0;
      bool sipmsAtSide1;
      int  widthDirection, thicknessDirection, lengthDirection;
      double      counterHalfLength;
      double      posSide0, posSide1;
      int         PEthreshold;
      double      maxTimeDifferenceAdjacentPulses;
      double      maxTimeDifference;
      double      minSlope, maxSlope, maxSlopeDifference;
      int         coincidenceLayers;
      double      minClusterPEs;
      double      initialClusterMaxDistance;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;

    struct sectorTypeProperties
    {
      double maxTimeDifference;
      double maxHalfLength;
      double posSide0;
      double posSide1;
      bool   sipmsAtSide0;
      bool   sipmsAtSide1;
    };
    std::map<int,sectorTypeProperties> _sectorTypeMap;

    struct CrvHit
    {
      art::Ptr<CrvRecoPulse> _crvRecoPulse;
      CLHEP::Hep3Vector      _pos;
      double _x, _y;
      double _time;
      double _PEs;
      int    _crvSector;
      int    _layer;
      int    _counter;
      int    _SiPM;
      int    _PEthreshold;
      double _maxTimeDifferenceAdjacentPulses;
      double _maxTimeDifference;
      double _minSlope, _maxSlope, _maxSlopeDifference;
      int    _coincidenceLayers;
      double _minClusterPEs;
      mutable double _maxDistance; //initially set to initialClusterMaxDistance, which is just an estimate
                                   //used for the initial clustering process (to keep the number of
                                   //hit combinations down that need to be checked for coincidendes).
                                   //_maxDistance is updated when the coindidences are checked
                                   //(used for the final clustering process)

      CrvHit(const art::Ptr<CrvRecoPulse> crvRecoPulse, const CLHEP::Hep3Vector &pos,
             double x, double y, double time, double PEs, int crvSector, int layer, int counter, int SiPM, int PEthreshold,
             double maxTimeDifferenceAdjacentPulses, double maxTimeDifference,
             double minSlope, double maxSlope, double maxSlopeDifference, int coincidenceLayers, double minClusterPEs, double maxDistance) :
               _crvRecoPulse(crvRecoPulse), _pos(pos),
               _x(x), _y(y), _time(time), _PEs(PEs), _crvSector(crvSector), _layer(layer), _counter(counter), _SiPM(SiPM), _PEthreshold(PEthreshold),
               _maxTimeDifferenceAdjacentPulses(maxTimeDifferenceAdjacentPulses), _maxTimeDifference(maxTimeDifference),
               _minSlope(minSlope), _maxSlope(maxSlope), _maxSlopeDifference(maxSlopeDifference), _coincidenceLayers(coincidenceLayers),
               _minClusterPEs(minClusterPEs), _maxDistance(maxDistance) {}
    };

    void clusterProperties(int crvSectorType, const std::vector<std::vector<CrvHit> > &clusters,
                           std::unique_ptr<CrvCoincidenceClusterCollection> &crvCoincidenceClusterCollection,
                           const art::Handle<CrvRecoPulseCollection> &crvRecoPulseCollection);
    void filterHits(const std::vector<CrvHit> &hits, std::list<CrvHit> &hitsFiltered, unsigned int readoutSide);
    void findClusters(std::list<CrvHit> &hits, std::vector<std::vector<CrvHit> > &clusters, double clusterMaxTimeDifference,
                      bool combineOppositeSides=false, double timeDifferenceOppositeSides=0);
    void checkCoincidence(const std::vector<CrvHit> &hits, std::list<CrvHit> &coincidenceHits);
    bool checkCombination(std::vector<CrvHit>::const_iterator layerIterators[], int n);

  };

  CrvCoincidenceFinder::CrvCoincidenceFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _verboseLevel(conf().verboseLevel()),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _sectorConfig(conf().sectorConfig()),
    _bigClusterThreshold(conf().bigClusterThreshold()),
    _fiberSignalSpeed(conf().fiberSignalSpeed()),
    _maxTimeInterval(conf().maxTimeInterval()),
    _timeOffset(conf().timeOffset()),
    _compensateChannelStatus(conf().compensateChannelStatus()),
    _totalEvents(0),
    _totalEventsCoincidence(0)
  {
    produces<CrvCoincidenceClusterCollection>();
  }

  void CrvCoincidenceFinder::beginJob()
  {
  }

  void CrvCoincidenceFinder::endJob()
  {
    if(_verboseLevel>0)
    {
      std::cout<<"SUMMARY "<<moduleDescription().moduleLabel()<<"    "<<_totalEventsCoincidence<<" / "<<_totalEvents<<" events satisfied coincidence requirements"<<std::endl;
    }
  }

  void CrvCoincidenceFinder::beginRun(art::Run &run)
  {
    //first get some properties of each CRV sector
    //-from the geometry
    //-from the fcl parameters
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<CRSScintillatorShield> &sectors = CRS->getCRSScintillatorShields();
    std::map<int,std::vector<sectorCoincidenceProperties> > sectorTypes;
    for(size_t i=0; i<sectors.size(); ++i)
    {
      sectorCoincidenceProperties s;
      s.precedingCounters=0;
      int precedingSector=sectors[i].getPrecedingSector();
      while(precedingSector!=-1)
      {
        int nModules=sectors[precedingSector].nModules();
        int nCountersPerModule=sectors[precedingSector].getCountersPerModule();
        s.precedingCounters+=nModules*nCountersPerModule;
        precedingSector=sectors[precedingSector].getPrecedingSector();
      }

      s.nCountersPerModule=sectors[i].getCountersPerModule();

      s.sectorType=sectors[i].getSectorType();

      s.sipmsAtSide0=sectors[i].getCRSScintillatorBarDetail().hasCMB(0);
      s.sipmsAtSide1=sectors[i].getCRSScintillatorBarDetail().hasCMB(1);

      s.widthDirection=sectors[i].getCRSScintillatorBarDetail().getWidthDirection();
      s.thicknessDirection=sectors[i].getCRSScintillatorBarDetail().getThicknessDirection();
      s.lengthDirection=sectors[i].getCRSScintillatorBarDetail().getLengthDirection();
      s.counterHalfLength=sectors[i].getCRSScintillatorBarDetail().getHalfLength();

      s.posSide0=sectors[i].getFirstBar().getPosition()[s.lengthDirection] - s.counterHalfLength;
      s.posSide1=sectors[i].getFirstBar().getPosition()[s.lengthDirection] + s.counterHalfLength;

      std::string sectorName = sectors[i].getName().substr(4); //removes the "CRV_" part
      std::vector<SectorConfig>::iterator sectorConfigIter=std::find_if(_sectorConfig.begin(), _sectorConfig.end(),
                                                                        [sectorName](const SectorConfig &s){return s.CRVSector()==sectorName;});
      if(sectorConfigIter==_sectorConfig.end())
        throw std::logic_error("CrvCoincidenceFinder: The geometry has a CRV sector for which no coincidence properties were defined in the fcl file.");

      s.PEthreshold                     = sectorConfigIter->PEthreshold();
      s.maxTimeDifferenceAdjacentPulses = sectorConfigIter->maxTimeDifferenceAdjacentPulses();
      s.maxTimeDifference               = sectorConfigIter->maxTimeDifference();
      s.minSlope                        = sectorConfigIter->minSlope();
      s.maxSlope                        = sectorConfigIter->maxSlope();
      s.maxSlopeDifference              = sectorConfigIter->maxSlopeDifference();
      s.coincidenceLayers               = sectorConfigIter->coincidenceLayers();
      s.minClusterPEs                   = sectorConfigIter->minClusterPEs();
      s.initialClusterMaxDistance       = sectorConfigIter->initialClusterMaxDistance();

      _sectorMap[i]=s;
      sectorTypes[s.sectorType].push_back(s);
    }

    for(const auto& [sectorType,sectorVector] : sectorTypes)
    {
      _sectorTypeMap[sectorType].posSide0=std::min_element(sectorVector.begin(),sectorVector.end(),
                                                           [](const sectorCoincidenceProperties &a, const sectorCoincidenceProperties &b)
                                                           {return a.posSide0<b.posSide0;})->posSide0;
      _sectorTypeMap[sectorType].posSide1=std::max_element(sectorVector.begin(),sectorVector.end(),
                                                           [](const sectorCoincidenceProperties &a, const sectorCoincidenceProperties &b)
                                                           {return a.posSide1<b.posSide1;})->posSide1;

      _sectorTypeMap[sectorType].maxHalfLength=(_sectorTypeMap[sectorType].posSide1-_sectorTypeMap[sectorType].posSide0)/2.0;

      _sectorTypeMap[sectorType].maxTimeDifference=std::max_element(sectorVector.begin(),sectorVector.end(),
                                                                    [](const sectorCoincidenceProperties &a, const sectorCoincidenceProperties &b)
                                                                    {return a.maxTimeDifference<b.maxTimeDifference;})->maxTimeDifference;

      _sectorTypeMap[sectorType].sipmsAtSide0=std::any_of(sectorVector.begin(),sectorVector.end(),
                                                          [](const sectorCoincidenceProperties &a){return a.sipmsAtSide0;});
      _sectorTypeMap[sectorType].sipmsAtSide1=std::any_of(sectorVector.begin(),sectorVector.end(),
                                                          [](const sectorCoincidenceProperties &a){return a.sipmsAtSide1;});
    }
  }

  void CrvCoincidenceFinder::produce(art::Event& event)
  {
    std::unique_ptr<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection(new CrvCoincidenceClusterCollection);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);

    if(crvRecoPulseCollection.product()==NULL) return;

    auto const& sipmStatus = _sipmStatus.get(event.id());

    //loop through all reco pulses
    //distribute them into the crv sector types
    std::map<int, std::vector<CrvHit> > sectorTypeMap;
    for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); ++recoPulseIndex)
    {
      const art::Ptr<CrvRecoPulse> crvRecoPulse(crvRecoPulseCollection, recoPulseIndex);

      //get information about the counter
      const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse->GetScintillatorBarIndex();
      int sectorNumber, moduleNumber, layerNumber, barNumber;
      CrvHelper::GetCrvCounterInfo(CRS, crvBarIndex, sectorNumber, moduleNumber, layerNumber, barNumber);

      //ignore pulses from modules that have a number of layers other than 4.
      if(CRS->getCRSScintillatorShield(sectorNumber).getModule(moduleNumber).nLayers()!=CRVId::nLayers) continue;

      //sector properties
      std::map<int,sectorCoincidenceProperties>::const_iterator sIter = _sectorMap.find(sectorNumber);
      if(sIter==_sectorMap.end()) throw std::logic_error("CrvCoincidenceFinder: Found a CRV hit at a CRV sector without properties.");
      const sectorCoincidenceProperties &sector = sIter->second;

      //counter number within the entire sector type (like CRV-T, CRV-R, ...)
      int counterNumber = sector.precedingCounters + sector.nCountersPerModule*moduleNumber + barNumber;

      CLHEP::Hep3Vector crvCounterPos=CrvHelper::GetCrvCounterPos(CRS, crvBarIndex);
      double x=crvCounterPos[sector.widthDirection];
      double y=crvCounterPos[sector.thicknessDirection];

      //get the reco pulses information
      int    SiPM=crvRecoPulse->GetSiPMNumber();
      double time=crvRecoPulse->GetPulseTime();
      float  PEs =crvRecoPulse->GetPEs();

      //adjust pulse time to move it to longitudinal position that is farthest outside (at the same SiPM side)
      //this is done to be able to compare hits that happen at two adjacent sectors with different counter lengths
      int    side=SiPM%CRVId::nSidesPerBar;
      if(side==0)
      {
        double posDiff=sector.posSide0-_sectorTypeMap[sector.sectorType].posSide0;  //will always be >=0
        time+=posDiff/_fiberSignalSpeed;
      }
      else if(side==1)
      {
        double posDiff=_sectorTypeMap[sector.sectorType].posSide1-sector.posSide1;  //will always be >=0
        time+=posDiff/_fiberSignalSpeed;
      }

      //compensate for dead or ignored channels
      size_t currentChannel = crvBarIndex.asUint()*CRVId::nChanPerBar + SiPM;
      size_t testChannel = currentChannel + (SiPM-2<0 ? 2 : -2);  //the other channel at the same side of the counter
      std::bitset<16> status(sipmStatus.status(testChannel));
      bool compensate=false;
      for(size_t iChannelStatus=0; iChannelStatus<_compensateChannelStatus.size(); ++iChannelStatus)
      {
        if(status.test(_compensateChannelStatus.at(iChannelStatus))) compensate=true;  //test only channel status bits that are mentioned in the fcl file.
      }
      if(compensate) PEs *= 2.0; //double the PE value of this channel to compensate for the dead or ignored neighbor channel

      //don't split counter sides for the purpose of finding clusters
      sectorTypeMap[sector.sectorType].emplace_back(crvRecoPulse, crvCounterPos,
                                                    x, y, time, PEs, sectorNumber,
                                                    layerNumber, counterNumber, SiPM, sector.PEthreshold,
                                                    sector.maxTimeDifferenceAdjacentPulses, sector.maxTimeDifference,
                                                    sector.minSlope, sector.maxSlope, sector.maxSlopeDifference, sector.coincidenceLayers,
                                                    sector.minClusterPEs,sector.initialClusterMaxDistance);
    }

    //loop through all crv sectors types
    std::map<int, std::vector<CrvHit> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); ++sectorTypeMapIter)
    {
      int crvSectorType = sectorTypeMapIter->first;
      const std::vector<CrvHit> &hitsUnfiltered = sectorTypeMapIter->second;

      //filter hits, i.e. remove all hits below PE threshold
      //separate into both readout sides
      std::list<CrvHit> hitsFiltered0;
      std::list<CrvHit> hitsFiltered1;
      filterHits(hitsUnfiltered, hitsFiltered0, 0);
      filterHits(hitsUnfiltered, hitsFiltered1, 1);

      //maximum coincidence time difference between hits for this sector type
      //used to keep hits together in a cluster that could potentially form coincidences
      double clusterMaxTimeDifference=_sectorTypeMap.at(crvSectorType).maxTimeDifference;

      //distribute the hits into clusters
      //initial clustering is done to keep the number of hit combinations down that need to be checked for coincidences
      std::vector<std::vector<CrvHit> > clusters0;
      std::vector<std::vector<CrvHit> > clusters1;

      findClusters(hitsFiltered0, clusters0, clusterMaxTimeDifference);
      findClusters(hitsFiltered1, clusters1, clusterMaxTimeDifference);

      //all hits belonging to a coincidence group are collected in a new list
      std::list<CrvHit> coincidenceHits;

      //loop through all clusters and check whether they have coincidence groups
      //add each hit of all coincidence groups to a list of coincidence hits
      for(size_t iCluster=0; iCluster<clusters0.size(); ++iCluster) checkCoincidence(clusters0.at(iCluster),coincidenceHits);
      for(size_t iCluster=0; iCluster<clusters1.size(); ++iCluster) checkCoincidence(clusters1.at(iCluster),coincidenceHits);

      //maximum counter halflength and readout sides of this sector type
      float maxHalfLength=_sectorTypeMap.at(crvSectorType).maxHalfLength;
      bool  readoutSide0=_sectorTypeMap.at(crvSectorType).sipmsAtSide0;
      bool  readoutSide1=_sectorTypeMap.at(crvSectorType).sipmsAtSide1;

      //maximum time difference between hits of opposite ends for this sector type
      //based on the signal travel time through the max counter length
      //used to keep hits at opposite ends in one cluster (but only, if there are two readout sides)
      bool combineOppositeSides=false;
      if(readoutSide0 && readoutSide1) combineOppositeSides=true;
      double timeDifferenceOppositeSides = 2.0*maxHalfLength/_fiberSignalSpeed;

      //create new clusters based only on coincidence hits
      std::vector<std::vector<CrvHit> > coincidenceClusters;
      findClusters(coincidenceHits, coincidenceClusters, clusterMaxTimeDifference, combineOppositeSides, timeDifferenceOppositeSides);

      clusterProperties(crvSectorType, coincidenceClusters, crvCoincidenceClusterCollection, crvRecoPulseCollection);
    }//loop over all sector types

    ++_totalEvents;
    if(crvCoincidenceClusterCollection->size()>0) ++_totalEventsCoincidence;

    if(_verboseLevel>1)
    {
      std::cout<<moduleDescription().moduleLabel()<<"   run "<<event.id().run()<<"  subrun "<<event.id().subRun()<<"  event "<<event.id().event();
      std::cout<<"   Coincidence clusters found: "<<crvCoincidenceClusterCollection->size()<<std::endl;
    }

    event.put(std::move(crvCoincidenceClusterCollection));
  } // end produce


  void CrvCoincidenceFinder::clusterProperties(int crvSectorType, const std::vector<std::vector<CrvHit> > &clusters,
                                               std::unique_ptr<CrvCoincidenceClusterCollection> &crvCoincidenceClusterCollection,
                                               const art::Handle<CrvRecoPulseCollection> &crvRecoPulseCollection)
  {
    //loop through all clusters
    for(size_t iCluster=0; iCluster<clusters.size(); ++iCluster)
    {
      const std::vector<CrvHit> &cluster=clusters[iCluster];

      std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses;
      double startTime=cluster.front()._time;
      double endTime=cluster.front()._time;
      double PEs=0;
      CLHEP::Hep3Vector avgHitPos;  //PE-weighted average position
      std::set<int> layerSet;
      double sumX =0;
      double sumY =0;
      double sumYY=0;
      double sumXY=0;
      double minClusterPEs=cluster.front()._minClusterPEs;  //find the minimum of all minClusterPEs of the cluster hits
      for(auto hit=cluster.begin(); hit!=cluster.end(); ++hit)
      {
        crvRecoPulses.push_back(hit->_crvRecoPulse);

        PEs+=hit->_PEs;
        avgHitPos+=hit->_pos*hit->_PEs;
        layerSet.insert(hit->_layer);
        sumX +=hit->_PEs*hit->_x;
        sumY +=hit->_PEs*hit->_y;
        sumYY+=hit->_PEs*hit->_y*hit->_y;
        sumXY+=hit->_PEs*hit->_x*hit->_y;

        if(startTime>hit->_time) startTime=hit->_time;
        if(endTime<hit->_time) endTime=hit->_time;
        if(minClusterPEs>hit->_minClusterPEs) minClusterPEs=hit->_minClusterPEs;
      } //loop over hits of the cluster

      assert(PEs>0);
      assert(layerSet.size()>1);

      //average counter position (PE weighted), slope, layers
      avgHitPos/=PEs;
      double slope=(PEs*sumXY-sumX*sumY)/(PEs*sumYY-sumY*sumY);
      std::vector<int> layers(layerSet.begin(), layerSet.end());

      //don't store clusters that are below the minimum number of PEs for this sector (or sectors, if the cluster involves multiple sectors)
      if(PEs<minClusterPEs) continue;

      //some values we want to store for the cluster
      std::array<size_t, CRVId::nSidesPerBar> sideHits{0,0};
      std::array<float, CRVId::nSidesPerBar>  sidePEs{0,0};
      std::array<double, CRVId::nSidesPerBar> sideTimes{0,0};
      std::array<double, CRVId::nSidesPerBar> sideTimesMin{0,0};
      std::array<double, CRVId::nSidesPerBar> sideTimesMax{0,0};
      int lengthDirection=_sectorMap.at(cluster.begin()->_crvSector).lengthDirection;   //same for all counters of this cluster
      for(auto hit=cluster.begin(); hit!=cluster.end(); ++hit)
      {
        if(hit->_PEs<=0.0) continue;  //avoid dealing with hits that have no PEs

        int side=hit->_SiPM%CRVId::nSidesPerBar;

        sideHits[side]++;
        sidePEs[side]+=hit->_PEs;
        sideTimes[side]+=hit->_time*hit->_PEs; //PE-weighted times
        if(sideHits[side]==1) {sideTimesMin[side]=hit->_time; sideTimesMax[side]=hit->_time;}
        if(sideTimesMin[side]>hit->_time) sideTimesMin[side]=hit->_time;
        if(sideTimesMax[side]<hit->_time) sideTimesMax[side]=hit->_time;
      }
      if(sidePEs[0]>0) sideTimes[0]/=sidePEs[0]; //find avg PE-weighted time by dividing by total PEs
      if(sidePEs[1]>0) sideTimes[1]/=sidePEs[1];

      bool hitPosAndTimeCalculated=false;
      double avgHitTime=0.0;
      if(sideHits[0]>0 && sideHits[1]>0)  //both readout sides have hits
      {
        double posSide0=_sectorTypeMap[crvSectorType].posSide0;
        double posSide1=_sectorTypeMap[crvSectorType].posSide1;

        double posOrigin=0.5*(_fiberSignalSpeed*(sideTimes[0]-sideTimes[1])+(posSide0+posSide1));

        double timeInterval0=sideTimesMax[0]-sideTimesMin[0];
        double timeInterval1=sideTimesMax[1]-sideTimesMin[1];

        if(posOrigin>posSide0 && posOrigin<posSide1 && posSide0<posSide1 &&    //hit origin is within the counter length and positions of
           timeInterval0<_maxTimeInterval && timeInterval1<_maxTimeInterval)   //readout sides are not inverted (should not happen in our geometry)
                                                                               //and hit time intervals at both sides are within a limit
        {
          hitPosAndTimeCalculated=true;
          avgHitPos[lengthDirection]=posOrigin; //correct the longitudinal coordinate of the avgHitPos
          avgHitTime=sideTimes[0]-(posOrigin-posSide0)/_fiberSignalSpeed; //time of hit origin
          //can also be calculated from the other side in the following way:
          //avgHitTime = sideTimes[1]-(posSide1-posOrigin)/_fiberSignalSpeed
          //both results are identical.
        }
        else //hits at both readout don't seem to be correlated (and would result in a hit origin outside the counter)
        {
          //assume the hit origin is at the center (since no other information is available)
          posOrigin=avgHitPos[lengthDirection];
          //calculate the time of the hit origin from both sides
          double avgHitTime0=sideTimes[0]-(posOrigin-posSide0)/_fiberSignalSpeed;
          double avgHitTime1=sideTimes[1]-(posSide1-posOrigin)/_fiberSignalSpeed;
          //and take the average. that's the best we can do
          //(other than breaking this cluster apart. we may do this later.)
          avgHitTime=0.5*(avgHitTime0+avgHitTime1);
        }
      }
      else //only one readout side has hits
      {
        //assume the hit origin is at the center (since no other information is available)
        double posOrigin=avgHitPos[lengthDirection];
        if(sideHits[0]>0)
        {
          double posSide0=_sectorTypeMap[crvSectorType].posSide0;
          avgHitTime=sideTimes[0]-(posOrigin-posSide0)/_fiberSignalSpeed;
        }
        else if(sideHits[1]>0)
        {
          double posSide1=_sectorTypeMap[crvSectorType].posSide1;
          avgHitTime=sideTimes[1]-(posSide1-posOrigin)/_fiberSignalSpeed;
        }
      }

      //remove additional time due to electronics response and processes such as scintillation/WLS decay times
      avgHitTime-=_timeOffset;   //we could use a non-const parameter by taking the pulse widths of this cluster into consideration

      //insert the cluster information into the vector of the crv coincidence clusters
      crvCoincidenceClusterCollection->emplace_back(crvSectorType, startTime, endTime, PEs, crvRecoPulses, slope, layers,
                                                    sideHits, sidePEs, sideTimes, avgHitTime, avgHitPos,
                                                    hitPosAndTimeCalculated);
    } //loop through all clusters
  } //end cluster properies


  //remove hits below the threshold and separate into both readout sides
  void CrvCoincidenceFinder::filterHits(const std::vector<CrvHit> &hits, std::list<CrvHit> &hitsFiltered, unsigned int readoutSide)
  {
    std::vector<CrvHit>::const_iterator iterHit;
    for(iterHit=hits.begin(); iterHit!=hits.end(); ++iterHit)
    {
      int    SiPM=iterHit->_SiPM;
      if(SiPM%CRVId::nSidesPerBar!=readoutSide) continue;  //not the readout side we are looking for

      int    layer=iterHit->_layer;
      int    counter=iterHit->_counter;  //counter number in one layer counted from the beginning of the counter type
      double time=iterHit->_time;

      int    PEthreshold=iterHit->_PEthreshold;
      double maxTimeDifferenceAdjacentPulses=iterHit->_maxTimeDifferenceAdjacentPulses;

      //check other SiPM and the SiPMs at the adjacent counters
      double PEs_thisCounter=0;
      double PEs_adjacentCounter1=0;
      double PEs_adjacentCounter2=0;
      std::vector<CrvHit>::const_iterator iterHitAdjacent;
      for(iterHitAdjacent=hits.begin(); iterHitAdjacent!=hits.end(); ++iterHitAdjacent)
      {
        //use hits of the same readout side only
        if(iterHitAdjacent->_SiPM%CRVId::nSidesPerBar!=readoutSide) continue;

        //use hits of the same layer only
        if(iterHitAdjacent->_layer!=layer) continue;

        //use hits within a certain time window only
        if(fabs(iterHitAdjacent->_time-time)>maxTimeDifferenceAdjacentPulses) continue;

        //collect all PEs of this and the adjacent counters
        int counterDiff=iterHitAdjacent->_counter-counter;
        if(counterDiff==0) PEs_thisCounter+=iterHitAdjacent->_PEs;   //add PEs from the same counter (i.e. the "other" SiPM),
                                                                     //if the "other" hit is within a certain time window (5ns)
                                                                     //this will include PEs from the current pulse
        if(counterDiff==-1) PEs_adjacentCounter1+=iterHitAdjacent->_PEs;  //add PEs from an adjacent counter,
                                                                          //if these hits are within a certain time window (5ns)
        if(counterDiff==1) PEs_adjacentCounter2+=iterHitAdjacent->_PEs;   //add PEs from an adjacent counter,
                                                                          //if these hits are within a certain time window (5ns)
      }

      //if the number of PEs of this hit (plus the number of PEs of the same or one of the adjacent counter, if their time
      //difference is small enough) is above the PE threshold, add this hit to vector of filtered hits
      if(PEs_thisCounter+PEs_adjacentCounter1>=PEthreshold || PEs_thisCounter+PEs_adjacentCounter2>=PEthreshold)
         hitsFiltered.push_back(*iterHit);
    }
  } //end filter hits


  void CrvCoincidenceFinder::findClusters(std::list<CrvHit> &hits, std::vector<std::vector<CrvHit> > &clusters,
                                          double clusterMaxTimeDifference, bool combineOppositeSides, double timeDifferenceOppositeSides)
  {
    while(!hits.empty()) //run through clustering processes until all hits are distributed into clusters
    {
      clusters.resize(clusters.size()+1); //add a new cluster
      std::vector<CrvHit> &cluster = clusters.back();

      //need to loop several times to check the unused hits until the cluster size remains stable
      size_t lastClusterSize=0;
      do
      {
        lastClusterSize=cluster.size();

        for(auto hitsIter=hits.begin(); hitsIter!=hits.end(); )
        {
          if(cluster.empty())
          {
            //first element of current cluster
            //gets moved from list of hits to the current cluster
            cluster.push_back(*hitsIter);
            hitsIter=hits.erase(hitsIter);
          }
          else
          {
            //check whether current hit satisfies time and distance condition w.r.t. to any hit of current cluster
            bool erasedHit=false;
            for(auto clusterIter=cluster.begin(); clusterIter!=cluster.end(); ++clusterIter)
            {
              double maxDistance = std::max(hitsIter->_maxDistance,clusterIter->_maxDistance);
              double maxTimeDifference = clusterMaxTimeDifference;
              if(combineOppositeSides)
              {
                if((hitsIter->_SiPM%CRVId::nSidesPerBar)!=(clusterIter->_SiPM%CRVId::nSidesPerBar)) //both hits are at opposite sides --> use additional time difference
                  maxTimeDifference+=timeDifferenceOppositeSides;
              }
              if((std::fabs(hitsIter->_x-clusterIter->_x)<=maxDistance) &&
                 (std::fabs(hitsIter->_time-clusterIter->_time)<maxTimeDifference))
              {
                //this hit satisfied the conditions
                //move it from list of hits to the current cluster
                cluster.push_back(*hitsIter);
                hitsIter=hits.erase(hitsIter);
                erasedHit=true;
                break;  //no need for more comparisons with other hits in current cluster, go to the next hit in list of hits
              }
            } //loop over all hits in the cluster (for comparison with current hit)

            if(!erasedHit) ++hitsIter;

          } //cluster not empty
        } //loop over all undistributed hits

      } while(lastClusterSize!=cluster.size()); //loop until cluster does not change anymore

    } //loop until all hits in list of hits are distributed into clusters
  } //end finder clusters


  void CrvCoincidenceFinder::checkCoincidence(const std::vector<CrvHit> &hits, std::list<CrvHit> &coincidenceHits)
  {
    if(hits.empty()) return;

    std::vector<CrvHit> hitsLayers[CRVId::nLayers];  //separated by layers
    std::vector<CrvHit>::const_iterator iterHit;
    for(iterHit=hits.begin(); iterHit!=hits.end(); ++iterHit)
    {
      int    layer=iterHit->_layer;
      hitsLayers[layer].push_back(*iterHit);
    }

    //we want to collect all hits belonging to coincidence groups,
    //but avoid collecting hits multiple times, if they belong to different coincidence groups.
    //can be done by placing the hit interator into a set.
    auto setComp = [](const std::vector<CrvHit>::const_iterator &a, const std::vector<CrvHit>::const_iterator &b) {return a->_crvRecoPulse < b->_crvRecoPulse;};
    std::set<std::vector<CrvHit>::const_iterator,decltype(setComp)> coincidenceHitSet(setComp);

    int minCoincidenceLayers = std::min_element(hits.begin(),hits.end(),
                               [](const CrvHit &a, const CrvHit &b){return a._coincidenceLayers < b._coincidenceLayers;})->_coincidenceLayers;
    int maxCoincidenceLayers = std::max_element(hits.begin(),hits.end(),
                               [](const CrvHit &a, const CrvHit &b){return a._coincidenceLayers < b._coincidenceLayers;})->_coincidenceLayers;

    if(hits.size()>_bigClusterThreshold)
    {
      //this cluster has so many hits that it makes no sense anymore to search for individual coincidences.
      //we still need to check that the minimum number of layers were hit to skip the coincidence check.
      int nonEmptyLayers=0;
      for(size_t iLayer=0; iLayer<CRVId::nLayers; ++iLayer)
      {
        if(!hitsLayers[iLayer].empty()) ++nonEmptyLayers;
      }
      if(nonEmptyLayers>=minCoincidenceLayers)
      {
        for(auto iterHit=hits.begin(); iterHit!=hits.end(); ++iterHit) coincidenceHits.push_back(*iterHit);
        return;
      }
    }

    //***************************************************
    //find coincidences using 2/4 coincidence requirement
    if(minCoincidenceLayers==2)
    {
      std::vector<CrvHit>::const_iterator layerIterators[2];

      for(int layer1=0; layer1<4; ++layer1)
      for(int layer2=layer1+1; layer2<4; ++layer2)
      {
        const std::vector<CrvHit> &layer1Hits=hitsLayers[layer1];
        const std::vector<CrvHit> &layer2Hits=hitsLayers[layer2];
        std::vector<CrvHit>::const_iterator layer1Iter;
        std::vector<CrvHit>::const_iterator layer2Iter;

        //it will loop only, if both layers have hits
        //loops will exit, if all three hits require a 3/4 or 4/4 coincidence
        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); ++layer1Iter)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); ++layer2Iter)
        {
          if(layer1Iter->_coincidenceLayers>2 && layer2Iter->_coincidenceLayers>2) continue; //all hits require at least a 3/4 coincidence

          //add all both hits into one array
          layerIterators[0]=layer1Iter;
          layerIterators[1]=layer2Iter;
          if(checkCombination(layerIterators,2))
          {
            coincidenceHitSet.insert(layer1Iter);
            coincidenceHitSet.insert(layer2Iter);
          }
        }
      }
    }  //two layer coincidences

    //***************************************************
    //find coincidences using 3/4 coincidence requirement
    if(minCoincidenceLayers<=3 && maxCoincidenceLayers>=3)
    {
      std::vector<CrvHit>::const_iterator layerIterators[3];

      for(int layer1=0; layer1<4; ++layer1)
      for(int layer2=layer1+1; layer2<4; ++layer2)
      for(int layer3=layer2+1; layer3<4; ++layer3)
      {
        const std::vector<CrvHit> &layer1Hits=hitsLayers[layer1];
        const std::vector<CrvHit> &layer2Hits=hitsLayers[layer2];
        const std::vector<CrvHit> &layer3Hits=hitsLayers[layer3];
        std::vector<CrvHit>::const_iterator layer1Iter;
        std::vector<CrvHit>::const_iterator layer2Iter;
        std::vector<CrvHit>::const_iterator layer3Iter;

        //it will loop only, if all 3 layers have hits
        //loops will exit, if all three hits require a 4/4 coincidence
        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); ++layer1Iter)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); ++layer2Iter)
        for(layer3Iter=layer3Hits.begin(); layer3Iter!=layer3Hits.end(); ++layer3Iter)
        {
          if(layer1Iter->_coincidenceLayers>3 && layer2Iter->_coincidenceLayers>3 && layer3Iter->_coincidenceLayers>3) continue; //all hits require at 4/4 coincidence

          //add all 3 hits into one array
          layerIterators[0]=layer1Iter;
          layerIterators[1]=layer2Iter;
          layerIterators[2]=layer3Iter;
          if(checkCombination(layerIterators,3))
          {
            coincidenceHitSet.insert(layer1Iter);
            coincidenceHitSet.insert(layer2Iter);
            coincidenceHitSet.insert(layer3Iter);
          }
        }
      }
    }  //three layer coincidences

    //***************************************************
    //find coincidences using 4/4 coincidence requirement
    if(maxCoincidenceLayers==4)
    {
      std::vector<CrvHit>::const_iterator layerIterators[4];

      const std::vector<CrvHit> &layer0Hits=hitsLayers[0];
      const std::vector<CrvHit> &layer1Hits=hitsLayers[1];
      const std::vector<CrvHit> &layer2Hits=hitsLayers[2];
      const std::vector<CrvHit> &layer3Hits=hitsLayers[3];
      std::vector<CrvHit>::const_iterator layer0Iter;
      std::vector<CrvHit>::const_iterator layer1Iter;
      std::vector<CrvHit>::const_iterator layer2Iter;
      std::vector<CrvHit>::const_iterator layer3Iter;

      //it will loop only, if all 4 layers have hits
      for(layer0Iter=layer0Hits.begin(); layer0Iter!=layer0Hits.end(); ++layer0Iter)
      for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); ++layer1Iter)
      for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); ++layer2Iter)
      for(layer3Iter=layer3Hits.begin(); layer3Iter!=layer3Hits.end(); ++layer3Iter)
      {
        //add all 4 hits into one array
        layerIterators[0]=layer0Iter;
        layerIterators[1]=layer1Iter;
        layerIterators[2]=layer2Iter;
        layerIterators[3]=layer3Iter;
        if(checkCombination(layerIterators,4))
        {
          coincidenceHitSet.insert(layer0Iter);
          coincidenceHitSet.insert(layer1Iter);
          coincidenceHitSet.insert(layer2Iter);
          coincidenceHitSet.insert(layer3Iter);
        }
      }
    } // four layer coincidences

    //move the set of coincidence hit iterators to a list of hits
    for(auto iterHit=coincidenceHitSet.begin(); iterHit!=coincidenceHitSet.end(); ++iterHit) coincidenceHits.push_back(**iterHit);

  } //end check coincidence

  bool CrvCoincidenceFinder::checkCombination(std::vector<CrvHit>::const_iterator layerIterators[], int n)
  {
    typedef const std::vector<CrvHit>::const_iterator L;

    double maxTimeDifference = (*std::max_element(layerIterators,layerIterators+n,
                               [](L &a, L &b){return a->_maxTimeDifference < b->_maxTimeDifference;}))->_maxTimeDifference;
    double timeMax = (*std::max_element(layerIterators,layerIterators+n,
                     [](L &a, L &b){return a->_time < b->_time;}))->_time;
    double timeMin = (*std::min_element(layerIterators,layerIterators+n,
                     [](L &a, L &b){return a->_time < b->_time;}))->_time;
    if(timeMax-timeMin>maxTimeDifference) return false;  //hits don't fall within the time window

    double minSlope = (*std::min_element(layerIterators,layerIterators+n,
                      [](L &a, L &b){return a->_minSlope < b->_minSlope;}))->_minSlope;
    double maxSlope = (*std::max_element(layerIterators,layerIterators+n,
                      [](L &a, L &b){return a->_maxSlope < b->_maxSlope;}))->_maxSlope;
    double maxSlopeDifference = (*std::max_element(layerIterators,layerIterators+n,
                                [](L &a, L &b){return a->_maxSlopeDifference < b->_maxSlopeDifference;}))->_maxSlopeDifference;
    std::array<double,CRVId::nLayers-1> slopes;
    for(int d=0; d<n-1; ++d)
    {
      //slope = width direction / thickness direction
      slopes[d]=(layerIterators[d+1]->_x-layerIterators[d]->_x)/(layerIterators[d+1]->_y-layerIterators[d]->_y);
      if(slopes[d]<minSlope || slopes[d]>maxSlope) return false; //slopes need to be within minSlope and maxSlope
    }

    if(n>2)
    {
      if(fabs(slopes[0]-slopes[1])>maxSlopeDifference) return false;
      if(n>3)
      {
        if(fabs(slopes[0]-slopes[2])>maxSlopeDifference) return false;
        if(fabs(slopes[1]-slopes[2])>maxSlopeDifference) return false;
      }
    }

    //need distances between subsequent positions, when ordered by position
    std::sort(layerIterators,layerIterators+n,[](L &a, L &b){return a->_x < b->_x;});
    for(int i=0; i<n; ++i)
    {
      //find the max distances between hits belonging to a coincidence group.
      //these max distances are used later when new clusters based on coincidence hits are created,
      //because the distance between hits belonging to a coincidence group
      //may be greater than the original clusterMaxDistance (if it was chosen very small).
      //this is done to avoid breaking coincidence groups apart during the next clustering process
      if(i>0)
      {
        double distance=layerIterators[i]->_x-layerIterators[i-1]->_x;
        if(distance>layerIterators[i]->_maxDistance) layerIterators[i]->_maxDistance=distance;
      }
      if(i+1<n)
      {
        double distance=layerIterators[i+1]->_x-layerIterators[i]->_x;
        if(distance>layerIterators[i]->_maxDistance) layerIterators[i]->_maxDistance=distance;
      }
    }

    return true; //all coincidence criteria for this combination satisfied
  } // end checkCombination

} // end namespace mu2e

using mu2e::CrvCoincidenceFinder;
DEFINE_ART_MODULE(CrvCoincidenceFinder)
