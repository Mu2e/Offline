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
      fhicl::Atom<double> minOverlapTimeAdjacentPulses{Name("minOverlapTimeAdjacentPulses"), Comment("minimum overlap time between pulses of adjacent channels to be considered for coincidences")};
      fhicl::Atom<double> minOverlapTime{Name("minOverlapTime"), Comment("minimum overlap time between pulses of a coincidence hit combination")};
      fhicl::Atom<double> minSlope{Name("minSlope"), Comment("minimum slope allowed for a coincidence")};
      fhicl::Atom<double> maxSlope{Name("maxSlope"), Comment("maximum slope allowed for a coincidence")};
      fhicl::Atom<double> maxSlopeDifference{Name("maxSlopeDifference"), Comment("maximum slope difference between layers allowed for a coincidence")};
      fhicl::Atom<int> coincidenceLayers{Name("coincidenceLayers"), Comment("number of layers required for a coincidence")};
      fhicl::Atom<double> minClusterPEs{Name("minClusterPEs"), Comment("minimum number of PEs in a cluster required for storage")};
    };
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> verboseLevel{Name("verboseLevel"), Comment("verbose level")};
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("module label of the input CrvRecoPulses")};
      //cluster settings
      fhicl::Atom<double> initialClusterMaxDistance{Name("initialClusterMaxDistance"), Comment("maximum distances between hits to be considered for a hit cluster at initial clustering process")};
      fhicl::Atom<double> clusterMaxTimeDifference{Name("clusterMaxTimeDifference"), Comment("maximum time difference between hits to be considered for a hit cluster (when overlap option is not used)")};
      fhicl::Atom<double> clusterMinOverlapTime{Name("clusterMinOverlapTime"), Comment("minimum overlap time between hits to be considered for a hit cluster (when overlap option is used)")};
      //sector-specific coincidence settings
      fhicl::Sequence<fhicl::Table<SectorConfig> > sectorConfig{Name("sectorConfig"), Comment("sector-specific settings")};
      //coincidence settings for overlap option
      fhicl::Atom<bool> usePulseOverlaps{Name("usePulseOverlaps"), Comment("use pulse overlaps instead of peak times to determine coincidences")};
      //other settings
      fhicl::Atom<bool> useNoFitReco{Name("useNoFitReco"), Comment("use pulse reco results not based on a Gumbel fit")};
      fhicl::Atom<bool> usePEsPulseHeight{Name("usePEsPulseHeight"), Comment("use PEs determined by pulse height instead of pulse area")};
      fhicl::Atom<int> bigClusterThreshold{Name("bigClusterThreshold"), Comment("no coincidence check for clusters with a number of hits above this threshold")};
      fhicl::Atom<double> fiberSignalSpeed{Name("fiberSignalSpeed"), Comment("effective speed of signals inside the CRV fibers in mm/ns")};
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

    double                   _initialClusterMaxDistance;
    double                   _initialClusterMaxTimeDifference;
    double                   _initialClusterMinOverlapTime;
    double                   _clusterMaxTimeDifference;
    double                   _clusterMinOverlapTime;

    std::vector<SectorConfig> _sectorConfig;

    bool        _usePulseOverlaps;
    bool        _useNoFitReco;
    bool        _usePEsPulseHeight;
    size_t      _bigClusterThreshold;
    double      _fiberSignalSpeed;
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
      int         PEthreshold;
      double      maxTimeDifferenceAdjacentPulses;
      double      maxTimeDifference;
      double      minOverlapTimeAdjacentPulses;
      double      minOverlapTime;
      double      minSlope, maxSlope, maxSlopeDifference;
      int         coincidenceLayers;
      double      minClusterPEs;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;

    struct CrvHit
    {
      art::Ptr<CrvRecoPulse> _crvRecoPulse;
      CLHEP::Hep3Vector      _pos;
      double _x, _y;
      double _time, _timePulseStart, _timePulseEnd;
      double _PEs;
      int    _crvSector;
      int    _layer;
      int    _counter;
      int    _SiPM;
      int    _PEthreshold;
      double _maxTimeDifferenceAdjacentPulses;
      double _maxTimeDifference;
      double _minOverlapTimeAdjacentPulses;
      double _minOverlapTime;
      double _minSlope, _maxSlope, _maxSlopeDifference;
      int    _coincidenceLayers;
      double _minClusterPEs;
      mutable double _maxDistance; //initially set to initialClusterMaxDistance, which is just an estimate
                                   //used for the initial clustering process (to keep the number of
                                   //hit combinations down that need to be checked for coincidendes).
                                   //_maxDistance is updated when the coindidences are checked
                                   //(used for the final clustering process)

      CrvHit(const art::Ptr<CrvRecoPulse> crvRecoPulse, const CLHEP::Hep3Vector &pos,
             double x, double y, double time, double timePulseStart, double timePulseEnd,
             double PEs, int crvSector, int layer, int counter, int SiPM, int PEthreshold,
             double maxTimeDifferenceAdjacentPulses, double maxTimeDifference,
             double minOverlapTimeAdjacentPulses, double minOverlapTime,
             double minSlope, double maxSlope, double maxSlopeDifference, int coincidenceLayers, double minClusterPEs, double maxDistance) :
               _crvRecoPulse(crvRecoPulse), _pos(pos),
               _x(x), _y(y), _time(time), _timePulseStart(timePulseStart), _timePulseEnd(timePulseEnd),
               _PEs(PEs), _crvSector(crvSector), _layer(layer), _counter(counter), _SiPM(SiPM), _PEthreshold(PEthreshold),
               _maxTimeDifferenceAdjacentPulses(maxTimeDifferenceAdjacentPulses), _maxTimeDifference(maxTimeDifference),
               _minOverlapTimeAdjacentPulses(minOverlapTimeAdjacentPulses), _minOverlapTime(minOverlapTime),
               _minSlope(minSlope), _maxSlope(maxSlope), _maxSlopeDifference(maxSlopeDifference), _coincidenceLayers(coincidenceLayers),
               _minClusterPEs(minClusterPEs), _maxDistance(maxDistance) {}
    };

    void clusterProperties(int crvSectorType, const std::vector<std::vector<CrvHit> > &clusters,
                           std::unique_ptr<CrvCoincidenceClusterCollection> &crvCoincidenceClusterCollection,
                           const art::Handle<CrvRecoPulseCollection> &crvRecoPulseCollection);
    void filterHits(const std::vector<CrvHit> &hits, std::list<CrvHit> &hitsFiltered);
    void findClusters(std::list<CrvHit> &hits, std::vector<std::vector<CrvHit> > &clusters,
                      double clusterMaxTimeDifference, double clusterMinOverlapTime);
    void checkCoincidence(const std::vector<CrvHit> &hits, std::list<CrvHit> &coincidenceHits);
    bool checkCombination(std::vector<CrvHit>::const_iterator layerIterators[], int n);

  };

  CrvCoincidenceFinder::CrvCoincidenceFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _verboseLevel(conf().verboseLevel()),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _initialClusterMaxDistance(conf().initialClusterMaxDistance()),
    _clusterMaxTimeDifference(conf().clusterMaxTimeDifference()),
    _clusterMinOverlapTime(conf().clusterMinOverlapTime()),
    _sectorConfig(conf().sectorConfig()),
    _usePulseOverlaps(conf().usePulseOverlaps()),
    _useNoFitReco(conf().useNoFitReco()),
    _usePEsPulseHeight(conf().usePEsPulseHeight()),
    _bigClusterThreshold(conf().bigClusterThreshold()),
    _fiberSignalSpeed(conf().fiberSignalSpeed()),
    _timeOffset(conf().timeOffset()),
    _compensateChannelStatus(conf().compensateChannelStatus()),
    _totalEvents(0),
    _totalEventsCoincidence(0)
  {
    produces<CrvCoincidenceClusterCollection>();
    //get initial cluster time parameters from coincidence parameters
    //initial clustering is done to keep the number of hit combinations down that need to be checked for coincidences
    _initialClusterMaxTimeDifference=std::max_element(_sectorConfig.begin(), _sectorConfig.end(),
                              [](const SectorConfig &a, const SectorConfig &b)
                              {return a.maxTimeDifference() < b.maxTimeDifference();})->maxTimeDifference();
    _initialClusterMinOverlapTime=std::min_element(_sectorConfig.begin(), _sectorConfig.end(),
                              [](const SectorConfig &a, const SectorConfig &b)
                              {return a.minOverlapTime() < b.minOverlapTime();})->minOverlapTime();
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

      std::string sectorName = sectors[i].getName().substr(4); //removes the "CRV_" part
      std::vector<SectorConfig>::iterator sectorConfigIter=std::find_if(_sectorConfig.begin(), _sectorConfig.end(),
                                                                        [sectorName](const SectorConfig &s){return s.CRVSector()==sectorName;});
      if(sectorConfigIter==_sectorConfig.end())
        throw std::logic_error("CrvCoincidenceFinder: The geometry has a CRV sector for which no coincidence properties were defined in the fcl file.");

      s.PEthreshold                     = sectorConfigIter->PEthreshold();
      s.maxTimeDifferenceAdjacentPulses = sectorConfigIter->maxTimeDifferenceAdjacentPulses();
      s.maxTimeDifference               = sectorConfigIter->maxTimeDifference();
      s.minOverlapTimeAdjacentPulses    = sectorConfigIter->minOverlapTimeAdjacentPulses();
      s.minOverlapTime                  = sectorConfigIter->minOverlapTime();
      s.minSlope                        = sectorConfigIter->minSlope();
      s.maxSlope                        = sectorConfigIter->maxSlope();
      s.maxSlopeDifference              = sectorConfigIter->maxSlopeDifference();
      s.coincidenceLayers               = sectorConfigIter->coincidenceLayers();
      s.minClusterPEs                   = sectorConfigIter->minClusterPEs();

      _sectorMap[i]=s;
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

      if(_usePulseOverlaps && crvRecoPulse->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::duplicateNoFitPulse)) continue;

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
      int SiPM = crvRecoPulse->GetSiPMNumber();

      double time=crvRecoPulse->GetPulseTime();
      double timePulseStart=crvRecoPulse->GetPulseStart();
      double timePulseEnd=crvRecoPulse->GetPulseEnd();
      float  PEs =crvRecoPulse->GetPEs();
      if(_usePEsPulseHeight) PEs=crvRecoPulse->GetPEsPulseHeight();
      if(_useNoFitReco) PEs=crvRecoPulse->GetPEsNoFit();
      if(_useNoFitReco) time=crvRecoPulse->GetPulseTimeNoFit();
      if(_usePulseOverlaps) PEs=crvRecoPulse->GetPEsNoFit();

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
                                                    x, y, time, timePulseStart, timePulseEnd, PEs, sectorNumber,
                                                    layerNumber, counterNumber, SiPM, sector.PEthreshold,
                                                    sector.maxTimeDifferenceAdjacentPulses, sector.maxTimeDifference,
                                                    sector.minOverlapTimeAdjacentPulses, sector.minOverlapTime,
                                                    sector.minSlope, sector.maxSlope, sector.maxSlopeDifference, sector.coincidenceLayers,
                                                    sector.minClusterPEs,_initialClusterMaxDistance);
    }

    //loop through all crv sectors types
    std::map<int, std::vector<CrvHit> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); ++sectorTypeMapIter)
    {
      int crvSectorType = sectorTypeMapIter->first;
      const std::vector<CrvHit> &hitsUnfiltered = sectorTypeMapIter->second;

      //filter hits, i.e. remove all hits below PE threshold
      std::list<CrvHit> hitsFiltered;
      filterHits(hitsUnfiltered, hitsFiltered);

      //distribute the hits into clusters
      //initial clustering is done to keep the number of hit combinations down that need to be checked for coincidences
      std::vector<std::vector<CrvHit> > clusters;
      findClusters(hitsFiltered, clusters, _initialClusterMaxTimeDifference, _initialClusterMinOverlapTime);

      //all hits belonging to a coincidence group are collected in a new list
      std::list<CrvHit> coincidenceHits;

      //loop through all clusters
      for(size_t iCluster=0; iCluster<clusters.size(); ++iCluster)
      {
        const std::vector<CrvHit> &cluster=clusters[iCluster];
        std::vector<CrvHit> cluster0;  //cluster for SiPMs at negative end
        std::vector<CrvHit> cluster1;  //cluster for SiPMs at positive end
        for(auto hit=cluster.begin(); hit!=cluster.end(); ++hit)
        {
          if(hit->_SiPM%CRVId::nSidesPerBar==0) cluster0.push_back(*hit); else cluster1.push_back(*hit);
        }

        //check whether this hit cluster has coincidences
        //(separately for both readout sides)
        checkCoincidence(cluster0,coincidenceHits);
        checkCoincidence(cluster1,coincidenceHits);
      }//loop over all cluster in sector type

      //create new clusters based only on coincidence hits
      std::vector<std::vector<CrvHit> > coincidenceClusters;
      findClusters(coincidenceHits, coincidenceClusters, _clusterMaxTimeDifference, _clusterMinOverlapTime);

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
      if(_usePulseOverlaps)
      {
        startTime=cluster.front()._timePulseStart;
        endTime=cluster.front()._timePulseEnd;
      }
      double PEs=0;
      CLHEP::Hep3Vector avgCounterPos;  //PE-weighted average position
      std::set<int> layerSet;
      double sumX =0;
      double sumY =0;
      double sumYY=0;
      double sumXY=0;
      double minClusterPEs=cluster.front()._minClusterPEs;  //find the minimum of all minClusterPEs of the cluster hits
      std::map<int,std::vector<std::vector<CrvHit>::const_iterator> > recoTimes;  //collect hits of each counter
      for(auto hit=cluster.begin(); hit!=cluster.end(); ++hit)
      {
        crvRecoPulses.push_back(hit->_crvRecoPulse);
        if(_usePulseOverlaps)
        {
          //duplicate nofit pulses are removed in the usePulseOverlap option, but should be included in the list of reco pulses
          for(size_t recoPulseIndex=hit->_crvRecoPulse.key()+1; recoPulseIndex<crvRecoPulseCollection->size(); ++recoPulseIndex)
          {
            const art::Ptr<CrvRecoPulse> crvRecoPulse(crvRecoPulseCollection, recoPulseIndex);
            if(!crvRecoPulse->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::duplicateNoFitPulse)) break;
            crvRecoPulses.push_back(crvRecoPulse);
          }
        }

        PEs+=hit->_PEs;
        avgCounterPos+=hit->_pos*hit->_PEs;
        layerSet.insert(hit->_layer);
        sumX +=hit->_PEs*hit->_x;
        sumY +=hit->_PEs*hit->_y;
        sumYY+=hit->_PEs*hit->_y*hit->_y;
        sumXY+=hit->_PEs*hit->_x*hit->_y;
        if(_usePulseOverlaps)
        {
          if(startTime>hit->_timePulseStart) startTime=hit->_timePulseStart;
          if(endTime<hit->_timePulseEnd) endTime=hit->_timePulseEnd;
        }
        else
        {
          if(startTime>hit->_time) startTime=hit->_time;
          if(endTime<hit->_time) endTime=hit->_time;
        }

        if(minClusterPEs>hit->_minClusterPEs) minClusterPEs=hit->_minClusterPEs;

        //collect hits of individual counters
        recoTimes[hit->_crvRecoPulse->GetScintillatorBarIndex().asInt()].push_back(hit);
      } //loop over hits of the cluster

      assert(PEs>0);
      assert(layerSet.size()>1);

      //average counter position (PE weighted), slope, layers
      avgCounterPos/=PEs;
      double slope=(PEs*sumXY-sumX*sumY)/(PEs*sumYY-sumY*sumY);
      std::vector<int> layers(layerSet.begin(), layerSet.end());

      //don't store clusters that are below the minimum number of PEs for this sector (or sectors, if the cluster involves multiple sectors)
      if(PEs<minClusterPEs) continue;

      //find average hit times and positions
      std::array<size_t, CRVId::nSidesPerBar> sideHits{0,0};
      std::array<float, CRVId::nSidesPerBar>  sidePEs{0,0};
      std::array<double, CRVId::nSidesPerBar> sideTimes{0,0};  //this value becomes meaningless, if a coincidence cluster spans more than one sector
                                                               //with different counter lengths. (not a problem for extrated position)
      double                                  avgHalfLength{0};  //needed, if only one readout side is available, and the center of the counter is used
      double                                  avgHitTime{0};
      CLHEP::Hep3Vector                       avgHitPos;
      double                                  PEbothSides{0};  //number of PEs for counters with readouts on both sides (needed to PE-weight the counter hits)

      //calculate average times and positions for each counter separately
      for(auto recoTimeIter=recoTimes.begin(); recoTimeIter!=recoTimes.end(); ++recoTimeIter)
      {
        int crvSector=recoTimeIter->second.at(0)->_crvSector;
        CLHEP::Hep3Vector counterPos=recoTimeIter->second.at(0)->_pos;

        //separate hits for each SiPM
        std::map<int,std::vector<CrvHit>::const_iterator> sipmHits;
        const std::vector<std::vector<CrvHit>::const_iterator> &counterHits=recoTimeIter->second;
        for(auto counterHit=counterHits.begin(); counterHit!=counterHits.end(); ++counterHit)
        {
          auto sipmHit = sipmHits.find((*counterHit)->_SiPM);
          if(sipmHit==sipmHits.end()) sipmHits.emplace((*counterHit)->_SiPM,*counterHit);
          else
          {
            //if more than one hit per SiPM, use the hit with the highest PE
            if(sipmHit->second->_PEs<(*counterHit)->_PEs) sipmHit->second=*counterHit;
          }
        }

        //calculate average side times for this counter
        std::array<size_t, CRVId::nSidesPerBar> sideHitsCounter{0,0};
        std::array<float, CRVId::nSidesPerBar>  sidePEsCounter{0,0};
        std::array<double, CRVId::nSidesPerBar> sideTimesCounter{0,0};
        for(auto sipmHit=sipmHits.begin(); sipmHit!=sipmHits.end(); ++sipmHit)
        {
          int side=sipmHit->first%CRVId::nSidesPerBar;
          //use fit values and not the no-fit option
          sideHitsCounter[side]++;
          sidePEsCounter[side]+=sipmHit->second->_crvRecoPulse->GetPEs();
          sideTimesCounter[side]+=sipmHit->second->_crvRecoPulse->GetPulseTime()*sipmHit->second->_crvRecoPulse->GetPEs();  //PE-weighted pulse time average
        }
        int validSides=0;
        for(size_t side=0; side<CRVId::nSidesPerBar; ++side)
        {
          if(sidePEsCounter[side]>0) {sideTimesCounter[side]/=sidePEsCounter[side]; ++validSides;} else sideTimesCounter[side]=0.0;
        }

        //calculate some global numbers
        for(size_t side=0; side<CRVId::nSidesPerBar; ++side)
        {
          sideHits[side]+=sideHitsCounter[side];
          sidePEs[side]+=sidePEsCounter[side];
          sideTimes[side]+=sideTimesCounter[side]*sidePEsCounter[side];
        }
        double halfLength=_sectorMap.at(crvSector).counterHalfLength;
        avgHalfLength+=halfLength*(sidePEsCounter[0]+sidePEsCounter[1]);

        //calculate hit times and hit positions for this counter
        if(validSides==CRVId::nSidesPerBar)  //can calculate the longitudinal position only, if there are hits on both readout sides
        {
          double timeDifference=sideTimesCounter[0]-sideTimesCounter[1];
          double distanceFromCounterCenter=0.5*timeDifference*_fiberSignalSpeed;
          CLHEP::Hep3Vector offsetFromCounterCenter;
          offsetFromCounterCenter[_sectorMap.at(crvSector).lengthDirection]=distanceFromCounterCenter;
          CLHEP::Hep3Vector avgHitPosCounter=counterPos+offsetFromCounterCenter;

          double hitTime0=sideTimesCounter[0] - (halfLength+distanceFromCounterCenter)/_fiberSignalSpeed;
          double hitTime1=sideTimesCounter[1] - (halfLength-distanceFromCounterCenter)/_fiberSignalSpeed;
          double avgHitTimeCounter=0.5*(hitTime0+hitTime1);  //both hit times and the avg hit should be identical

          PEbothSides+=sidePEsCounter[0]+sidePEsCounter[1];
          avgHitPos+=avgHitPosCounter*(sidePEsCounter[0]+sidePEsCounter[1]);
          avgHitTime+=avgHitTimeCounter*(sidePEsCounter[0]+sidePEsCounter[1]);
        }
      }//hit time/pos for one counter

      for(size_t side=0; side<CRVId::nSidesPerBar; ++side)
      {
        if(sidePEs[side]>0) sideTimes[side]/=sidePEs[side]; else sideTimes[side]=0.0;
      }
      bool twoReadoutSides=false;
      if(PEbothSides>0)
      {
        twoReadoutSides=true;
        avgHitPos/=PEbothSides;
        avgHitTime/=PEbothSides;
      }
      else
      {
        //hits only at one readout side.
        //longitudinal position of hit cannot be calculated.
        //use the average counter position instead, i.e. use the center of the counter as longitudinal position (avgCounterPos)
        //and adjust the hit time
        avgHitPos=avgCounterPos;
        if(sidePEs[0]>0)
        {
          avgHalfLength/=sidePEs[0];
          avgHitTime=sideTimes[0]-avgHalfLength/_fiberSignalSpeed;
        }
        else if(sidePEs[1]>0)
        {
          avgHalfLength/=sidePEs[1];
          avgHitTime=sideTimes[1]-avgHalfLength/_fiberSignalSpeed;
        }
      }
      avgHitTime-=_timeOffset;  //remove additional time due to electronics response and processes such as scintillation/WLS decay times

      //insert the cluster information into the vector of the crv coincidence clusters
      crvCoincidenceClusterCollection->emplace_back(crvSectorType, startTime, endTime, PEs, crvRecoPulses, slope, layers,
                                                    sideHits, sidePEs, sideTimes, twoReadoutSides, avgHitTime, avgHitPos);
    } //loop through all clusters
  } //end cluster properies


  //remove hits below the threshold
  void CrvCoincidenceFinder::filterHits(const std::vector<CrvHit> &hits, std::list<CrvHit> &hitsFiltered)
  {
    std::vector<CrvHit>::const_iterator iterHit;
    for(iterHit=hits.begin(); iterHit!=hits.end(); ++iterHit)
    {
      int    layer=iterHit->_layer;
      int    counter=iterHit->_counter;  //counter number in one layer counted from the beginning of the counter type
      double time=iterHit->_time;
      double timePulseStart=iterHit->_timePulseStart;
      double timePulseEnd=iterHit->_timePulseEnd;

      int    PEthreshold=iterHit->_PEthreshold;
      double maxTimeDifferenceAdjacentPulses=iterHit->_maxTimeDifferenceAdjacentPulses;
      double minOverlapTimeAdjacentPulses=iterHit->_minOverlapTimeAdjacentPulses;

      //check other SiPM and the SiPMs at the adjacent counters
      double PEs_thisCounter=0;
      double PEs_adjacentCounter1=0;
      double PEs_adjacentCounter2=0;
      std::vector<CrvHit>::const_iterator iterHitAdjacent;
      for(iterHitAdjacent=hits.begin(); iterHitAdjacent!=hits.end(); ++iterHitAdjacent)
      {
        //use hits of the same layer only
        if(iterHitAdjacent->_layer!=layer) continue;

        //use hits within a certain time window only
        if(!_usePulseOverlaps)
        {
          if(fabs(iterHitAdjacent->_time-time)>maxTimeDifferenceAdjacentPulses) continue;
        }
        else
        {
          double overlapTime=std::min(iterHitAdjacent->_timePulseEnd,timePulseEnd)-std::max(iterHitAdjacent->_timePulseStart,timePulseStart);
          if(overlapTime<minOverlapTimeAdjacentPulses) continue; //no overlap or overlap time too short
        }

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
                                          double clusterMaxTimeDifference, double clusterMinOverlapTime)
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
              if(_usePulseOverlaps)
              {
                if((std::fabs(hitsIter->_x-clusterIter->_x)<=maxDistance) &&
                   (hitsIter->_timePulseEnd-clusterIter->_timePulseStart>clusterMinOverlapTime) &&
                   (clusterIter->_timePulseEnd-hitsIter->_timePulseStart>clusterMinOverlapTime))
                {
                  //this hit satisfied the conditions
                  //move it from list of hits to the current cluster
                  cluster.push_back(*hitsIter);
                  hitsIter=hits.erase(hitsIter);
                  erasedHit=true;
                  break;  //no need for more comparisons with other hits in current cluster, go to the next hit in list of hits
                }
              }
              else
              {
                if((std::fabs(hitsIter->_x-clusterIter->_x)<=maxDistance) &&
                   (std::fabs(hitsIter->_time-clusterIter->_time)<clusterMaxTimeDifference))
                {
                  //this hit satisfied the conditions
                  //move it from list of hits to the current cluster
                  cluster.push_back(*hitsIter);
                  hitsIter=hits.erase(hitsIter);
                  erasedHit=true;
                  break;  //no need for more comparisons with other hits in current cluster, go to the next hit in list of hits
                }
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

    if(!_usePulseOverlaps)
    {
      double maxTimeDifference = (*std::max_element(layerIterators,layerIterators+n,
                                 [](L &a, L &b){return a->_maxTimeDifference < b->_maxTimeDifference;}))->_maxTimeDifference;
      double timeMax = (*std::max_element(layerIterators,layerIterators+n,
                       [](L &a, L &b){return a->_time < b->_time;}))->_time;
      double timeMin = (*std::min_element(layerIterators,layerIterators+n,
                       [](L &a, L &b){return a->_time < b->_time;}))->_time;
      if(timeMax-timeMin>maxTimeDifference) return false;  //hits don't fall within the time window
    }
    else
    {
      double minOverlapTime = (*std::min_element(layerIterators,layerIterators+n,
                                 [](L &a, L &b){return a->_minOverlapTime < b->_minOverlapTime;}))->_minOverlapTime;
      double timeMaxPulseStart = (*std::max_element(layerIterators,layerIterators+n,
                                 [](L &a, L &b){return a->_timePulseStart < b->_timePulseStart;}))->_timePulseStart;
      double timeMinPulseEnd   = (*std::min_element(layerIterators,layerIterators+n,
                                 [](L &a, L &b){return a->_timePulseEnd < b->_timePulseEnd;}))->_timePulseEnd;
      if(timeMinPulseEnd-timeMaxPulseStart<minOverlapTime) return false;  //pulses don't overlap, or overlap time too short
    }

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
