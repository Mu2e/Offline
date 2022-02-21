//
// A module to find clusters of coincidences of CRV pulses
//
// 
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/CRVResponse/inc/CrvHelper.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Sequence.h"

#include <string>

namespace mu2e 
{
  class CrvCoincidenceFinder : public art::EDProducer 
  {
    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> verboseLevel{Name("verboseLevel"), Comment("verbose level")};
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("module label of CrvRecoPulses")};
      //cluster settings
      fhicl::Atom<double> maxDistance{Name("maxDistance"), Comment("maximum distance between hits to be considered for a hit cluster")};
      fhicl::Atom<double> maxTimeDifference{Name("maxTimeDifference"), Comment("maximum time difference between hits to be considered for a hit cluster")};
      //sector specific coincidence settings
      fhicl::Sequence<std::string> CRVSectors{Name("CRVSectors"), Comment("CRV sectors")};
      fhicl::Sequence<int> PEthresholds{Name("PEthresholds"), Comment("PE thresholds required for a coincidence")};
      fhicl::Sequence<double> maxTimeDifferencesAdjacentPulses{Name("maxTimeDifferencesAdjacentPulses"), Comment("maximum time differences of pulses of adjacent channels considered for coincidences")};
      fhicl::Sequence<double> maxTimeDifferences{Name("maxTimeDifferences"), Comment("maximum time differences of a coincidence hit combination")};
      fhicl::Sequence<int> coincidenceLayers{Name("coincidenceLayers"), Comment("number of layers required for a coincidence")};
      //coincidence settings for overlap option
      fhicl::Atom<bool> usePulseOverlaps{Name("usePulseOverlaps"), Comment("use pulse overlaps instead of peak times to determine coincidences")};
      fhicl::Atom<double> minOverlapTimeAdjacentPulses{Name("minOverlapTimeAdjacentPulses"), Comment("minimum overlap time between pulses of adjacent channels to be considered for coincidences")};
      fhicl::Atom<double> minOverlapTime{Name("minOverlapTime"), Comment("minimum overlap time between pulses of a coincidence hit combination")};
      //other coincidence settings
      fhicl::Atom<bool> useNoFitReco{Name("useNoFitReco"), Comment("use pulse reco results not based on a Gumbel fit")};
      fhicl::Atom<bool> usePEsPulseHeight{Name("usePEsPulseHeight"), Comment("use PEs determined by pulse height instead of pulse area")};
      fhicl::Atom<double> maxSlope{Name("maxSlope"), Comment("maximum slope allowed for a coincidence")};
      fhicl::Atom<double> maxSlopeDifference{Name("maxSlopeDifference"), Comment("maximum slope differences between layers allowed for a coincidence")};
    };

    typedef art::EDProducer::Table<Config> Parameters;

    explicit CrvCoincidenceFinder(const Parameters& config);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    int         _verboseLevel;
    std::string _crvRecoPulsesModuleLabel;

    double      _maxDistance;
    double      _maxTimeDifference;

    std::vector<std::string> _CRVSectors;
    std::vector<int>         _PEthresholds;
    std::vector<double>      _maxTimeDifferencesAdjacentPulses;
    std::vector<double>      _maxTimeDifferences;
    std::vector<int>         _coincidenceLayers;

    bool        _usePulseOverlaps;
    double      _minOverlapTimeAdjacentPulses;
    double      _minOverlapTime;
    bool        _useNoFitReco;
    bool        _usePEsPulseHeight;
    double      _maxSlope;
    double      _maxSlopeDifference;

    int         _totalEvents;
    int         _totalEventsCoincidence;

    struct sectorCoincidenceProperties
    {
      int  precedingCounters;
      int  nCountersPerModule;
      int  sectorType;
      bool sipmsAtSide0;
      bool sipmsAtSide1;
      int  widthDirection, thicknessDirection;
      std::string name;
      int         PEthreshold;
      double      maxTimeDifferenceAdjacentPulses;
      double      maxTimeDifference;
      int         coincidenceLayers;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;

    struct CrvHit
    {
      art::Ptr<CrvRecoPulse> _crvRecoPulse; 
      CLHEP::Hep3Vector      _pos;
      float  _x, _y;
      double _time, _timePulseStart, _timePulseEnd;
      double _PEs;
      int    _layer;
      int    _counter;
      int    _SiPM;
      int    _PEthreshold;
      double _maxTimeDifferenceAdjacentPulses;
      double _maxTimeDifference;
      int    _coincidenceLayers;
      CrvHit(const art::Ptr<CrvRecoPulse> crvRecoPulse, const CLHEP::Hep3Vector &pos,
             double x, double y, double time, double timePulseStart, double timePulseEnd,
             double PEs, int layer, int counter, int SiPM,
             int PEthreshold, double maxTimeDifferenceAdjacentPulses, double maxTimeDifference, int coincidenceLayers) :
               _crvRecoPulse(crvRecoPulse), _pos(pos),
               _x(x), _y(y), _time(time), _timePulseStart(timePulseStart), _timePulseEnd(timePulseEnd), 
               _PEs(PEs), _layer(layer), _counter(counter), _SiPM(SiPM), 
               _PEthreshold(PEthreshold), _maxTimeDifferenceAdjacentPulses(maxTimeDifferenceAdjacentPulses),
               _maxTimeDifference(maxTimeDifference), _coincidenceLayers(coincidenceLayers) {}
    };

    struct PosCompare
    {
      bool operator() (const CrvHit &lhs, const CrvHit &rhs) const 
      {
        return lhs._x < rhs._x; 
      }
    };

    struct TimeCompare
    {
      bool operator() (const CrvHit &lhs, const CrvHit &rhs) const 
      {
        return lhs._time < rhs._time; 
      }
    };


    bool checkCoincidence(const std::vector<CrvHit> &hits);
    bool checkCombination(const std::vector<CrvHit>::const_iterator layerIterators[], int n);

    constexpr static const int nLayers = 4;
  };

  CrvCoincidenceFinder::CrvCoincidenceFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _verboseLevel(conf().verboseLevel()),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _maxDistance(conf().maxDistance()),
    _maxTimeDifference(conf().maxTimeDifference()),
    _CRVSectors(conf().CRVSectors()),
    _PEthresholds(conf().PEthresholds()),
    _maxTimeDifferencesAdjacentPulses(conf().maxTimeDifferencesAdjacentPulses()),
    _maxTimeDifferences(conf().maxTimeDifferences()),
    _coincidenceLayers(conf().coincidenceLayers()),
    _usePulseOverlaps(conf().usePulseOverlaps()),
    _minOverlapTimeAdjacentPulses(conf().minOverlapTimeAdjacentPulses()),
    _minOverlapTime(conf().minOverlapTime()),
    _useNoFitReco(conf().useNoFitReco()),
    _usePEsPulseHeight(conf().usePEsPulseHeight()),
    _maxSlope(conf().maxSlope()),
    _maxSlopeDifference(conf().maxSlopeDifference()),
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
    for(unsigned int i=0; i<sectors.size(); ++i)
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
      s.name=sectors[i].getName();

      s.sipmsAtSide0=sectors[i].getCRSScintillatorBarDetail().hasCMB(0);
      s.sipmsAtSide1=sectors[i].getCRSScintillatorBarDetail().hasCMB(1);

      s.widthDirection=sectors[i].getCRSScintillatorBarDetail().getWidthDirection();
      s.thicknessDirection=sectors[i].getCRSScintillatorBarDetail().getThicknessDirection();

      std::string sectorName = s.name.substr(4); //removes the "CRV_" part
      std::vector<std::string>::iterator userPropertyIter=std::find(_CRVSectors.begin(), _CRVSectors.end(), sectorName);
      if(userPropertyIter==_CRVSectors.end())
        throw std::logic_error("CrvCoincidenceFinder: The geometry has a CRV sector for which no coincidence properties were defined in the fcl file.");
      int userPropertyPosition = std::distance(_CRVSectors.begin(),userPropertyIter);  //that's the position of the vector in the fcl file which sets PE thresholds, time differences, etc.

      s.PEthreshold = _PEthresholds[userPropertyPosition];
      s.maxTimeDifferenceAdjacentPulses = _maxTimeDifferencesAdjacentPulses[userPropertyPosition];
      s.maxTimeDifference = _maxTimeDifferences[userPropertyPosition];
      s.coincidenceLayers = _coincidenceLayers[userPropertyPosition];

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

    //loop through all reco pulses
    //distribute them into the crv sector types
    //and put them into a (x-pos ordered) multiset
    std::map<int, std::multiset<CrvHit,PosCompare> > sectorTypeMap;
    for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); ++recoPulseIndex)
    {
      const art::Ptr<CrvRecoPulse> crvRecoPulse(crvRecoPulseCollection, recoPulseIndex);

      //get information about the counter
      const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse->GetScintillatorBarIndex();
      int sectorNumber, moduleNumber, layerNumber, barNumber;
      CrvHelper::GetCrvCounterInfo(CRS, crvBarIndex, sectorNumber, moduleNumber, layerNumber, barNumber);

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

      //ignore SiPMs on counter sides which don't have SiPMs according to the geometry file
      int counterSide=SiPM%2;
      if(counterSide==0 && !sector.sipmsAtSide0) continue;
      if(counterSide==1 && !sector.sipmsAtSide1) continue;

      double time=crvRecoPulse->GetPulseTime();
      double timePulseStart=crvRecoPulse->GetPulseStart();
      double timePulseEnd=crvRecoPulse->GetPulseEnd();
      float  PEs =crvRecoPulse->GetPEs();
      if(_usePEsPulseHeight) PEs=crvRecoPulse->GetPEsPulseHeight();
      if(_useNoFitReco) PEs=crvRecoPulse->GetPEsNoFit();
      if(_useNoFitReco) time=crvRecoPulse->GetPulseTimeNoFit();

      //don't split counter sides for the purpose of finding clusters
      sectorTypeMap[sector.sectorType].emplace(crvRecoPulse, crvCounterPos,
                                               x, y, time, timePulseStart, timePulseEnd, PEs,
                                               layerNumber, counterNumber, SiPM,
                                               sector.PEthreshold, sector.maxTimeDifferenceAdjacentPulses,
                                               sector.maxTimeDifference, sector.coincidenceLayers);
    }

    //loop through all crv sectors types
    std::map<int, std::multiset<CrvHit,PosCompare> >::const_iterator sectorTypeMapIter;
    for(sectorTypeMapIter=sectorTypeMap.begin(); sectorTypeMapIter!=sectorTypeMap.end(); ++sectorTypeMapIter)
    {
      //each entry is a x-pos ordered multiset of hits
      int crvSectorType = sectorTypeMapIter->first;
      const std::multiset<CrvHit,PosCompare> &posOrderedHits = sectorTypeMapIter->second; 

      //loop through the multiset of x-pos ordered hits for this particular crv sector type
      std::multiset<CrvHit,PosCompare>::const_iterator hp=posOrderedHits.begin();
      while(hp!=posOrderedHits.end())
      {
        //find clusters of x positions
        //and fill their hits into a time ordered multiset for this cluster
        std::multiset<CrvHit,TimeCompare> posCluster;  //time ordered hits of this x-position cluster
        posCluster.insert(*hp);
        double lastPos=hp->_x;
        while(++hp!=posOrderedHits.end())
        {
          if(hp->_x-lastPos>_maxDistance) break;  //this x-position does not belong to this cluster anymore
          posCluster.insert(*hp);
          lastPos=hp->_x;
        }
        //filled one full x-position cluster

        //loop through the vector of time ordered hits for this x-position cluster
        std::multiset<CrvHit,TimeCompare>::const_iterator ht=posCluster.begin();
        while(ht!=posCluster.end())
        {
          //find time cluster for this x-pos cluster
          std::vector<CrvHit> posAndTimeCluster;   //cluster in position and time
          std::vector<CrvHit> posAndTimeCluster0;  //cluster in position and time for SiPMs at negative end
          std::vector<CrvHit> posAndTimeCluster1;  //cluster in position and time for SiPMs at positive end
          posAndTimeCluster.push_back(*ht);
          if(ht->_SiPM%2==0) posAndTimeCluster0.push_back(*ht); else posAndTimeCluster1.push_back(*ht);
          double lastTime=ht->_time;
          while(++ht!=posCluster.end())
          {
            if(ht->_time-lastTime>_maxTimeDifference) break;
            posAndTimeCluster.push_back(*ht);
            if(ht->_SiPM%2==0) posAndTimeCluster0.push_back(*ht); else posAndTimeCluster1.push_back(*ht);
            lastTime=ht->_time;
          }
          //filled one full position&time cluster

          //check whether this position&time cluster has at least one coincidence
          //if not, try to find the next cluster
          if(!checkCoincidence(posAndTimeCluster0))
          {
            if(!checkCoincidence(posAndTimeCluster1)) continue;
          }

          //calulcate properties of this cluster
          std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses;
          double startTime=posAndTimeCluster.front()._time;
          double endTime=posAndTimeCluster.back()._time;
          double PEs=0;
          CLHEP::Hep3Vector avgCounterPos;  //PE-weighted average position
          int nHits=posAndTimeCluster.size();
          std::set<int> layerSet;
          double sumX =0;
          double sumY =0;
          double sumYY=0;
          double sumXY=0;
          for(size_t i=0; i<posAndTimeCluster.size(); ++i)
          {
            const CrvHit &chit = posAndTimeCluster[i];
            crvRecoPulses.push_back(chit._crvRecoPulse);
            PEs+=chit._PEs;
            avgCounterPos+=chit._pos*chit._PEs;
            layerSet.insert(chit._layer);
            sumX +=chit._x;
            sumY +=chit._y;
            sumYY+=chit._y*chit._y;
            sumXY+=chit._x*chit._y;
          }

          //average counter position (PE weighted), slope, layers
          avgCounterPos/=PEs;
          double slope=(nHits*sumXY-sumX*sumY)/(nHits*sumYY-sumY*sumY);
          std::vector<int> layers(layerSet.begin(), layerSet.end());

          //insert the cluster information into the vector of the crv coincidence clusters
          crvCoincidenceClusterCollection->emplace_back(crvSectorType, avgCounterPos, startTime, endTime, PEs, crvRecoPulses, slope, layers);

        }//loop over all hits in position cluster
      }//loop over all hits in sector
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

  bool CrvCoincidenceFinder::checkCoincidence(const std::vector<CrvHit> &hits)
  {
    //remove hits below the threshold
    std::vector<CrvHit> hitsFiltered[nLayers];  //separated by layers
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
          if(overlapTime<_minOverlapTimeAdjacentPulses) continue; //no overlap or overlap time too short
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
         hitsFiltered[layer].push_back(*iterHit);
    }

    //***************************************************
    //find coincidences using 2/4 coincidence requirement
    {
      std::vector<CrvHit>::const_iterator layerIterators[2];

      for(int layer1=0; layer1<4; ++layer1)
      for(int layer2=layer1+1; layer2<4; ++layer2)
      {
        const std::vector<CrvHit> &layer1Hits=hitsFiltered[layer1];
        const std::vector<CrvHit> &layer2Hits=hitsFiltered[layer2];
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
          //if one hit combination of this cluster satisfies the coincidence criteria, return true
          //no need to check more combinations for this cluster
          if(checkCombination(layerIterators,2)) return true;
        }
      }
    }  //two layer coincidences

    //***************************************************
    //find coincidences using 3/4 coincidence requirement
    {
      std::vector<CrvHit>::const_iterator layerIterators[3];

      for(int layer1=0; layer1<4; ++layer1)
      for(int layer2=layer1+1; layer2<4; ++layer2)
      for(int layer3=layer2+1; layer3<4; ++layer3)
      {
        const std::vector<CrvHit> &layer1Hits=hitsFiltered[layer1];
        const std::vector<CrvHit> &layer2Hits=hitsFiltered[layer2];
        const std::vector<CrvHit> &layer3Hits=hitsFiltered[layer3];
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
          //if one hit combination of this cluster satisfies the coincidence criteria, return true
          //no need to check more combinations for this cluster
          if(checkCombination(layerIterators,3)) return true;
        }
      }
    }  //three layer coincidences

    //***************************************************
    //find coincidences using 4/4 coincidence requirement
    {
      std::vector<CrvHit>::const_iterator layerIterators[4];

      const std::vector<CrvHit> &layer0Hits=hitsFiltered[0];
      const std::vector<CrvHit> &layer1Hits=hitsFiltered[1];
      const std::vector<CrvHit> &layer2Hits=hitsFiltered[2];
      const std::vector<CrvHit> &layer3Hits=hitsFiltered[3];
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
        //if one hit combination of this cluster satisfies the coincidence criteria, return true
        //no need to check more combinations for this cluster
        if(checkCombination(layerIterators,4)) return true;
      }
    } // four layer coincidences

    return false;  //no coincidences found
  } //end check coincidence

  bool CrvCoincidenceFinder::checkCombination(const std::vector<CrvHit>::const_iterator layerIterators[], int n) 
  {
    double maxTimeDifferences[n];
    double times[n];
    double timesPulseStart[n], timesPulseEnd[n];
    double x[n], y[n];
    for(int i=0; i<n; ++i)
    {
      const std::vector<CrvHit>::const_iterator &iter=layerIterators[i];
      maxTimeDifferences[i]=iter->_maxTimeDifference;
      times[i]=iter->_time;
      timesPulseStart[i]=iter->_timePulseStart;
      timesPulseEnd[i]=iter->_timePulseEnd;
      x[i]=iter->_x;
      y[i]=iter->_y;
    }

    if(!_usePulseOverlaps)
    {
      double maxTimeDifference=*std::max_element(maxTimeDifferences,maxTimeDifferences+n);
      double timeMin = *std::min_element(times,times+n);
      double timeMax = *std::max_element(times,times+n);
      if(timeMax-timeMin>maxTimeDifference) return false;  //hits don't fall within the time window
    }
    else
    {
      double timeMaxPulseStart = *std::max_element(timesPulseStart,timesPulseStart+n);
      double timeMinPulseEnd = *std::min_element(timesPulseEnd,timesPulseEnd+n);
      if(timeMinPulseEnd-timeMaxPulseStart<_minOverlapTime) return false;  //pulses don't overlap, or overlap time too short
    }

    std::vector<float> slopes;
    for(int d=0; d<n-1; ++d)
    {
      slopes.push_back((x[d+1]-x[d])/(y[d+1]-y[d]));   //width direction / thickness direction
      if(fabs(slopes.back())>_maxSlope) return false;  //not more than maxSlope allowed for coincidence;
    }

    if(n>2)
    {
      //slope must not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)
      if(fabs(slopes[0]-slopes[1])>_maxSlopeDifference) return false;
      if(n>3)
      {
        if(fabs(slopes[0]-slopes[2])>_maxSlopeDifference) return false;
        if(fabs(slopes[1]-slopes[2])>_maxSlopeDifference) return false;
      }
    }

    return true; //all coincidence criteria for this combination satisfied
  } // end checkCombination

} // end namespace mu2e

using mu2e::CrvCoincidenceFinder;
DEFINE_ART_MODULE(CrvCoincidenceFinder)
