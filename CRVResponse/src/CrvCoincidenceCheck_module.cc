//
// A module to check for coincidences of CRV pulses
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/CRVResponse/inc/CrvHelper.hh"

#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidence.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>
#include <vector>

#include <TMath.h>
#include <TH2D.h>

namespace mu2e
{
  class CrvCoincidenceCheck : public art::EDProducer
  {
    public:
    explicit CrvCoincidenceCheck(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    int                      _verboseLevel;
    std::string              _crvRecoPulsesModuleLabel;
    std::vector<std::string> _CRVSectors;
    std::vector<int>         _PEthresholds;
    std::vector<double>      _adjacentPulseTimeDifferences;
    std::vector<double>      _maxTimeDifferences;
    std::vector<int>         _coincidenceLayers;
    bool        _usePulseOverlaps;
    bool        _useNoFit;
    double      _minOverlapTime;
    bool        _usingPEsPulseHeight;
    double      _maxSlope;
    double      _maxSlopeDifference;
    double      _timeWindowStart;
    double      _timeWindowEnd;

    //the following variable are only used to print out results and summaries
    int         _totalEvents;
    int         _totalEventsCoincidence;
    std::string _moduleLabel;  //for this instance of the CrvCoincidenceCheck module
                               //to distinguish the output from other instances of this module, if there are more than one instances

    struct CrvHit
    {
      art::Ptr<CrvRecoPulse>        _crvRecoPulse;
      double                        _time;
      double                        _timePulseStart, _timePulseEnd;
      float                         _PEs, _PEsNoFit;
      int                           _layer, _counter;
      double                        _x, _y;
      int                           _PEthreshold;
      double                        _adjacentPulseTimeDifference;
      double                        _maxTimeDifference;
      int                           _coincidenceLayers;
      CrvHit(const art::Ptr<CrvRecoPulse> &crvRecoPulse, double time, double timePulseStart, double timePulseEnd, float PEs,
             int layer, int counter, double x, double y, int PEthreshold, double adjacentPulseTimeDifference, double maxTimeDifference,
             int coincidenceLayers):
             _crvRecoPulse(crvRecoPulse), _time(time), _timePulseStart(timePulseStart), _timePulseEnd(timePulseEnd), _PEs(PEs),
             _layer(layer), _counter(counter), _x(x), _y(y), _PEthreshold(PEthreshold), _adjacentPulseTimeDifference(adjacentPulseTimeDifference),
             _maxTimeDifference(maxTimeDifference), _coincidenceLayers(coincidenceLayers) {}
      void Print(int sectorType) const
      {
        std::cout<<"sectorType: "<<sectorType<<"   layer: "<<_layer<<"   counter: "<<_counter<<"  SiPM: "<<_crvRecoPulse->GetSiPMNumber()<<"      ";
        std::cout<<"  PEs: "<<_PEs<<"   time: "<<_time<<"   x: "<<_x<<"   y: "<<_y<<"         "<<_crvRecoPulse->GetScintillatorBarIndex()<<std::endl;
      }
    };

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
      double      adjacentPulseTimeDifference;
      double      maxTimeDifference;
      int         coincidenceLayers;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;

    void checkCombination(const std::vector<std::vector<CrvHit>::const_iterator> &layerIterators,
                          int sectorType, std::unique_ptr<CrvCoincidenceCollection> &crvCoincidenceCollection);
  };

  CrvCoincidenceCheck::CrvCoincidenceCheck(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _verboseLevel(pset.get<int>("verboseLevel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _CRVSectors(pset.get<std::vector<std::string> >("CRVSectors")),
    _PEthresholds(pset.get<std::vector<int> >("PEthresholds")),
    _adjacentPulseTimeDifferences(pset.get<std::vector<double> >("adjacentPulseTimeDifferences")),
    _maxTimeDifferences(pset.get<std::vector<double> >("maxTimeDifferences")),
    _coincidenceLayers(pset.get<std::vector<int> >("coincidenceLayers")),
    _usePulseOverlaps(pset.get<bool>("usePulseOverlaps")),
    _useNoFit(pset.get<bool>("useNoFit")),
    _minOverlapTime(pset.get<double>("minOverlapTime")),
    _usingPEsPulseHeight(pset.get<bool>("usingPEsPulseHeight")),
    _maxSlope(pset.get<double>("maxSlope")),
    _maxSlopeDifference(pset.get<double>("maxSlopeDifference")),
    _timeWindowStart(pset.get<double>("timeWindowStart")),
    _timeWindowEnd(pset.get<double>("timeWindowEnd"))
  {
    produces<CrvCoincidenceCollection>();
    _totalEvents=0;
    _totalEventsCoincidence=0;
  }

  void CrvCoincidenceCheck::beginJob()
  {
  }

  void CrvCoincidenceCheck::endJob()
  {
    if(_verboseLevel>0)
    {
      std::cout<<"SUMMARY "<<_moduleLabel<<"    "<<_totalEventsCoincidence<<" / "<<_totalEvents<<" events satisfied coincidence requirements"<<std::endl;
    }
  }

  void CrvCoincidenceCheck::beginRun(art::Run &run)
  {
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<CRSScintillatorShield> &sectors = CRS->getCRSScintillatorShields();
    for(unsigned int i=0; i<sectors.size(); i++)
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
      s.adjacentPulseTimeDifference = _adjacentPulseTimeDifferences[userPropertyPosition];
      s.maxTimeDifference = _maxTimeDifferences[userPropertyPosition];
      s.coincidenceLayers = _coincidenceLayers[userPropertyPosition];

      _sectorMap[i]=s;
    }
  }

  void CrvCoincidenceCheck::produce(art::Event& event)
  {
    std::unique_ptr<CrvCoincidenceCollection> crvCoincidenceCollection(new CrvCoincidenceCollection);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);

    //collect crvHits
    std::map<int, std::vector<CrvHit> > crvHits;    //hits are separated by sector type (like CRV-T, CRV-R, ...)
                                                    //the key is -sector type for sipms at side 0
                                                    //the key is +sector type for sipms at side 1
                                                    //(sector types start at 1)

    //loop over reco pulse collection (=loop over counters)
    for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
    {
      const art::Ptr<CrvRecoPulse> crvRecoPulse(crvRecoPulseCollection, recoPulseIndex);

      //get information about the counter
      const CRSScintillatorBarIndex &crvBarIndex = crvRecoPulse->GetScintillatorBarIndex();
      int sectorNumber, moduleNumber, layerNumber, barNumber;
      CrvHelper::GetCrvCounterInfo(CRS, crvBarIndex, sectorNumber, moduleNumber, layerNumber, barNumber);

      //need to find the counter number within the entire sector type (like CRV-T, CRV-R, ...)
      std::map<int,sectorCoincidenceProperties>::const_iterator sIter = _sectorMap.find(sectorNumber);
      if(sIter==_sectorMap.end()) throw std::logic_error("CrvCoincidenceFinder: Found a CRV hit at a CRV sector without properties.");
      const sectorCoincidenceProperties &sector = sIter->second;

      int counterNumber = sector.precedingCounters + sector.nCountersPerModule*moduleNumber + barNumber;

      CLHEP::Hep3Vector crvCounterPos=CrvHelper::GetCrvCounterPos(CRS, crvBarIndex);
//      CLHEP::Hep3Vector crvCounterPos;
      double x=crvCounterPos[sector.widthDirection];
      double y=crvCounterPos[sector.thicknessDirection];

      //get the reco pulses information
      int SiPM = crvRecoPulse->GetSiPMNumber();

      //ignore SiPMs on counter sides which don't have SiPMs according to the geometry file
      int counterSide=SiPM%2;
      if(counterSide==0 && !sector.sipmsAtSide0) continue;
      if(counterSide==1 && !sector.sipmsAtSide1) continue;

      //find the hit mapkey (positive numbers for side 1, negative numbers for side 0)
      int sectorType=sector.sectorType;
      if(counterSide==0) sectorType*=-1;

      double time=crvRecoPulse->GetPulseTime();
      double timePulseStart=crvRecoPulse->GetPulseStart();
      double timePulseEnd=crvRecoPulse->GetPulseEnd();
      float  PEs =crvRecoPulse->GetPEs();
      if(_usingPEsPulseHeight) PEs=crvRecoPulse->GetPEsPulseHeight();
      if(_useNoFit) PEs=crvRecoPulse->GetPEsNoFit();
      if(_useNoFit) time=crvRecoPulse->GetPulseTimeNoFit();

      if(_verboseLevel>4)
      {
        std::cout<<"sector: "<<sectorNumber<<"  module: "<<moduleNumber<<"  layer: "<<layerNumber<<"  counter: "<<counterNumber<<"  SiPM: "<<SiPM<<"      ";
        std::cout<<"  PEs: "<<PEs<<"   time: "<<time<<"   x: "<<x<<"   y: "<<y<<"       bar: "<<barNumber<<std::endl;
      }

      //check whether this reco pulses is within the time window
      //(don't check the PE threshold yet, since a hit which doesn't reach the threshold can be combined with a hit of an adjacent counter)
      if(crvRecoPulse->GetPulseTime()>=_timeWindowStart && crvRecoPulse->GetPulseTime()<=_timeWindowEnd)
      {
        //get the right set of hits based on the hitmap key, and insert a new hit
        crvHits[sectorType].emplace_back(crvRecoPulse, time, timePulseStart, timePulseEnd, PEs, layerNumber, counterNumber, x,y,
                                         sector.PEthreshold, sector.adjacentPulseTimeDifference, sector.maxTimeDifference, sector.coincidenceLayers);
        if(_verboseLevel==4) crvHits[sectorType].back().Print(sectorType);
      }//loop over SiPM
    }//loop over reco pulse collection


    //find coincidences for each sector type and side (=hitmap key)
    std::map<int,std::vector<CrvHit> >::const_iterator iterHitMap;
    for(iterHitMap = crvHits.begin(); iterHitMap!=crvHits.end(); iterHitMap++)
    {
      //this is the collection for which a coincidence needs to be found
      int sectorType=iterHitMap->first;
      const std::vector<CrvHit> &crvHitsOfSectorType = iterHitMap->second;

      //remove hits below the threshold
      std::vector< std::vector<CrvHit> > crvHitsFiltered;  //separated by layers
      crvHitsFiltered.resize(4); //for 4 layers
      std::vector<CrvHit>::const_iterator iterHit;
      for(iterHit=crvHitsOfSectorType.begin(); iterHit!=crvHitsOfSectorType.end(); iterHit++)
      {
        int layer=iterHit->_layer;
        int counter=iterHit->_counter;
        float  PEs=iterHit->_PEs;
        double time=iterHit->_time;
        double timePulseStart=iterHit->_timePulseStart;
        double timePulseEnd=iterHit->_timePulseEnd;

        int    PEthreshold=iterHit->_PEthreshold;
        double adjacentPulseTimeDifference=iterHit->_adjacentPulseTimeDifference;

        //check other SiPM and the SiPMs at the adjacent counters
        float PEs_thisCounter=PEs;
        float PEs_adjacentCounter1=0;
        float PEs_adjacentCounter2=0;
        std::vector<CrvHit>::const_iterator iterHitAdjacent;
        for(iterHitAdjacent=crvHitsOfSectorType.begin(); iterHitAdjacent!=crvHitsOfSectorType.end(); iterHitAdjacent++)
        {
          if(std::distance(iterHitAdjacent,iterHit)==0) continue;  //don't compare with itself
          if(iterHitAdjacent->_layer!=layer) continue;             //compare hits of the same layer only

          //compare hits within a certain time window only
          if(!_usePulseOverlaps)
          {
            if(fabs(iterHitAdjacent->_time-time)>adjacentPulseTimeDifference) continue;
          }
          else
          {
            double overlapTime=std::min(iterHitAdjacent->_timePulseEnd,timePulseEnd)-std::max(iterHitAdjacent->_timePulseStart,timePulseStart);
            if(overlapTime<_minOverlapTime) continue; //no overlap or overlap time too short
          }

          int counterDiff=iterHitAdjacent->_counter-counter;
          if(counterDiff==0) PEs_thisCounter+=iterHitAdjacent->_PEs;   //add PEs from the same counter (i.e. the "other" SiPM),
                                                                       //if the "other" hit is within a certain time window (5ns)
          if(counterDiff==-1) PEs_adjacentCounter1+=iterHitAdjacent->_PEs;  //add PEs from an adjacent counter,
                                                                            //if these hits are within a certain time window (5ns)
          if(counterDiff==1) PEs_adjacentCounter2+=iterHitAdjacent->_PEs;   //add PEs from an adjacent counter,
                                                                            //if these hits are within a certain time window (5ns)
        }
        //check, if the number of PEs of the adjacent counter added to the current hit's PE number
        //brings this hit above the threshold, add this hit to vector of hits
        if(PEs_thisCounter+PEs_adjacentCounter1>=PEthreshold) crvHitsFiltered[layer].push_back(*iterHit);
        else {if(PEs_thisCounter+PEs_adjacentCounter2>=PEthreshold) crvHitsFiltered[layer].push_back(*iterHit);}
      }

      //***************************************************
      //find coincidences using 4/4 coincidence requirement
      {
        const std::vector<CrvHit> &layer0Hits=crvHitsFiltered[0];
        const std::vector<CrvHit> &layer1Hits=crvHitsFiltered[1];
        const std::vector<CrvHit> &layer2Hits=crvHitsFiltered[2];
        const std::vector<CrvHit> &layer3Hits=crvHitsFiltered[3];
        std::vector<CrvHit>::const_iterator layer0Iter;
        std::vector<CrvHit>::const_iterator layer1Iter;
        std::vector<CrvHit>::const_iterator layer2Iter;
        std::vector<CrvHit>::const_iterator layer3Iter;

        //it will loop only, if all 4 layers have hits
        for(layer0Iter=layer0Hits.begin(); layer0Iter!=layer0Hits.end(); layer0Iter++)
        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); layer1Iter++)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); layer2Iter++)
        for(layer3Iter=layer3Hits.begin(); layer3Iter!=layer3Hits.end(); layer3Iter++)
        {
          const std::vector<std::vector<CrvHit>::const_iterator> layerIterators={layer0Iter, layer1Iter, layer2Iter, layer3Iter};
          checkCombination(layerIterators, sectorType, crvCoincidenceCollection);
        }
      } // four layer coincidences

      //***************************************************
      //find coincidences using 3/4 coincidence requirement
      for(int layer1=0; layer1<4; layer1++)
      for(int layer2=layer1+1; layer2<4; layer2++)
      for(int layer3=layer2+1; layer3<4; layer3++)
      {

        const std::vector<CrvHit> &layer1Hits=crvHitsFiltered[layer1];
        const std::vector<CrvHit> &layer2Hits=crvHitsFiltered[layer2];
        const std::vector<CrvHit> &layer3Hits=crvHitsFiltered[layer3];
        std::vector<CrvHit>::const_iterator layer1Iter;
        std::vector<CrvHit>::const_iterator layer2Iter;
        std::vector<CrvHit>::const_iterator layer3Iter;

        //it will loop only, if all 3 layers have hits
        //loops will exit, if all three hits require a 4/4 coincidence
        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); layer1Iter++)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); layer2Iter++)
        for(layer3Iter=layer3Hits.begin(); layer3Iter!=layer3Hits.end(); layer3Iter++)
        {
          if(layer1Iter->_coincidenceLayers>3 && layer2Iter->_coincidenceLayers>3 && layer3Iter->_coincidenceLayers>3) continue; //all hits require at 4/4 coincidence

          const std::vector<std::vector<CrvHit>::const_iterator> layerIterators={layer1Iter, layer2Iter, layer3Iter};
          checkCombination(layerIterators, sectorType, crvCoincidenceCollection);
        }
      }  //three layer coincidences

      //***************************************************
      //find coincidences using 2/4 coincidence requirement
      for(int layer1=0; layer1<4; layer1++)
      for(int layer2=layer1+1; layer2<4; layer2++)
      {

        const std::vector<CrvHit> &layer1Hits=crvHitsFiltered[layer1];
        const std::vector<CrvHit> &layer2Hits=crvHitsFiltered[layer2];
        std::vector<CrvHit>::const_iterator layer1Iter;
        std::vector<CrvHit>::const_iterator layer2Iter;

        //it will loop only, if both layers have hits
        //loops will exit, if all three hits require a 3/4 or 4/4 coincidence
        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); layer1Iter++)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); layer2Iter++)
        {
          if(layer1Iter->_coincidenceLayers>2 && layer2Iter->_coincidenceLayers>2) continue; //all hits require at least a 3/4 coincidence

          const std::vector<std::vector<CrvHit>::const_iterator> layerIterators={layer1Iter, layer2Iter};
          checkCombination(layerIterators, sectorType, crvCoincidenceCollection);
        }
      }  //two layer coincidences
    }

    _totalEvents++;
    if(crvCoincidenceCollection->size()>0) _totalEventsCoincidence++;
    _moduleLabel = moduleDescription().moduleLabel();

    if(_verboseLevel>0)
    {
      std::cout<<_moduleLabel<<"   run "<<event.id().run()<<"  subrun "<<event.id().subRun()<<"  event "<<event.id().event()<<"    ";
      std::cout<<(crvCoincidenceCollection->size()>0?"Coincidence satisfied":"No coincidence found")<<std::endl;
    }

/*
if(crvCoincidenceCollection->size()==0)
{
  std::cout<<"============="<<std::endl;
  for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
  {
    const art::Ptr<CrvRecoPulse> crvRecoPulse(crvRecoPulseCollection, recoPulseIndex);
    const CRSScintillatorBarIndex &barIndex = crvRecoPulse->GetScintillatorBarIndex();
    const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);
    const CRSScintillatorBarId &barId = CRSbar.id();
    int sectorNumber=barId.getShieldNumber();
    int moduleNumber=barId.getModuleNumber();
    int layerNumber=barId.getLayerNumber();
    int barNumber=barId.getBarNumber();

    //need to find the counter number within the entire sector type (like CRV-T, CRV-R, ...)
    std::map<int,sectorCoincidenceProperties>::const_iterator sIter = _sectorMap.find(sectorNumber);
    if(sIter==_sectorMap.end()) throw std::logic_error("CrvCoincidenceFinder: Found a CRV hit at a CRV sector without properties.");
    const sectorCoincidenceProperties &sector = sIter->second;

    int counterNumber = sector.precedingCounters + sector.nCountersPerModule*moduleNumber + barNumber;

    double x=CRSbar.getPosition()[sector.widthDirection];
    double y=CRSbar.getPosition()[sector.thicknessDirection];

    //get the reco pulses information
    int SiPM = crvRecoPulse->GetSiPMNumber();

    //ignore SiPMs on counter sides which don't have SiPMs according to the geometry file
    int counterSide=SiPM%2;
    if(counterSide==0 && !sector.sipmsAtSide0) continue;
    if(counterSide==1 && !sector.sipmsAtSide1) continue;

    //find the hit mapkey (positive numbers for side 1, negative numbers for side 0)
    int sectorType=sector.sectorType;
    if(counterSide==0) sectorType*=-1;

    double time=crvRecoPulse->GetPulseTime();
    double LEtime=crvRecoPulse->GetLEtime();
    int    PEs =crvRecoPulse->GetPEs();

    std::cout<<"sector: "<<sectorNumber<<"  module: "<<moduleNumber<<"  layer: "<<layerNumber<<"  counter: "<<counterNumber<<"  SiPM: "<<SiPM<<"      ";
    std::cout<<"  PEs: "<<PEs<<"   time: "<<time<<"   LE: "<<LEtime<<"   x: "<<x<<"   y: "<<y<<"       bar: "<<barNumber<<std::endl;
  }
}
*/

    event.put(std::move(crvCoincidenceCollection));

  } // end produce

  void CrvCoincidenceCheck::checkCombination(const std::vector<std::vector<CrvHit>::const_iterator> &layerIterators,
                                             int sectorType, std::unique_ptr<CrvCoincidenceCollection> &crvCoincidenceCollection)
  {
    size_t n=layerIterators.size();
    if(n==0) return;

    std::vector<double> maxTimeDifferences(n);
    std::vector<double> times(n);
    std::vector<double> timesPulseStart(n), timesPulseEnd(n);
    std::vector<double> x(n), y(n);
    std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses;
    std::vector<int> layers;
    for(size_t i=0; i<n; ++i)
    {
      const std::vector<CrvHit>::const_iterator &iter=layerIterators[i];
      maxTimeDifferences[i]=iter->_maxTimeDifference;
      times[i]=iter->_time;
      timesPulseStart[i]=iter->_timePulseStart;
      timesPulseEnd[i]=iter->_timePulseEnd;
      x[i]=iter->_x;
      y[i]=iter->_y;
      crvRecoPulses.push_back(iter->_crvRecoPulse);
      layers.push_back(iter->_layer);
    }

    if(!_usePulseOverlaps)
    {
      double maxTimeDifference=*std::max_element(maxTimeDifferences.begin(),maxTimeDifferences.end());
      double timeMin = *std::min_element(times.begin(),times.end());
      double timeMax = *std::max_element(times.begin(),times.end());
      if(timeMax-timeMin>maxTimeDifference) return;  //hits don't fall within the time window
    }
    else
    {
      double timeMaxPulseStart = *std::max_element(timesPulseStart.begin(),timesPulseStart.end());
      double timeMinPulseEnd = *std::min_element(timesPulseEnd.begin(),timesPulseEnd.end());
      if(timeMinPulseEnd-timeMaxPulseStart<_minOverlapTime) return;  //pulses don't overlap, or overlap time too short
    }

    std::vector<float> slopes;
    for(size_t d=0; d<n-1; ++d)
    {
      slopes.push_back((x[d+1]-x[d])/(y[d+1]-y[d]));   //width direction / thickness direction
      if(fabs(slopes.back())>_maxSlope) return;  //not more than maxSlope allowed for coincidence;
    }

    if(n>2)
    {
      //slope must not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)
      if(fabs(slopes[0]-slopes[1])>_maxSlopeDifference) return;
      if(n>3)
      {
        if(fabs(slopes[0]-slopes[2])>_maxSlopeDifference) return;
        if(fabs(slopes[1]-slopes[2])>_maxSlopeDifference) return;
      }
    }

    crvCoincidenceCollection->emplace_back(crvRecoPulses, sectorType, slopes, layers);
  } // end checkCombination

} // end namespace mu2e

using mu2e::CrvCoincidenceCheck;
DEFINE_ART_MODULE(CrvCoincidenceCheck)
