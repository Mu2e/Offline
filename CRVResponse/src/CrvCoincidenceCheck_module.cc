//
// A module to check for coincidences of CRV pulses
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
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"

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
  class CrvCoincidenceCheck : public art::EDProducer 
  {
    public:
    explicit CrvCoincidenceCheck(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    int         _verboseLevel;
    std::string _crvRecoPulsesModuleLabel;
    std::vector<std::string> _CRVSectors;
    std::vector<int>         _PEthresholds;
    std::vector<double>      _adjacentPulseTimeDifferences;
    std::vector<double>      _maxTimeDifferences;
    double      _maxSlope;
    double      _maxSlopeDifference;
    bool        _acceptThreeAdjacentCounters;
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _microBunchPeriod;

    //the following variable are only used to print out results and summaries
    double      _leadingVetoTime;
    double      _trailingVetoTime;
    double      _totalTime;
    double      _totalDeadTime;
    int         _totalEvents;
    int         _totalEventsCoincidence;
    std::string _moduleLabel;  //for this instance of the CrvCoincidenceCheck module
                               //to distinguish the output from other instances of this module, if there are more than one instances

    struct CrvHit
    {
      double                        _time;
      int                           _PEs;
      mu2e::CRSScintillatorBarIndex _barIndex;
      int                           _layer, _counter;
      int                           _SiPM;
      double                        _x, _y;
      int                           _PEthreshold;
      double                        _adjacentPulseTimeDifference;
      double                        _maxTimeDifference;
      CrvHit(double time, int PEs, mu2e::CRSScintillatorBarIndex barIndex, int layer, int counter, int SiPM, double x, double y, 
                                                          int PEthreshold, double adjacentPulseTimeDifference, double maxTimeDifference):
                                                          _time(time), _PEs(PEs), _barIndex(barIndex), 
                                                          _layer(layer), _counter(counter), _SiPM(SiPM), _x(x), _y(y),
                                                          _PEthreshold(PEthreshold), _adjacentPulseTimeDifference(adjacentPulseTimeDifference), 
                                                          _maxTimeDifference(maxTimeDifference) {}
      void Print(int sectorType) const
      {
        std::cout<<"sectorType: "<<sectorType<<"   layer: "<<_layer<<"   counter: "<<_counter<<"  SiPM: "<<_SiPM<<"      ";
        std::cout<<"  PEs: "<<_PEs<<"   LE: "<<_time<<"   x: "<<_x<<"   y: "<<_y<<"         "<<_barIndex<<std::endl;
      }
    };

    void AddCoincidence(std::unique_ptr<CrvCoincidenceCheckResult> &crvCoincidenceCheckResult, int sectorType, const CrvHit &h1, const CrvHit &h2, const CrvHit &h3);

    struct sectorCoincidenceProperties
    {
      int precedingCounters;
      int nCountersPerModule;
      int sectorType;
      bool sipmsAtSide0;
      bool sipmsAtSide1;
      int widthDirection, thicknessDirection;
      std::string name;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;
  };

  CrvCoincidenceCheck::CrvCoincidenceCheck(fhicl::ParameterSet const& pset) :
    _verboseLevel(pset.get<int>("verboseLevel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _CRVSectors(pset.get<std::vector<std::string> >("CRVSectors")),
    _PEthresholds(pset.get<std::vector<int> >("PEthresholds")),
    _adjacentPulseTimeDifferences(pset.get<std::vector<double> >("adjacentPulseTimeDifferences")),
    _maxTimeDifferences(pset.get<std::vector<double> >("maxTimeDifferences")),
    _maxSlope(pset.get<double>("maxSlope")),
    _maxSlopeDifference(pset.get<double>("maxSlopeDifference")),
    _acceptThreeAdjacentCounters(pset.get<bool>("acceptThreeAdjacentCounters")),
    _timeWindowStart(pset.get<double>("timeWindowStart")),
    _timeWindowEnd(pset.get<double>("timeWindowEnd")),
    _leadingVetoTime(pset.get<double>("leadingVetoTime")),
    _trailingVetoTime(pset.get<double>("trailingVetoTime"))
  {
    produces<CrvCoincidenceCheckResult>();
    _totalTime=0;
    _totalDeadTime=0;
    _totalEvents=0;
    _totalEventsCoincidence=0;
  }

  void CrvCoincidenceCheck::beginJob()
  {
  }

  void CrvCoincidenceCheck::endJob()
  {
    std::cout<<"SUMMARY "<<_moduleLabel<<"    "<<_totalEventsCoincidence<<" / "<<_totalEvents<<" events satisfied coincidence requirements"<<std::endl;
  }

  void CrvCoincidenceCheck::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;

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

      _sectorMap[i]=s;
    }
  }

  void CrvCoincidenceCheck::produce(art::Event& event) 
  {
    std::unique_ptr<CrvCoincidenceCheckResult> crvCoincidenceCheckResult(new CrvCoincidenceCheckResult);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    //collect crvHits
    std::map<int, std::vector<CrvHit> > crvHits;    //hits are separated by sector type (like CRV-T, CRV-R, ...)
                                                    //the key is -sector type for sipms at side 0
                                                    //the key is +sector type for sipms at side 1
                                                    //(sector types start at 1)

    //loop over reco pulse collection (=loop over counters)
    for(CrvRecoPulsesCollection::const_iterator iter=crvRecoPulsesCollection->begin(); 
        iter!=crvRecoPulsesCollection->end(); iter++)
    {
      //get information about the counter
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);
      const CRSScintillatorBarId &barId = CRSbar.id();
      int sectorNumber=barId.getShieldNumber();
      int moduleNumber=barId.getModuleNumber();
      int layerNumber=barId.getLayerNumber();
      int barNumber=barId.getBarNumber();   //counter number within this sector, module, and layer

      //need to find the counter number within the entire sector type (like CRV-T, CRV-R, ...)
      std::map<int,sectorCoincidenceProperties>::const_iterator sIter = _sectorMap.find(sectorNumber);
      if(sIter==_sectorMap.end()) throw std::logic_error("CrvCoincidenceFinder: Found a CRV hit at a CRV sector without properties.");
      const sectorCoincidenceProperties &sector = sIter->second;

      //find the sector specific "user properties" for this hit (as defined in the fcl file)
      std::string sectorName = sector.name.substr(4); //removes the "CRV_" part
      std::vector<std::string>::iterator userPropertyIter=std::find(_CRVSectors.begin(), _CRVSectors.end(), sectorName);
      if(userPropertyIter==_CRVSectors.end()) continue; //found a CRV hit at a CRV sector which has no user properties in the fcl file, 
                                                        //which propably means that the user doesn't want to check for coincidences in this sector
      int    userPropertyPosition = std::distance(_CRVSectors.begin(),userPropertyIter);  //that's the position of the vector in the fcl file which sets PE thresholds, time differences, etc.
      int    PEthreshold = _PEthresholds[userPropertyPosition];
      double adjacentPulseTimeDifference = _adjacentPulseTimeDifferences[userPropertyPosition];
      double maxTimeDifference = _maxTimeDifferences[userPropertyPosition];

      int counterNumber = sector.precedingCounters + sector.nCountersPerModule*moduleNumber + barNumber;

      double x=CRSbar.getPosition()[sector.widthDirection];
      double y=CRSbar.getPosition()[sector.thicknessDirection];

      //get the reco pulses for this counter
      const CrvRecoPulses &crvRecoPulses = iter->second;
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        //ignore SiPMs on counter sides which don't have SiPMs according to the geometry file
        int counterSide=SiPM%2;
        if(counterSide==0 && !sector.sipmsAtSide0) continue;
        if(counterSide==1 && !sector.sipmsAtSide1) continue;

        //find the hit mapkey (positive numbers for side 1, negative numbers for side 0)
        int sectorType=sector.sectorType;
        if(counterSide==0) sectorType*=-1;

        //get the reco pulses for one SiPM
        const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &pulseVector = crvRecoPulses.GetRecoPulses(SiPM);
        for(unsigned int i = 0; i<pulseVector.size(); i++) 
        {
          //get one reco pulse
          const CrvRecoPulses::CrvSingleRecoPulse &pulse = pulseVector[i];
          int PEs=pulse._PEs;
//          double time=pulse._peakTime;  //TODO
          double time=pulse._leadingEdge;
          if(_verboseLevel>4)
          {
            std::cout<<"sector: "<<sectorNumber<<"  module: "<<moduleNumber<<"  layer: "<<layerNumber<<"  counter: "<<counterNumber<<"  SiPM: "<<SiPM<<"      ";
            std::cout<<"  PEs: "<<PEs<<"   LE: "<<time<<"   x: "<<x<<"   y: "<<y<<"       bar: "<<barNumber<<std::endl;
          }

          //check whether this reco pulses is within the time window
          //(don't check the PE threshold yet, since a hit which doesn't reach the threshold can be combined with a hit of an adjacent counter)
          if(time>=_timeWindowStart && time<=_timeWindowEnd)
          {
            //get the right set of hits based on the hitmap key, and insert a new hit
            crvHits[sectorType].emplace_back(time,PEs,barIndex,layerNumber,counterNumber,SiPM, x,y, PEthreshold, adjacentPulseTimeDifference, maxTimeDifference);
            if(_verboseLevel==4) crvHits[sectorType].back().Print(sectorType);
          }

          //allow coincidences between consecutive microbunches (if the time window extends over the microbunch end)
          time+=_microBunchPeriod;
          if(time>=_timeWindowStart && time<=_timeWindowEnd)
          {
            //get the right set of hits based on the hitmap key, and insert a new hit
            crvHits[sectorType].emplace_back(time,PEs,barIndex,layerNumber,counterNumber,SiPM, x,y, PEthreshold, adjacentPulseTimeDifference, maxTimeDifference);
            if(_verboseLevel==4) crvHits[sectorType].back().Print(sectorType);
          }
        }
      }//loop over SiPM
    }//loop over reco pulse collection (=loop over counters)


    //find coincidences for each sector type and side (=hitmap key)
    std::map<int,std::vector<CrvHit> >::const_iterator iterHitMap;
    for(iterHitMap = crvHits.begin(); iterHitMap!=crvHits.end(); iterHitMap++)
    {
      //this is the collection for which a coincidence needs to be found
      const std::vector<CrvHit> &crvHitsOfSectorType = iterHitMap->second; 

      //remove hits below the threshold 
      std::vector< std::vector<CrvHit> > crvHitsFiltered;  //separated by layers
      crvHitsFiltered.resize(4); //for 4 layers
      std::vector<CrvHit>::const_iterator iterHit;
      for(iterHit=crvHitsOfSectorType.begin(); iterHit!=crvHitsOfSectorType.end(); iterHit++)
      {
        int layer=iterHit->_layer;
        int counter=iterHit->_counter;
        int PEs=iterHit->_PEs;
        int time=iterHit->_time;

        int    PEthreshold=iterHit->_PEthreshold;
        double adjacentPulseTimeDifference=iterHit->_adjacentPulseTimeDifference;

        //check other SiPM and the SiPMs at the adjacent counters
        int PEs_thisCounter=PEs;
        int PEs_adjacentCounter1=0;
        int PEs_adjacentCounter2=0;
        std::vector<CrvHit>::const_iterator iterHitAdjacent;
        for(iterHitAdjacent=crvHitsOfSectorType.begin(); iterHitAdjacent!=crvHitsOfSectorType.end(); iterHitAdjacent++)
        {
          if(std::distance(iterHitAdjacent,iterHit)==0) continue;  //don't compare with itself
          if(iterHitAdjacent->_layer!=layer) continue;             //compare hits of the same layer only
          if(fabs(iterHitAdjacent->_time-time)>adjacentPulseTimeDifference) continue; //compare hits within a certain time window only

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

      //find coincidences using 3 hits in 3 layers
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

        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); layer1Iter++)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); layer2Iter++)
        for(layer3Iter=layer3Hits.begin(); layer3Iter!=layer3Hits.end(); layer3Iter++)
        {
          double maxTimeDifferences[3]={layer1Iter->_maxTimeDifference,layer2Iter->_maxTimeDifference,layer3Iter->_maxTimeDifference};
          double maxTimeDifference=*std::max_element(maxTimeDifferences,maxTimeDifferences+3); 

          if(fabs(layer1Iter->_time-layer2Iter->_time)>maxTimeDifference) break;  //no need to check any triplets containing the current pair of layer1 and layer2

          double times[3]={layer1Iter->_time,layer2Iter->_time,layer3Iter->_time};
          double timeMin = *std::min_element(times,times+3);
          double timeMax = *std::max_element(times,times+3);
          if(timeMax-timeMin>maxTimeDifference) continue;  //hits don't fall within the time window
      
          double x[3]={layer1Iter->_x,layer2Iter->_x,layer3Iter->_x};
          double y[3]={layer1Iter->_y,layer2Iter->_y,layer3Iter->_y};

          bool coincidenceFound=true;
          double slope[2];
          for(int d=0; d<2; d++)
          {
            slope[d]=(x[d+1]-x[d])/(y[d+1]-y[d]);
            if(fabs(slope[d])>_maxSlope) coincidenceFound=false;   //not more than maxSlope allowed for coincidence;
          }

          if(fabs(slope[0])>_maxSlope) break;  //no need to check any triplets containing the current pair of layer1 and layer2

          if(fabs(slope[0]-slope[1])>_maxSlopeDifference) coincidenceFound=false;   //slope most not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)

          if(coincidenceFound) AddCoincidence(crvCoincidenceCheckResult,iterHitMap->first,*layer1Iter,*layer2Iter,*layer3Iter);
        }
      }

      //find coincidences using 3 hits in adjacent counters in one layer
      if(_acceptThreeAdjacentCounters)
      {
        for(int layer=0; layer<4; layer++) 
        {
          const std::vector<CrvHit> &layerHits=crvHitsFiltered[layer];
          if(layerHits.size()<3) continue; //less than three hits in this layer

          std::vector<CrvHit>::const_iterator i1;
          std::vector<CrvHit>::const_iterator i2;
          std::vector<CrvHit>::const_iterator i3;
          for(i1=layerHits.begin(); i1!=layerHits.end(); i1++)
          for(i2=i1+1; i2!=layerHits.end(); i2++)
          for(i3=i2+1; i3!=layerHits.end(); i3++)
          {
            double times[3]={i1->_time,i2->_time,i3->_time};
            double timeMin = *std::min_element(times,times+3);
            double timeMax = *std::max_element(times,times+3);

            double maxTimeDifferences[3]={i1->_maxTimeDifference,i2->_maxTimeDifference,i3->_maxTimeDifference};
            double maxTimeDifference=*std::max_element(maxTimeDifferences,maxTimeDifferences+3); 

            if(timeMax-timeMin>maxTimeDifference) continue;  //hits don't fall within the time window
      
            std::set<int> counters{i1->_counter,i2->_counter,i3->_counter};
            bool coincidenceFound=true;
            if(counters.size()<3) coincidenceFound=false;
            if(*counters.rbegin()-*counters.begin()!=2) coincidenceFound=false;

            if(coincidenceFound) AddCoincidence(crvCoincidenceCheckResult,iterHitMap->first,*i1,*i2,*i3);
          }
        }
      }
    }

    _totalEvents++;
    if(crvCoincidenceCheckResult->CoincidenceFound()) _totalEventsCoincidence++;
    _moduleLabel = *this->currentContext()->moduleLabel();

/*************************************************/
//This section is used only to print out results

    if(_verboseLevel>0)
    {
      std::cout<<_moduleLabel<<"   run "<<event.id().run()<<"  subrun "<<event.id().subRun()<<"  event "<<event.id().event()<<"    ";
      std::cout<<(crvCoincidenceCheckResult->CoincidenceFound()?"Coincidence satisfied":"No coincidence found")<<std::endl;

      std::vector<CrvCoincidenceCheckResult::DeadTimeWindow> deadTimeWindows;
      deadTimeWindows = crvCoincidenceCheckResult->GetDeadTimeWindows(_leadingVetoTime,_trailingVetoTime);

      double deadTime = 0;
      for(unsigned int i=0; i < deadTimeWindows.size(); i++)
      {
        double t1 = deadTimeWindows[i]._startTime;
        double t2 = deadTimeWindows[i]._endTime;
        if(t1<_timeWindowStart) t1=_timeWindowStart;
        if(t2<_timeWindowStart) continue;
        if(t1>_timeWindowEnd)   continue;
        if(t2>_timeWindowEnd)   t2=_timeWindowEnd;
        deadTime = t2 - t1;
        std::cout << "   Found Dead time: " << deadTime << " (" << deadTimeWindows[i]._startTime << " ... " << deadTimeWindows[i]._endTime << ")" << std::endl;
        _totalDeadTime += deadTime;
        if(_verboseLevel>1)
        {
          const std::vector<CrvCoincidenceCheckResult::CoincidenceHit> &hits = deadTimeWindows[i]._hits;
          for(unsigned int j=0; j<hits.size(); j++)
          {
            std::cout<<"time: "<<hits[j]._time<<"   PEs: "<<hits[j]._PEs<<"   bar index: "<<hits[j]._counter<<"   SiPM: "<<hits[j]._SiPM<<std::endl;
          }
        }
      }
      _totalTime += _timeWindowEnd - _timeWindowStart;
      double fractionDeadTime = _totalDeadTime / _totalTime;
      std::cout << "Dead time so far: " << _totalDeadTime << " / " << _totalTime << " = " << fractionDeadTime*100 << "%" << std::endl;
    }

/*************************************************/

    event.put(std::move(crvCoincidenceCheckResult));

  } // end produce


  void CrvCoincidenceCheck::AddCoincidence(std::unique_ptr<CrvCoincidenceCheckResult> &crvCoincidenceCheckResult, int sectorType, const CrvHit &h1, const CrvHit &h2, const CrvHit &h3)
  {
    std::vector<const CrvHit*> hits{&h1,&h2,&h3};
    CrvCoincidenceCheckResult::CoincidenceCombination combination;
    for(int k=0; k<3; k++) 
    {
      combination._time[k] = hits[k]->_time;
      combination._PEs[k]  = hits[k]->_PEs;
      combination._counters[k] = hits[k]->_barIndex;
      combination._SiPMs[k] = hits[k]->_SiPM;
      if(_verboseLevel>2) hits[k]->Print(sectorType);
    }
    crvCoincidenceCheckResult->GetCoincidenceCombinations().push_back(combination);
  }

} // end namespace mu2e

using mu2e::CrvCoincidenceCheck;
DEFINE_ART_MODULE(CrvCoincidenceCheck)
