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

#include "art/Persistency/Common/Ptr.h"
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
    int         _PEthreshold;
    double      _adjacentPulseTimeDifference;
    double      _maxTimeDifference;
    double      _maxSlope;
    double      _maxSlopeDifference;
    bool        _acceptThreeAdjacentCounters;
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _microBunchPeriod;
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
      CrvHit(double time, int PEs, mu2e::CRSScintillatorBarIndex barIndex, int layer, int counter, int SiPM, double x, double y): 
                                                          _time(time), _PEs(PEs), _barIndex(barIndex), 
                                                          _layer(layer), _counter(counter), _SiPM(SiPM), _x(x), _y(y) {}
    };

    void AddCoincidence(std::unique_ptr<CrvCoincidenceCheckResult> &crvCoincidenceCheckResult, int sectorType, const CrvHit &h1, const CrvHit &h2, const CrvHit &h3);

    std::string                   _coincidencePropertiesFile;
    int                           _nSectors;
    std::vector<std::string>      _crvSectorNames;

    struct sectorCoincidenceProperties
    {
      int precedingCounters;
      int nCountersPerModule;
      int sectorType;
      bool sipmsAtSide0;
      bool sipmsAtSide1;
      int widthDirection, thicknessDirection;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;
  };

  CrvCoincidenceCheck::CrvCoincidenceCheck(fhicl::ParameterSet const& pset) :
    _verboseLevel(pset.get<int>("verboseLevel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _PEthreshold(pset.get<int>("PEthreshold")),
    _adjacentPulseTimeDifference(pset.get<double>("adjacentPulseTimeDifference")),
    _maxTimeDifference(pset.get<double>("maxTimeDifference")),
    _maxSlope(pset.get<double>("maxSlope")),
    _maxSlopeDifference(pset.get<double>("maxSlopeDifference")),
    _acceptThreeAdjacentCounters(pset.get<bool>("acceptThreeAdjacentCounters")),
    _timeWindowStart(pset.get<double>("timeWindowStart")),
    _timeWindowEnd(pset.get<double>("timeWindowEnd")),
    _coincidencePropertiesFile(pset.get<std::string>("coincidencePropertiesFile"))
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

    SimpleConfig config(_coincidencePropertiesFile, true, false, false);
    _nSectors = config.getInt("crs.nSectors");
    config.getVectorString("crs.sectorNames",_crvSectorNames,_nSectors);
    for(int i=0; i<_nSectors; i++) 
    {
      sectorCoincidenceProperties s;
      s.precedingCounters=0;
      int precedingSector=config.getInt("crs.precedingSectorFor"+_crvSectorNames[i]);
      while(precedingSector!=-1)
      {
        int nModules=config.getInt("crs.nModules"+_crvSectorNames[precedingSector]);
        int nCountersPerModule=config.getInt("crs.nCountersPerModule"+_crvSectorNames[precedingSector]);
        s.precedingCounters+=nModules*nCountersPerModule;
        precedingSector=config.getInt("crs.precedingSectorFor"+_crvSectorNames[precedingSector]);
      }

      s.nCountersPerModule=config.getInt("crs.nCountersPerModule"+_crvSectorNames[i]);

      s.sectorType=config.getInt("crs.sectorType"+_crvSectorNames[i]);
      if(s.sectorType==0) throw std::logic_error("CrvCoincidenceFinder: Sector type cannot be 0.");

      s.sipmsAtSide0=config.getBool("crs.sipmsAtSide0"+_crvSectorNames[i]);
      s.sipmsAtSide1=config.getBool("crs.sipmsAtSide1"+_crvSectorNames[i]);

      CLHEP::Hep3Vector gapDirection   = config.getHep3Vector("crs.gapDirection"+_crvSectorNames[i]);
      CLHEP::Hep3Vector layerDirection = config.getHep3Vector("crs.layerDirection"+_crvSectorNames[i]);
      for(int D=0; D<3; D++)
      {
        if(gapDirection[D]!=0) s.widthDirection=D;
        if(layerDirection[D]!=0) s.thicknessDirection=D;
      }

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
          double time=pulse._leadingEdge;
          if(_verboseLevel>4)
          {
            std::cout<<"sector: "<<sectorNumber<<"  module: "<<moduleNumber<<"  layer: "<<layerNumber<<"  bar: "<<barNumber<<"  SiPM: "<<SiPM<<"      ";
            std::cout<<"  PEs: "<<PEs<<"   LE: "<<time<<"   x: "<<x<<"   y: "<<y<<std::endl;
          }

          //check whether this reco pulses is within the time window
          //(don't check the PE threshold yet, since a hit which doesn't reach the threshold can be combined with a hit of an adjacent counter)
          if(time>=_timeWindowStart && time<=_timeWindowEnd)
          {
            //get the right set of hits based on the hitmap key, and insert a new hit
            crvHits[sectorType].emplace_back(time,PEs,barIndex,layerNumber,counterNumber,SiPM, x,y);
            if(_verboseLevel==4)
            {
              std::cout<<"sectorType: "<<sectorType<<"   layer: "<<layerNumber<<"   counter: "<<counterNumber<<"  SiPM: "<<SiPM<<"      ";
              std::cout<<"  PEs: "<<PEs<<"   LE: "<<time<<"   x: "<<x<<"   y: "<<y<<std::endl;
            }
          }

          //allow coincidences between consecutive microbunches (if the time window extends over the microbunch end)
          time+=_microBunchPeriod;
          if(time>=_timeWindowStart && time<=_timeWindowEnd)
          {
            //get the right set of hits based on the hitmap key, and insert a new hit
            crvHits[sectorType].emplace_back(time,PEs,barIndex,layerNumber,counterNumber,SiPM, x,y);
            if(_verboseLevel==4)
            {
              std::cout<<"sectorType: "<<sectorType<<"   layer: "<<layerNumber<<"   counter: "<<counterNumber<<"  SiPM: "<<SiPM<<"      ";
              std::cout<<"  PEs: "<<PEs<<"   LE: "<<time<<"   x: "<<x<<"   y: "<<y<<std::endl;
            }
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
        if(PEs>=_PEthreshold) crvHitsFiltered[layer].push_back(*iterHit);
        else  //check adjacent counters
        {
          std::vector<CrvHit>::const_iterator iterHitAdjacent;
          for(iterHitAdjacent=crvHitsOfSectorType.begin(); iterHitAdjacent!=crvHitsOfSectorType.end(); iterHitAdjacent++)
          {
            if(iterHitAdjacent->_layer==layer && 
               abs(iterHitAdjacent->_counter-counter)<=1 &&
               fabs(iterHitAdjacent->_time-time)<=_adjacentPulseTimeDifference &&
               iterHitAdjacent->_PEs<=PEs) //another hit at an adjacent counter within a certain window 
                                           //with the same number of PEs or less
            {
              //check, if the number of PEs of the adjacent counter added to the current hit's PE number 
              //brings this hit above the threshold, add this hit to vector of hits
              if(PEs+iterHitAdjacent->_PEs>=_PEthreshold) crvHitsFiltered[layer].push_back(*iterHit);
            }
          }
        }
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
          if(fabs(layer1Iter->_time-layer2Iter->_time)>_maxTimeDifference) break;  //no need to check any triplets containing the current pair of layer1 and layer2

          double times[3]={layer1Iter->_time,layer2Iter->_time,layer3Iter->_time};
          double timeMin = *std::min_element(times,times+3);
          double timeMax = *std::max_element(times,times+3);
          if(timeMax-timeMin>_maxTimeDifference) continue;  //hits don't fall within the time window
      
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
            if(timeMax-timeMin>_maxTimeDifference) continue;  //hits don't fall within the time window
      
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

    if(_verboseLevel>0)
    {
      std::cout<<_moduleLabel<<"   run "<<event.id().run()<<"  subrun "<<event.id().subRun()<<"  event "<<event.id().event()<<"    ";
      std::cout<<(crvCoincidenceCheckResult->CoincidenceFound()?"Coincidence satisfied":"No coincidence found")<<std::endl;

      std::vector<CrvCoincidenceCheckResult::DeadTimeWindow> deadTimeWindows;
      deadTimeWindows = crvCoincidenceCheckResult->GetDeadTimeWindows(25,125);  //TODO: Don't hardcode these numbers

      double deadTime = 0;
      for(unsigned int i=0; i < deadTimeWindows.size(); i++)
      {
        deadTime = deadTimeWindows[i]._endTime - deadTimeWindows[i]._startTime;
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
      _totalTime += _microBunchPeriod;
      double fractionDeadTime = _totalDeadTime / _totalTime;
      std::cout << "Dead time so far: " << _totalDeadTime << " / " << _totalTime << " = " << fractionDeadTime*100 << "%" << std::endl;
    }

    event.put(std::move(crvCoincidenceCheckResult));

  } // end produce


  void CrvCoincidenceCheck::AddCoincidence(std::unique_ptr<CrvCoincidenceCheckResult> &crvCoincidenceCheckResult, int sectorType, const CrvHit &h1, const CrvHit &h2, const CrvHit &h3)
  {
    std::vector<const CrvHit*> hits{&h1,&h2,&h3};
    CrvCoincidenceCheckResult::CoincidenceCombination combination;
    if(_verboseLevel>2) std::cout<<"Coincidence sectorType/layers/counters/times/counters/SiPMs/PEs: "<<std::endl;
    for(int k=0; k<3; k++) 
    {
      combination._time[k] = hits[k]->_time;
      combination._PEs[k]  = hits[k]->_PEs;
      combination._counters[k] = hits[k]->_barIndex;
      combination._SiPMs[k] = hits[k]->_SiPM;
      if(_verboseLevel>2)
      {
        std::cout<<"   "<<sectorType<<" / "<<hits[k]->_layer<<" / "<<hits[k]->_counter<<" / "<<hits[k]->_SiPM
                 <<" / "<<hits[k]->_PEs<<" / "<<hits[k]->_time<<std::endl;
      }
    }
    crvCoincidenceCheckResult->GetCoincidenceCombinations().push_back(combination);
  }

} // end namespace mu2e

using mu2e::CrvCoincidenceCheck;
DEFINE_ART_MODULE(CrvCoincidenceCheck)
