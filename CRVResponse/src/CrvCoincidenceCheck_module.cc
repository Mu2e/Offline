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
#include "CRVResponse/inc/CrvHelper.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCollection.hh"

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
    int                      _verboseLevel;
    std::string              _crvRecoPulsesModuleLabel;
    std::vector<std::string> _CRVSectors;
    std::vector<int>         _PEthresholds;
    std::vector<double>      _adjacentPulseTimeDifferences;
    std::vector<double>      _maxTimeDifferences;
    std::vector<bool>        _useFourLayers;
    bool        _usingPEsPulseHeight;
    double      _maxSlope;
    double      _maxSlopeDifference;
    bool        _acceptThreeAdjacentCounters;
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _microBunchPeriod;

    //the following variable are only used to print out results and summaries
    int         _totalEvents;
    int         _totalEventsCoincidence;
    std::string _moduleLabel;  //for this instance of the CrvCoincidenceCheck module
                               //to distinguish the output from other instances of this module, if there are more than one instances

    // these variables are used for efficiency checks with overlayed background, where the coincidence is accepted only, if it happens around the time of the muon
    bool        _muonsOnly;
    double      _muonMinTime, _muonMaxTime;
    std::string _genParticleModuleLabel;

    struct CrvHit
    {
      art::Ptr<CrvRecoPulse>        _crvRecoPulse;
      double                        _time;
      int                           _PEs;
      int                           _layer, _counter;
      double                        _x, _y;
      int                           _PEthreshold;
      double                        _adjacentPulseTimeDifference;
      double                        _maxTimeDifference;
      bool                          _useFourLayers;
      CrvHit(const art::Ptr<CrvRecoPulse> &crvRecoPulse, double time, int PEs, int layer, int counter, double x, double y,
                                                          int PEthreshold, double adjacentPulseTimeDifference, double maxTimeDifference, bool useFourLayers):
                                                          _crvRecoPulse(crvRecoPulse), _time(time), _PEs(PEs), _layer(layer), _counter(counter), _x(x), _y(y),
                                                          _PEthreshold(PEthreshold), _adjacentPulseTimeDifference(adjacentPulseTimeDifference),
                                                          _maxTimeDifference(maxTimeDifference), _useFourLayers(useFourLayers) {}
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
      bool        useFourLayers;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;
  };

  CrvCoincidenceCheck::CrvCoincidenceCheck(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _verboseLevel(pset.get<int>("verboseLevel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _CRVSectors(pset.get<std::vector<std::string> >("CRVSectors")),
    _PEthresholds(pset.get<std::vector<int> >("PEthresholds")),
    _adjacentPulseTimeDifferences(pset.get<std::vector<double> >("adjacentPulseTimeDifferences")),
    _maxTimeDifferences(pset.get<std::vector<double> >("maxTimeDifferences")),
    _useFourLayers(pset.get<std::vector<bool> >("useFourLayers")),
    _usingPEsPulseHeight(pset.get<bool>("usingPEsPulseHeight")),
    _maxSlope(pset.get<double>("maxSlope")),
    _maxSlopeDifference(pset.get<double>("maxSlopeDifference")),
    _acceptThreeAdjacentCounters(pset.get<bool>("acceptThreeAdjacentCounters")),
    _timeWindowStart(pset.get<double>("timeWindowStart")),
    _timeWindowEnd(pset.get<double>("timeWindowEnd")),
    _muonsOnly(pset.get<bool>("muonsOnly",false))
  {
    produces<CrvCoincidenceCollection>();
    _totalEvents=0;
    _totalEventsCoincidence=0;
    if(_muonsOnly)
    {
      _muonMinTime=pset.get<double>("muonMinTime");
      _muonMaxTime=pset.get<double>("muonMaxTime");
      _genParticleModuleLabel=pset.get<std::string>("genParticleModuleLabel");
    }
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

      std::string sectorName = s.name.substr(4); //removes the "CRV_" part
      std::vector<std::string>::iterator userPropertyIter=std::find(_CRVSectors.begin(), _CRVSectors.end(), sectorName);
      if(userPropertyIter==_CRVSectors.end())
        throw std::logic_error("CrvCoincidenceFinder: The geometry has a CRV sector for which no coincidence properties were defined in the fcl file.");
      int userPropertyPosition = std::distance(_CRVSectors.begin(),userPropertyIter);  //that's the position of the vector in the fcl file which sets PE thresholds, time differences, etc.

      s.PEthreshold = _PEthresholds[userPropertyPosition];
      s.adjacentPulseTimeDifference = _adjacentPulseTimeDifferences[userPropertyPosition];
      s.maxTimeDifference = _maxTimeDifferences[userPropertyPosition];
      s.useFourLayers = _useFourLayers[userPropertyPosition];

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
      int    PEs =crvRecoPulse->GetPEs();
      if(_usingPEsPulseHeight) PEs=crvRecoPulse->GetPEsPulseHeight();

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
        crvHits[sectorType].emplace_back(crvRecoPulse, time, PEs, layerNumber, counterNumber, x,y,
                                         sector.PEthreshold, sector.adjacentPulseTimeDifference, sector.maxTimeDifference, sector.useFourLayers);
        if(_verboseLevel==4) crvHits[sectorType].back().Print(sectorType);
      }//loop over SiPM
    }//loop over reco pulse collection


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

      {
      //find coincidences using 4 hits in 4 layers
        const std::vector<CrvHit> &layer0Hits=crvHitsFiltered[0];
        const std::vector<CrvHit> &layer1Hits=crvHitsFiltered[1];
        const std::vector<CrvHit> &layer2Hits=crvHitsFiltered[2];
        const std::vector<CrvHit> &layer3Hits=crvHitsFiltered[3];
        std::vector<CrvHit>::const_iterator layer0Iter;
        std::vector<CrvHit>::const_iterator layer1Iter;
        std::vector<CrvHit>::const_iterator layer2Iter;
        std::vector<CrvHit>::const_iterator layer3Iter;

        for(layer0Iter=layer0Hits.begin(); layer0Iter!=layer0Hits.end(); layer0Iter++)
        for(layer1Iter=layer1Hits.begin(); layer1Iter!=layer1Hits.end(); layer1Iter++)
        for(layer2Iter=layer2Hits.begin(); layer2Iter!=layer2Hits.end(); layer2Iter++)
        for(layer3Iter=layer3Hits.begin(); layer3Iter!=layer3Hits.end(); layer3Iter++)
        {
          double maxTimeDifferences[4]={layer0Iter->_maxTimeDifference,layer1Iter->_maxTimeDifference,layer2Iter->_maxTimeDifference,layer3Iter->_maxTimeDifference};
          double maxTimeDifference=*std::max_element(maxTimeDifferences,maxTimeDifferences+4);

          double times[4]={layer0Iter->_time,layer1Iter->_time,layer2Iter->_time,layer3Iter->_time};
          double timeMin = *std::min_element(times,times+4);
          double timeMax = *std::max_element(times,times+4);
          if(timeMax-timeMin>maxTimeDifference) continue;  //hits don't fall within the time window

          double x[4]={layer0Iter->_x,layer1Iter->_x,layer2Iter->_x,layer3Iter->_x};
          double y[4]={layer0Iter->_y,layer1Iter->_y,layer2Iter->_y,layer3Iter->_y};

          bool coincidenceFound=true;
          double slope[3];
          for(int d=0; d<3; d++)
          {
            slope[d]=(x[d+1]-x[d])/(y[d+1]-y[d]);
            if(fabs(slope[d])>_maxSlope) coincidenceFound=false;   //not more than maxSlope allowed for coincidence;
          }

          if(fabs(slope[0]-slope[1])>_maxSlopeDifference) coincidenceFound=false;   //slope must not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)
          if(fabs(slope[0]-slope[2])>_maxSlopeDifference) coincidenceFound=false;   //slope must not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)
          if(fabs(slope[1]-slope[2])>_maxSlopeDifference) coincidenceFound=false;   //slope must not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)

          if(_muonsOnly)   //used for efficiency checks with overlayed background: accept coincidence only, if it happens within e.g. 20ns and 120ns
          {
            art::Handle<GenParticleCollection> genParticleCollection;
            event.getByLabel(_genParticleModuleLabel,"",genParticleCollection);
            double genTime = genParticleCollection->at(0).time();
            if(timeMax>genTime+_muonMaxTime || timeMin<genTime+_muonMinTime) coincidenceFound=false;
          }

          if(coincidenceFound)
          {
            std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses{layer0Iter->_crvRecoPulse,layer1Iter->_crvRecoPulse,layer2Iter->_crvRecoPulse,layer3Iter->_crvRecoPulse};
            int sectorType=iterHitMap->first;
            crvCoincidenceCollection->emplace_back(crvRecoPulses, sectorType);
          }
        }
      } // four layer coincidences

      //find coincidences using 3 hits in 3 layers (ignored, if all three hits have a useFourLayers flag)
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
          if(layer1Iter->_useFourLayers && layer2Iter->_useFourLayers && layer3Iter->_useFourLayers) continue; //all hits require a four layer coincidence

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

          if(fabs(slope[0]-slope[1])>_maxSlopeDifference) coincidenceFound=false;   //slope must not change more than 2mm over 1mm (which is a little bit more than 1 counter per layer)

          if(_muonsOnly)   //used for efficiency checks with overlayed background: accept coincidence only, if it happens within e.g. 20ns and 120ns
          {
            art::Handle<GenParticleCollection> genParticleCollection;
            event.getByLabel(_genParticleModuleLabel,"",genParticleCollection);
            double genTime = genParticleCollection->at(0).time();
            if(timeMax>genTime+_muonMaxTime || timeMin<genTime+_muonMinTime) coincidenceFound=false;
          }

          if(coincidenceFound)
          {
            std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses{layer1Iter->_crvRecoPulse,layer2Iter->_crvRecoPulse,layer3Iter->_crvRecoPulse};
            int sectorType=iterHitMap->first;
            crvCoincidenceCollection->emplace_back(crvRecoPulses, sectorType);
          }
        }
      }  //three layer coincidences


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

            if(_muonsOnly)   //used for efficiency checks with overlayed background: accept coincidence only, if it happens within e.g. 20ns and 120ns
            {
              art::Handle<GenParticleCollection> genParticleCollection;
              event.getByLabel(_genParticleModuleLabel,"",genParticleCollection);
              double genTime = genParticleCollection->at(0).time();
              if(timeMax>genTime+_muonMaxTime || timeMin<genTime+_muonMinTime) coincidenceFound=false;
            }

            if(coincidenceFound)
            {
              std::vector<art::Ptr<CrvRecoPulse> > crvRecoPulses{i1->_crvRecoPulse,i2->_crvRecoPulse,i3->_crvRecoPulse};
              int sectorType=iterHitMap->first;
              crvCoincidenceCollection->emplace_back(crvRecoPulses, sectorType);
            }
          }
        }
      } //accept three adjacent counters

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

} // end namespace mu2e

using mu2e::CrvCoincidenceCheck;
DEFINE_ART_MODULE(CrvCoincidenceCheck)
