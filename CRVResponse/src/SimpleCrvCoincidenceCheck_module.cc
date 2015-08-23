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
  class SimpleCrvCoincidenceCheck : public art::EDProducer 
  {
    public:
    explicit SimpleCrvCoincidenceCheck(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    bool        _checkLayers;
    int         _verboseLevel;
    std::string _crvRecoPulsesModuleLabel;
    int         _PEthreshold;
    double      _maxDistance;
    double      _maxTimeDifference;
    double      _timeWindowStart;
    double      _timeWindowEnd;
    double      _microBunchPeriod;
    double      _totalTime;
    double      _totalDeadTime;
    int         _totalEvents;
    int         _totalEventsCoincidence;
    std::string _moduleLabel;  //for this instance of the SimpleCrvCoincidenceCheck module
                               //to distinguish the output from other instances of this module, if there are more than one instances

    struct coincidenceStruct
    {
      double                  pos;
      CRSScintillatorBarIndex counter;
      int                     layer, moduleType;
      std::vector<double>     time;
      std::vector<int>        PEs;
      std::vector<int>        SiPMs;
    };

    std::string                   _coincidenceGroupsFile;
    int                           _nSectors;
    std::vector<std::string>      _crvSectorNames;

    struct sectorCoincidenceProperties
    {
      int coincidenceGroup[2];
      int moduleType;
      int layerDirection;
    };
    std::map<int,sectorCoincidenceProperties> _sectorMap;
  };

  SimpleCrvCoincidenceCheck::SimpleCrvCoincidenceCheck(fhicl::ParameterSet const& pset) :
    _checkLayers(pset.get<bool>("checkLayers")),
    _verboseLevel(pset.get<int>("verboseLevel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _PEthreshold(pset.get<int>("PEthreshold")),
    _maxDistance(pset.get<double>("maxDistance")),
    _maxTimeDifference(pset.get<double>("maxTimeDifference")),
    _timeWindowStart(pset.get<double>("timeWindowStart")),
    _timeWindowEnd(pset.get<double>("timeWindowEnd")),
    _coincidenceGroupsFile(pset.get<std::string>("coincidenceGroupsFile"))
  {
    produces<CrvCoincidenceCheckResult>();
    _totalTime=0;
    _totalDeadTime=0;
    _totalEvents=0;
    _totalEventsCoincidence=0;
  }

  void SimpleCrvCoincidenceCheck::beginJob()
  {
  }

  void SimpleCrvCoincidenceCheck::endJob()
  {
    std::cout<<"SUMMARY "<<_moduleLabel<<"    "<<_totalEventsCoincidence<<" / "<<_totalEvents<<" events satisfied coincidence requirements"<<std::endl;
  }

  void SimpleCrvCoincidenceCheck::beginRun(art::Run &run)
  {
    mu2e::ConditionsHandle<mu2e::AcceleratorParams> accPar("ignored");
    _microBunchPeriod = accPar->deBuncherPeriod;

    SimpleConfig config(_coincidenceGroupsFile, true, false, false);
    _nSectors = config.getInt("crs.nSectors");
    config.getVectorString("crs.sectorNames",_crvSectorNames,_nSectors);
    for(int i=0; i<_nSectors; i++) 
    {
      std::vector<int> tmp;
      config.getVectorInt("crs.sectorCoincidenceProperties"+_crvSectorNames[i], tmp, 4);
      sectorCoincidenceProperties s;
      s.coincidenceGroup[0]=tmp[0];
      s.coincidenceGroup[1]=tmp[1];
      s.moduleType=tmp[2];
      s.layerDirection=tmp[3];
      _sectorMap[i]=s;
    }
  }


  void SimpleCrvCoincidenceCheck::produce(art::Event& event) 
  {
    std::unique_ptr<CrvCoincidenceCheckResult> crvCoincidenceCheckResult(new CrvCoincidenceCheckResult);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    std::map<int, std::vector<coincidenceStruct> > coincidenceMap;

    for(CrvRecoPulsesCollection::const_iterator iter=crvRecoPulsesCollection->begin(); 
        iter!=crvRecoPulsesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);
      int   sectionNumber=CRSbar.id().getShieldNumber();
      int   layerNumber=CRSbar.id().getLayerNumber();

      const CrvRecoPulses &crvRecoPulses = iter->second;
      std::map<int,sectorCoincidenceProperties>::const_iterator sector=_sectorMap.find(sectionNumber);
      if(sector==_sectorMap.end()) continue;

      for(int SiPM=0; SiPM<2; SiPM++)
      {
        int coincidenceGroup = sector->second.coincidenceGroup[SiPM];
        if(coincidenceGroup==0) continue;

        coincidenceStruct c;
        c.counter=barIndex;
        c.layer=layerNumber;
        c.moduleType = sector->second.moduleType;
        int layerDirection = sector->second.layerDirection;
        c.pos = CRSbar.getPosition()[layerDirection];

        for(int SiPMtmp=SiPM; SiPMtmp<4; SiPMtmp+=2)  //the reco pulse times of both SiPMs on the same side are put into the same c to avoid that 
                                                      //pulses from the same counter can be used to satisfy the coincidence condition
        {
          const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &pulseVector = crvRecoPulses.GetRecoPulses(SiPMtmp);
          for(unsigned int i = 0; i<pulseVector.size(); i++) 
          {
            const CrvRecoPulses::CrvSingleRecoPulse &pulse = pulseVector[i];
            double time=pulse._leadingEdge;
            if(_verboseLevel>4)
            {
              std::cout<<"coincidence group: "<<coincidenceGroup<<"   barIndex: "<<barIndex<<"   module: "<<sectionNumber<<"  layer: "<<layerNumber<<"  SiPM: "<<SiPMtmp<<"      ";
              std::cout<<"  PEs: "<<pulse._PEs<<"   LE: "<<time<<"   pos: "<<c.pos<<std::endl;
            }
            if(pulse._PEs>=_PEthreshold && 
               time>=_timeWindowStart && 
               time<=_timeWindowEnd)
               {
                 c.time.push_back(time);
                 c.PEs.push_back(pulse._PEs);
                 c.SiPMs.push_back(SiPMtmp);
                 if(_verboseLevel==4)
                 {
                   std::cout<<"coincidence group: "<<coincidenceGroup<<"   barIndex: "<<barIndex<<"   module: "<<sectionNumber<<"  layer: "<<layerNumber<<"  SiPM: "<<SiPMtmp<<"      ";
                   std::cout<<"  PEs: "<<pulse._PEs<<"   LE: "<<time<<"   pos: "<<c.pos<<std::endl;
                 }
               }
            time+=_microBunchPeriod;
            if(pulse._PEs>=_PEthreshold && 
               time>=_timeWindowStart && 
               time<=_timeWindowEnd)
               {
                 c.time.push_back(time);
                 c.PEs.push_back(pulse._PEs);
                 c.SiPMs.push_back(SiPMtmp);
                 if(_verboseLevel==4)
                 {
                   std::cout<<"coincidence group: "<<coincidenceGroup<<"   barIndex: "<<barIndex<<"   module: "<<sectionNumber<<"  layer: "<<layerNumber<<"  SiPM: "<<SiPMtmp<<"      ";
                   std::cout<<"  PEs: "<<pulse._PEs<<"   LE: "<<time<<"   pos: "<<c.pos<<std::endl;
                 }
               }
          }
        }
        if(c.time.size()>0) coincidenceMap[coincidenceGroup].push_back(c);
      }
    }

    std::map<int, std::vector<coincidenceStruct> >::const_iterator iterC;
    for(iterC = coincidenceMap.begin(); iterC!=coincidenceMap.end(); iterC++)
    {
      const std::vector<coincidenceStruct> &vectorC = iterC->second;
      unsigned int n=vectorC.size();
      for(unsigned int i1=0; i1<n; i1++) 
      for(unsigned int i2=i1+1; i2<n; i2++) 
      for(unsigned int i3=i2+1; i3<n; i3++)
      {
        if(_checkLayers)
        {
          int layers[3]={vectorC[i1].layer,vectorC[i2].layer,vectorC[i3].layer};
          int moduleTypes[3]={vectorC[i1].moduleType,vectorC[i2].moduleType,vectorC[i3].moduleType};
          if(layers[0]==layers[1] && moduleTypes[0]==moduleTypes[1]) continue;
          if(layers[0]==layers[2] && moduleTypes[0]==moduleTypes[2]) continue;
          if(layers[1]==layers[2] && moduleTypes[1]==moduleTypes[2]) continue;
        }

        double pos[3]={vectorC[i1].pos,vectorC[i2].pos,vectorC[i3].pos};
        double posMin = *std::min_element(pos,pos+3);
        double posMax = *std::max_element(pos,pos+3);
        if(posMax-posMin>_maxDistance) continue;

        const std::vector<double> &vectorTime1 = vectorC[i1].time;
        const std::vector<double> &vectorTime2 = vectorC[i2].time;
        const std::vector<double> &vectorTime3 = vectorC[i3].time;
        for(unsigned int j1=0; j1<vectorTime1.size(); j1++)
        for(unsigned int j2=0; j2<vectorTime2.size(); j2++)
        for(unsigned int j3=0; j3<vectorTime3.size(); j3++)
        {
          double time[3]={vectorTime1[j1],vectorTime2[j2],vectorTime3[j3]};
          double timeMin = *std::min_element(time,time+3);
          double timeMax = *std::max_element(time,time+3);
          if(timeMax-timeMin<_maxTimeDifference) 
          {
              CrvCoincidenceCheckResult::CoincidenceCombination combination;
              for(int k=0; k<3; k++) combination._time[k] = time[k];
              combination._PEs[0] = vectorC[i1].PEs[j1];
              combination._PEs[1] = vectorC[i2].PEs[j2];
              combination._PEs[2] = vectorC[i3].PEs[j3];
              combination._counters[0] = vectorC[i1].counter;
              combination._counters[1] = vectorC[i2].counter;
              combination._counters[2] = vectorC[i3].counter;
              combination._SiPMs[0] = vectorC[i1].SiPMs[j1];
              combination._SiPMs[1] = vectorC[i2].SiPMs[j2];
              combination._SiPMs[2] = vectorC[i3].SiPMs[j3];
              crvCoincidenceCheckResult->GetCoincidenceCombinations().push_back(combination);
              if(_verboseLevel>2)
              {
                std::cout<<"Coincidence times/counters/SiPMs/PEs: "<<std::endl;
                for(int k=0; k<3; k++)
                {
                  std::cout<<"   "<<time[k]<<" / "<<combination._counters[k]<<" / "<<combination._SiPMs[k]<<" / "<<combination._PEs[k]<<std::endl;
                }
              }
          }
        }
      } 
    }

    _totalEvents++;
    if(crvCoincidenceCheckResult->CoincidenceFound()) _totalEventsCoincidence++;

    if(_verboseLevel>0)
    {
      _moduleLabel = *this->currentContext()->moduleLabel();
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

} // end namespace mu2e

using mu2e::SimpleCrvCoincidenceCheck;
DEFINE_ART_MODULE(SimpleCrvCoincidenceCheck)
