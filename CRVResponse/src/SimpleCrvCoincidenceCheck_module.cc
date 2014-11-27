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

#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "RecoDataProducts/inc/CRVRecoPulsesCollection.hh"
#include "RecoDataProducts/inc/CRVCoincidenceCheckResult.hh"

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

#include <TRandom3.h>
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
    void endJob();

    private:
    std::string _crvRecoPulsesModuleLabel;
    double      _PEthreshold;
    double      _maxDistance;
    double      _maxTimeDifference;

    struct coincidenceStruct
    {
      double pos;
      std::vector<double> time;
    };
  };

  SimpleCrvCoincidenceCheck::SimpleCrvCoincidenceCheck(fhicl::ParameterSet const& pset) :
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _PEthreshold(pset.get<double>("PEthreshold")),
    _maxDistance(pset.get<double>("maxDistance")),
    _maxTimeDifference(pset.get<double>("maxTimeDifference"))
  {
    produces<CRVCoincidenceCheckResult>();
  }

  void SimpleCrvCoincidenceCheck::beginJob()
  {
  }

  void SimpleCrvCoincidenceCheck::endJob()
  {
  }

  void SimpleCrvCoincidenceCheck::produce(art::Event& event) 
  {
    std::unique_ptr<CRVCoincidenceCheckResult> crvCoincidenceCheckResult(new CRVCoincidenceCheckResult);

    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CRVRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    std::map<int, std::vector<coincidenceStruct> > coincidenceMap;

    for(CRVRecoPulsesCollection::const_iterator iter=crvRecoPulsesCollection->begin(); 
        iter!=crvRecoPulsesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);
      int   moduleNumber=CRSbar.id().getShieldNumber();

      const CRVRecoPulses &crvRecoPulses = iter->second;
      for(int SiPM=0; SiPM<2; SiPM++)
      {
        coincidenceStruct c;
        int coincidenceGroup = 0;
        switch (moduleNumber)
        {
          case 0: 
          case 3: 
          case 4: 
          case 5: if(SiPM==0) {coincidenceGroup = 1; c.pos=CRSbar.getPosition().z();}
                  if(SiPM==1) {coincidenceGroup = 2; c.pos=CRSbar.getPosition().z();}
                  break;
          case 1: if(SiPM==1) {coincidenceGroup = 2; c.pos=CRSbar.getPosition().z();}
                  break;
          case 2: if(SiPM==0) {coincidenceGroup = 1; c.pos=CRSbar.getPosition().z();}
                  if(SiPM==1) {coincidenceGroup = 3; c.pos=CRSbar.getPosition().z();}
                  break;
          case 6: 
          case 7: 
          case 8: if(SiPM==0) {coincidenceGroup = 4; c.pos=CRSbar.getPosition().z();}
                  if(SiPM==1) {coincidenceGroup = 5; c.pos=CRSbar.getPosition().z();}
                  break;
          case 9: if(SiPM==0) {coincidenceGroup = 2; c.pos=CRSbar.getPosition().z();}
                  break;
         case 10: 
         case 11: 
         case 12: if(SiPM==0) {coincidenceGroup = 2; c.pos=CRSbar.getPosition().z();}
                  if(SiPM==1) {coincidenceGroup = 5; c.pos=CRSbar.getPosition().z();}
                  break;
         case 13: if(SiPM==0) {coincidenceGroup = 6; c.pos=CRSbar.getPosition().y();}
                  if(SiPM==1) {coincidenceGroup = 7; c.pos=CRSbar.getPosition().y();}
                  break;
         case 14: if(SiPM==0) {coincidenceGroup = 8; c.pos=CRSbar.getPosition().y();}
                  break;
         case 15: if(SiPM==0) {coincidenceGroup = 9; c.pos=CRSbar.getPosition().x();}
                  break;
         case 16: if(SiPM==0) {coincidenceGroup = 10; c.pos=CRSbar.getPosition().x();}
                  break;
         case 17: if(SiPM==0) {coincidenceGroup = 11; c.pos=CRSbar.getPosition().x();}
                  break;
        };
//std::cout<<"coincidence group: "<<coincidenceGroup<<std::endl;
        const std::vector<CRVRecoPulses::CRVSingleRecoPulse> &pulseVector1 = crvRecoPulses.GetRecoPulses(SiPM);
        for(unsigned int i = 0; i<pulseVector1.size(); i++) 
        {
          if(pulseVector1[i]._PEs>=_PEthreshold) c.time.push_back(pulseVector1[i]._leadingEdge);
//std::cout<<"PEs: "<<pulseVector1[i]._PEs<<"   LE: "<<pulseVector1[i]._leadingEdge<<"   pos: "<<c.pos<<std::endl;
        }
        const std::vector<CRVRecoPulses::CRVSingleRecoPulse> &pulseVector2 = crvRecoPulses.GetRecoPulses(SiPM+2);
        for(unsigned int i = 0; i<pulseVector2.size(); i++) 
        {
          if(pulseVector2[i]._PEs>=_PEthreshold) c.time.push_back(pulseVector2[i]._leadingEdge);
//std::cout<<"PEs: "<<pulseVector2[i]._PEs<<"   LE: "<<pulseVector2[i]._leadingEdge<<"   pos: "<<c.pos<<std::endl;
        }
        coincidenceMap[coincidenceGroup].push_back(c);
      }
    }

    bool foundCoincidence=false;
    std::map<int, std::vector<coincidenceStruct> >::const_iterator iterC;
    for(iterC = coincidenceMap.begin(); iterC!=coincidenceMap.end() && foundCoincidence==false; iterC++)
    {
      const std::vector<coincidenceStruct> &vectorC = iterC->second;
      unsigned int n=vectorC.size();
      for(unsigned int i1=0; i1<n && foundCoincidence==false; i1++) 
      for(unsigned int i2=i1+1; i2<n && foundCoincidence==false; i2++) 
      for(unsigned int i3=i2+1; i3<n && foundCoincidence==false; i3++)
      {
        double pos[3]={vectorC[i1].pos,vectorC[i2].pos,vectorC[i3].pos};
        double posMin = *std::min_element(pos,pos+3);
        double posMax = *std::max_element(pos,pos+3);
        if(posMax-posMin>_maxDistance) continue;

        const std::vector<double> &vectorTime1 = vectorC[i1].time;
        const std::vector<double> &vectorTime2 = vectorC[i2].time;
        const std::vector<double> &vectorTime3 = vectorC[i3].time;
        for(unsigned int j1=0; j1<vectorTime1.size() && foundCoincidence==false; j1++)
        for(unsigned int j2=0; j2<vectorTime2.size() && foundCoincidence==false; j2++)
        for(unsigned int j3=0; j3<vectorTime3.size() && foundCoincidence==false; j3++)
        {
          double time[3]={vectorTime1[j1],vectorTime2[j2],vectorTime3[j3]};
          double timeMin = *std::min_element(time,time+3);
          double timeMax = *std::max_element(time,time+3);
          if(timeMax-timeMin<_maxTimeDifference) foundCoincidence=true;
        }
      } 
    }
    crvCoincidenceCheckResult->SetCoincidence(foundCoincidence);
std::cout<<"run "<<event.id().run()<<"  subrun "<<event.id().subRun()<<"  event "<<event.id().event()<<"    ";
std::cout<<(foundCoincidence?"Coincidence satisfied":"No coincidence found")<<std::endl;
    event.put(std::move(crvCoincidenceCheckResult));
  } // end produce

} // end namespace mu2e

using mu2e::SimpleCrvCoincidenceCheck;
DEFINE_ART_MODULE(SimpleCrvCoincidenceCheck)
