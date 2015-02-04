//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
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
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
#include "MCDataProducts/inc/CrvWaveformsCollection.hh"
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

#include <TRandom3.h>
#include <TMath.h>
#include <TH2D.h>

namespace mu2e 
{
  class CRVTest : public art::EDProducer 
  {
    public:
    explicit CRVTest(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string _crvPhotonArrivalsModuleLabel;
    std::string _crvSiPMResponsesModuleLabel;
    std::string _crvWaveformsModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    std::string _crvCoincidenceCheckModuleLabel;
    double      _threshold;

    TH2D  *_PEvsIntegral, *_PEvsPulseHeight;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    _crvPhotonArrivalsModuleLabel(pset.get<std::string>("crvPhotonArrivalsModuleLabel")),
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _crvCoincidenceCheckModuleLabel(pset.get<std::string>("crvCoincidenceCheckModuleLabel")),
    _threshold(pset.get<double>("threshold",0.005))
  {
    _PEvsIntegral = new TH2D("PEvsIntegral","PEvsIntegral", 60,0,60, 100,0,1.2);
    _PEvsIntegral->SetXTitle("PEs");
    _PEvsIntegral->SetYTitle("Integral");
    _PEvsPulseHeight = new TH2D("PEvsPulseHeight","PEvsPulseHeight", 60,0,60, 100,0,0.4);
    _PEvsPulseHeight->SetXTitle("PEs");
    _PEvsPulseHeight->SetYTitle("PulseHeight [mV]");
  }

  void CRVTest::beginJob()
  {
  }

  void CRVTest::endJob()
  {
//    _PEvsIntegral->SaveAs("PEvsIntegral.C");
//    _PEvsPulseHeight->SaveAs("PEvsPulseHeight.C");
  }

  void CRVTest::produce(art::Event& event) 
  {
    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection;
    event.getByLabel(_crvPhotonArrivalsModuleLabel,"",crvPhotonArrivalsCollection);

    art::Handle<CrvSiPMResponsesCollection> crvSiPMResponsesCollection;
    event.getByLabel(_crvSiPMResponsesModuleLabel,"",crvSiPMResponsesCollection);

    art::Handle<CrvWaveformsCollection> crvWaveformsCollection;
    event.getByLabel(_crvWaveformsModuleLabel,"",crvWaveformsCollection);

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    art::Handle<CrvCoincidenceCheckResult> crvCoincidenceCheckResult;
    event.getByLabel(_crvCoincidenceCheckModuleLabel,"",crvCoincidenceCheckResult);

    for(CrvPhotonArrivalsCollection::const_iterator iter1=crvPhotonArrivalsCollection->begin(); 
        iter1!=crvPhotonArrivalsCollection->end(); iter1++)
    {
      int PEs[4]={0};
      int PEsFromPulses[4]={0};
      int nPulses[4]={0};
      double singlePulseHeight[4]={0};
      double singlePulseLeadingEdge[4]={0};
      double integral[4]={0};

      const CRSScintillatorBarIndex &barIndex = iter1->first;
      const CrvPhotonArrivals &crvPhotonArrivals = iter1->second;

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        PEs[SiPM] = crvPhotonArrivals.GetNumberOfPhotonArrivals(SiPM);
      }

      CrvRecoPulsesCollection::const_iterator iter2 = crvRecoPulsesCollection->find(barIndex);
      if(iter2!=crvRecoPulsesCollection->end())
      {
        const CrvRecoPulses &crvRecoPulses = iter2->second;
        for(int SiPM=0; SiPM<4; SiPM++)
        {
          const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &pulseVector = crvRecoPulses.GetRecoPulses(SiPM);
          for(unsigned int i = 0; i<pulseVector.size(); i++) PEsFromPulses[SiPM]+=pulseVector[i]._PEs;
          nPulses[SiPM]=pulseVector.size();
          if(pulseVector.size()==1) singlePulseHeight[SiPM]=pulseVector[0]._pulseHeight;
          else singlePulseHeight[SiPM]=NAN;
          if(pulseVector.size()>=1) singlePulseLeadingEdge[SiPM]=pulseVector[0]._leadingEdge;
          else singlePulseLeadingEdge[SiPM]=NAN;
        }
      }

      CrvWaveformsCollection::const_iterator iter3 = crvWaveformsCollection->find(barIndex);
      if(iter3!=crvWaveformsCollection->end())
      {
        const CrvWaveforms &crvWaveforms = iter3->second;
        for(int SiPM=0; SiPM<4; SiPM++)
        {
          const std::vector<double> waveform = crvWaveforms.GetWaveform(SiPM);
          for(unsigned bin=0; bin<waveform.size(); bin++)
          {
            double voltage = waveform[bin];
            if(voltage>_threshold) integral[SiPM] += voltage;
          }
        }
      }

      std::cout<<"CRV bar: "<<barIndex<<std::endl;
      const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);
      int   moduleNumber=CRSbar.id().getShieldNumber();
      double pos=NAN;
      switch (moduleNumber)
      {
          case 0: 
          case 1: 
          case 2: 
          case 3: 
          case 4: 
          case 5: 
          case 6: 
          case 7: 
         case  8: 
         case  9: 
         case 10: 
         case 11: 
         case 12: pos=CRSbar.getPosition().z();
                  break;
         case 13: 
         case 14: pos=CRSbar.getPosition().y();
                  break;
         case 15: 
         case 16: 
         case 17: pos=CRSbar.getPosition().x();
                  break;
      };
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        std::cout<<"SiPM: "<<SiPM;
        std::cout<<"   PEs: "<<PEs[SiPM];
        std::cout<<"   PEsFromPulses: "<<PEsFromPulses[SiPM];
        std::cout<<"   PEs/integral: "<<PEs[SiPM]/integral[SiPM];
        std::cout<<"   PEs/pulseHeight: "<<PEs[SiPM]/singlePulseHeight[SiPM];
        std::cout<<"   nPulses: "<<nPulses[SiPM];
        std::cout<<"   leadingeEdge: "<<singlePulseLeadingEdge[SiPM];
        std::cout<<"   pos: "<<pos;
        std::cout<<std::endl;
        if(nPulses[SiPM]==1)
        {
          _PEvsIntegral->Fill(PEs[SiPM],integral[SiPM]);
          _PEvsPulseHeight->Fill(PEs[SiPM],singlePulseHeight[SiPM]);
        }
      }

    }
    std::cout<<(crvCoincidenceCheckResult->GetCoincidence()?"Coincidence satisfied":"No coincidence found")<<std::endl;
  } // end produce

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
