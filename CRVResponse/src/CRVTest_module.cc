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
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/CrvPhotonsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMChargesCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>
#include <TNtuple.h>

namespace mu2e 
{
  class CRVTest : public art::EDAnalyzer
  {
    public:
    explicit CRVTest(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);
    void beginJob();
    void endJob();

    private:
//    std::string _crvStepsModuleLabel;
    std::string _crvSiPMChargesModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    std::string _genParticleModuleLabel;

    TNtuple  *_recoPulses;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
//    _crvStepsModuleLabel(pset.get<std::string>("crvStepsModuleLabel")),
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _genParticleModuleLabel(pset.get<std::string>("genParticleModuleLabel"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CrvSingleCounter");
    _recoPulses = tfdir.make<TNtuple>("RecoPulses", "RecoPulses", "event:startX:startY:startZ:barIndex:SiPM:nRecoPulses:recoPEs:recoPulseHeight:recoPulseWidth:recoPulseTime:recoLEtime:MCPEs:chi2");
//    _recoPulses = tfdir.make<TNtuple>("RecoPulses", "RecoPulses", "event:startX:startY:startZ:barIndex:SiPM:nRecoPulses:recoPEs:recoPulseHeight:recoPulseWidth:recoPulseTime:MCPEs:ionizingEnergy:nonIonizingEnergy:energyLoss");
  }

  void CRVTest::beginJob()
  {
  }

  void CRVTest::endJob()
  {
  }

  void CRVTest::analyze(const art::Event& event) 
  {
//    art::Handle<StepPointMCCollection> crvStepsCollection;
//    event.getByLabel(_crvStepsModuleLabel,"CRV",crvStepsCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    art::Handle<GenParticleCollection> genParticleCollection;
    event.getByLabel(_genParticleModuleLabel,"",genParticleCollection);

    int eventID = event.id().event();
    CLHEP::Hep3Vector startPos = genParticleCollection->at(0).position();

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter; 
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();

      CrvRecoPulsesCollection::const_iterator    iterRecoPulses    = crvRecoPulsesCollection->find(barIndex);
      CrvSiPMChargesCollection::const_iterator   iterSiPMCharges   = crvSiPMChargesCollection->find(barIndex);

/*
      double ionizingEnergy=0;
      double nonIonizingEnergy=0;
      double energyLoss=0;
      for(size_t istep=0; istep<crvStepsCollection->size(); istep++)
      {
        StepPointMC const& step(crvStepsCollection->at(istep));
        if(step.volumeId()==barIndex.asUint())
        {
          nonIonizingEnergy+=step.nonIonizingEDep();
          ionizingEnergy+=step.ionizingEdep();
          if(step.simParticle()->id().asUint()==1) energyLoss=step.simParticle()->startMomentum().e()-step.simParticle()->endMomentum().e();
        }
      }
*/

      for(int SiPM=0; SiPM<4; SiPM++) 
      {

        int    nRecoPulses=0;
        int    recoPEs=0;
        double recoPulseHeight=0;
        double recoPulseWidth=0;
        double recoPulseTime=0;
        double recoLEtime=0;
        double MCPEs=0;
        double chi2=0;

        if(iterRecoPulses!=crvRecoPulsesCollection->end()) 
        {
          const CrvRecoPulses &crvRecoPulses = iterRecoPulses->second;

          const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &singlePulses = crvRecoPulses.GetRecoPulses(SiPM);
          nRecoPulses = singlePulses.size();
          for(size_t i=0; i<singlePulses.size(); i++)
          {
            if(recoPEs<singlePulses[i]._PEs)  //record the largest pulse to remove noise hits, after pulses, ...
            {
              recoPEs            = singlePulses[i]._PEs;
              recoPulseHeight    = singlePulses[i]._pulseHeight;
              recoPulseWidth     = singlePulses[i]._pulseWidth;
              recoPulseTime      = singlePulses[i]._pulseTime;
              recoLEtime         = singlePulses[i]._LEtime;
              chi2               = singlePulses[i]._pulseFitChi2;
            }
          }
        }

        if(iterSiPMCharges!=crvSiPMChargesCollection->end()) 
        {
          const CrvSiPMCharges &crvSiPMCharges = iterSiPMCharges->second;

          const std::vector<CrvSiPMCharges::CrvSingleCharge> &singleSiPMCharges = crvSiPMCharges.GetSiPMCharges(SiPM);
          for(size_t i=0; i<singleSiPMCharges.size(); i++) 
          {
            double charge = singleSiPMCharges[i]._chargeInPEs;
            MCPEs+=charge; 
          }
        }

        _recoPulses->Fill(eventID,startPos.x(),startPos.y(),startPos.z(),barIndex.asInt(),SiPM,nRecoPulses,recoPEs,recoPulseHeight,recoPulseWidth,recoPulseTime,recoLEtime,MCPEs,chi2);
//        _recoPulses->Fill(eventID,startPos.x(),startPos.y(),startPos.z(),barIndex.asInt(),SiPM,nRecoPulses,recoPEs,recoPulseHeight,recoPulseWidth,recoPulseTime,MCPEs,ionizingEnergy,nonIonizingEnergy,energyLoss);
      }
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
