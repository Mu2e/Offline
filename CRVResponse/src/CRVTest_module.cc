//
// A module to extract number of PEs, arrival times, hit positions, etc. from the CRV waveforms
//
// 
// Original Author: Ralf Ehrlich

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/CrvPhotons.hh"
#include "MCDataProducts/inc/CrvSiPMCharges.hh"
#include "RecoDataProducts/inc/CrvRecoPulse.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
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
    std::string _crvStepsModuleLabel;
    std::string _crvSiPMChargesModuleLabel;
    std::string _crvRecoPulsesModuleLabel;

    TNtuple  *_recoPulses;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _crvStepsModuleLabel(pset.get<std::string>("crvStepsModuleLabel")),
    _crvSiPMChargesModuleLabel(pset.get<std::string>("crvSiPMChargesModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CrvSingleCounter");
    _recoPulses = tfdir.make<TNtuple>("RecoPulses", "RecoPulses", "event:barIndex:SiPM:nRecoPulses:recoPEs:recoPulseHeight:recoPulseWidth:recoPulseTime:recoLEtime:chi2:MCPEs:visibleEnergyDeposited");
  }

  void CRVTest::beginJob()
  {
  }

  void CRVTest::endJob()
  {
  }

  void CRVTest::analyze(const art::Event& event) 
  {
    art::Handle<CrvStepCollection> crvStepsCollection;
    event.getByLabel(_crvStepsModuleLabel,"",crvStepsCollection);

    art::Handle<CrvSiPMChargesCollection> crvSiPMChargesCollection;
    event.getByLabel(_crvSiPMChargesModuleLabel,"",crvSiPMChargesCollection);

    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection);

    int eventID = event.id().event();

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator iter; 
    for(iter=counters.begin(); iter!=counters.end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = (*iter)->index();

      double visibleEnergyDeposited=0;
      if(crvStepsCollection.isValid())
      {
        for(size_t istep=0; istep<crvStepsCollection->size(); istep++)
        {
          CrvStep const& step(crvStepsCollection->at(istep));
          if(step.barIndex()==barIndex) visibleEnergyDeposited+=step.visibleEDep();
        }
      }

      for(int SiPM=0; SiPM<4; SiPM++) 
      {

        int    nRecoPulses=0;
        int    recoPEs=0;
        double recoPulseHeight=0;
        double recoPulseBeta=0;
        double recoPulseTime=0;
        double recoLEtime=0;
        double MCPEs=0;
        double chi2=0;

        for(size_t recoPulseIndex=0; recoPulseIndex<crvRecoPulseCollection->size(); recoPulseIndex++)
        {
          const CrvRecoPulse &crvRecoPulse = crvRecoPulseCollection->at(recoPulseIndex);
          if(crvRecoPulse.GetScintillatorBarIndex()==barIndex && crvRecoPulse.GetSiPMNumber()==SiPM) 
          {
            nRecoPulses++;
            if(recoPEs<crvRecoPulse.GetPEs())  //record the largest pulse to remove noise hits, after pulses, ...
            {
              recoPEs            = crvRecoPulse.GetPEs();
              recoPulseHeight    = crvRecoPulse.GetPulseHeight();
              recoPulseBeta      = crvRecoPulse.GetPulseBeta();
              recoPulseTime      = crvRecoPulse.GetPulseTime();
              recoLEtime         = crvRecoPulse.GetLEtime();
              chi2               = crvRecoPulse.GetPulseFitChi2();
            }
          }
        }

        for(size_t chargeIndex=0; chargeIndex<crvSiPMChargesCollection->size(); chargeIndex++)
        {
          const CrvSiPMCharges &crvSiPMCharges = crvSiPMChargesCollection->at(chargeIndex);
          if(crvSiPMCharges.GetScintillatorBarIndex()==barIndex && crvSiPMCharges.GetSiPMNumber()==SiPM) 
          {
            const std::vector<CrvSiPMCharges::SingleCharge> &singleCharges = crvSiPMCharges.GetCharges();
            for(size_t i=0; i<singleCharges.size(); i++) 
            {
              double charge = singleCharges[i]._chargeInPEs;
              MCPEs+=charge; 
            }
          }
        }

        _recoPulses->Fill(eventID,barIndex.asInt(),SiPM,nRecoPulses,recoPEs,recoPulseHeight,recoPulseBeta,recoPulseTime,recoLEtime,chi2,MCPEs,visibleEnergyDeposited);
      }
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
