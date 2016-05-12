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
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
#include "MCDataProducts/inc/CrvWaveformsCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"

#include "art/Persistency/Common/Ptr.h"
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
    std::string _crvSiPMResponsesModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    std::string _genParticleModuleLabel;
    std::string _simParticleModuleLabel;

    TNtuple  *_recoPulses, *_leadingEdgesGlobal, *_leadingEdgesCounter;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _genParticleModuleLabel(pset.get<std::string>("genParticleModuleLabel"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CrvSingleCounter");
    _recoPulses = tfdir.make<TNtuple>( "RecoPulses",    "RecoPulses",  "event:startX:startY:startZ:barIndex:SiPM:nRecoPulses:recoPEs:recoPulseHeight:recoPulseIntegral:MCPEs" );
  }

  void CRVTest::beginJob()
  {
  }

  void CRVTest::endJob()
  {
  }

  void CRVTest::analyze(const art::Event& event) 
  {
    art::Handle<CrvSiPMResponsesCollection> crvSiPMResponsesCollection;
    event.getByLabel(_crvSiPMResponsesModuleLabel,"",crvSiPMResponsesCollection);

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
      CrvSiPMResponsesCollection::const_iterator iterSiPMResponses = crvSiPMResponsesCollection->find(barIndex);

      for(int SiPM=0; SiPM<4; SiPM++) 
      {

        int    nRecoPulses=0;
        int    recoPEs=0;
        double recoPulseHeights=0;
        double recoPulseIntegrals=0;
        double MCPEs=0;

        if(iterRecoPulses!=crvRecoPulsesCollection->end()) 
        {
          const CrvRecoPulses &crvRecoPulses = iterRecoPulses->second;

          const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &singlePulses = crvRecoPulses.GetRecoPulses(SiPM);
          nRecoPulses = singlePulses.size();
          if(singlePulses.size()>0)
          {
            recoPEs            = singlePulses[0]._PEs;
            recoPulseHeights   = singlePulses[0]._pulseHeight;
            recoPulseIntegrals = singlePulses[0]._integral;
          }
        }

        if(iterSiPMResponses!=crvSiPMResponsesCollection->end()) 
        {
          const CrvSiPMResponses &crvSiPMResponses = iterSiPMResponses->second;

          const std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &singleSiPMResponses = crvSiPMResponses.GetSiPMResponses(SiPM);
          for(size_t i=0; i<singleSiPMResponses.size(); i++) MCPEs += singleSiPMResponses[i]._charge;
        }

        _recoPulses->Fill(eventID,startPos.x(),startPos.y(),startPos.z(),barIndex.asInt(),SiPM,nRecoPulses,recoPEs,recoPulseHeights,recoPulseIntegrals,MCPEs);
      }
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
