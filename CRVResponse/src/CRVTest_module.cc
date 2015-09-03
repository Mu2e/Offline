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
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
#include "MCDataProducts/inc/CrvWaveformsCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"

#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
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
//    std::string _crvPhotonArrivalsModuleLabel;
//    std::string _crvSiPMResponsesModuleLabel;
//    std::string _crvWaveformsModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    std::string _genParticleModuleLabel;

    TNtuple  *_recoPulses;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
//    _crvPhotonArrivalsModuleLabel(pset.get<std::string>("crvPhotonArrivalsModuleLabel")),
//    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
//    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _genParticleModuleLabel(pset.get<std::string>("genParticleModuleLabel"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CrvSingleCounter");
    _recoPulses = tfdir.make<TNtuple>( "RecoPulses",    "RecoPulses",  "SiPM:startX:startZ:recoPEs" );
  }

  void CRVTest::beginJob()
  {
  }

  void CRVTest::endJob()
  {
  }

  void CRVTest::analyze(const art::Event& event) 
  {
//    art::Handle<CrvPhotonArrivalsCollection> crvPhotonArrivalsCollection;
//    event.getByLabel(_crvPhotonArrivalsModuleLabel,"",crvPhotonArrivalsCollection);

//    art::Handle<CrvSiPMResponsesCollection> crvSiPMResponsesCollection;
//    event.getByLabel(_crvSiPMResponsesModuleLabel,"",crvSiPMResponsesCollection);

//    art::Handle<CrvWaveformsCollection> crvWaveformsCollection;
//    event.getByLabel(_crvWaveformsModuleLabel,"",crvWaveformsCollection);

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    art::Handle<GenParticleCollection> genParticleCollection;
    event.getByLabel(_genParticleModuleLabel,"",genParticleCollection);

    if(crvRecoPulsesCollection->size()>1) throw cet::exception("analyze")<< "Should never occur: more than 1 crv counters\n";
    if(genParticleCollection->size()!=1) throw cet::exception("analyze")<< "Should never occur: more (or less) than 1 gen particles\n";

    double startX = genParticleCollection->at(0).position().x();
    double startZ = genParticleCollection->at(0).position().z();
    int recoPEs[4]={0};

    for(CrvRecoPulsesCollection::const_iterator iter=crvRecoPulsesCollection->begin(); //this is intended for only one counter
        iter!=crvRecoPulsesCollection->end(); iter++)                                  //so there should be only one entry
    {
      const CrvRecoPulses &crvRecoPulses = iter->second;

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &singlePulses = crvRecoPulses.GetRecoPulses(SiPM);
        for(size_t i=0; i<singlePulses.size(); i++) recoPEs[SiPM] += singlePulses[i]._PEs;
      }
    }

    for(int SiPM=0; SiPM<4; SiPM++)
    {
      _recoPulses->Fill(SiPM,startX,startZ,recoPEs[SiPM]);
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
