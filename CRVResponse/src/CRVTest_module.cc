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
//    std::string _crvPhotonArrivalsModuleLabel;
    std::string _crvSiPMResponsesModuleLabel;
//    std::string _crvWaveformsModuleLabel;
    std::string _crvRecoPulsesModuleLabel;
    std::string _genParticleModuleLabel;
    std::string _simParticleModuleLabel;

    TNtuple  *_recoPulses, *_leadingEdgesGlobal, *_leadingEdgesCounter;
  };

  CRVTest::CRVTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
//    _crvPhotonArrivalsModuleLabel(pset.get<std::string>("crvPhotonArrivalsModuleLabel")),
    _crvSiPMResponsesModuleLabel(pset.get<std::string>("crvSiPMResponsesModuleLabel")),
//    _crvWaveformsModuleLabel(pset.get<std::string>("crvWaveformsModuleLabel")),
    _crvRecoPulsesModuleLabel(pset.get<std::string>("crvRecoPulsesModuleLabel")),
    _genParticleModuleLabel(pset.get<std::string>("genParticleModuleLabel")),
    _simParticleModuleLabel(pset.get<std::string>("simParticleModuleLabel"))
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CrvSingleCounter");
    _recoPulses = tfdir.make<TNtuple>( "RecoPulses",    "RecoPulses",  "SiPM:startX:startZ:recoPEs:PEs:theta:phi" );
    _leadingEdgesGlobal = tfdir.make<TNtuple>( "LeadingEdgesGlobal",    "LeadingEdgesGlobal",  "timeDifferenceGlobal:trackLength" );
    _leadingEdgesCounter = tfdir.make<TNtuple>( "LeadingEdgesCounter",    "LeadingEdgesCounter",  "timeDifferenceCounter" );
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

    art::Handle<CrvSiPMResponsesCollection> crvSiPMResponsesCollection;
    event.getByLabel(_crvSiPMResponsesModuleLabel,"",crvSiPMResponsesCollection);

//    art::Handle<CrvWaveformsCollection> crvWaveformsCollection;
//    event.getByLabel(_crvWaveformsModuleLabel,"",crvWaveformsCollection);

    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulsesCollection);

    art::Handle<GenParticleCollection> genParticleCollection;
    event.getByLabel(_genParticleModuleLabel,"",genParticleCollection);

    art::Handle<SimParticleCollection> simParticleCollection;
    event.getByLabel(_simParticleModuleLabel,"",simParticleCollection);

    double startX = genParticleCollection->at(0).position().x();
    double startZ = genParticleCollection->at(0).position().z();

    int recoPEs[4]={0};
    double PEs[4]={0};

    CrvRecoPulsesCollection::const_iterator iterRecoPulses=crvRecoPulsesCollection->begin(); //this is intended for only one counter
    if(iterRecoPulses!=crvRecoPulsesCollection->end()) 
    {
      const CrvRecoPulses &crvRecoPulses = iterRecoPulses->second;

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &singlePulses = crvRecoPulses.GetRecoPulses(SiPM);
        for(size_t i=0; i<singlePulses.size(); i++) recoPEs[SiPM] += singlePulses[i]._PEs;
      }
    }

    CrvSiPMResponsesCollection::const_iterator iterSiPMResponses=crvSiPMResponsesCollection->begin(); //this is intended for only one counter
    if(iterSiPMResponses!=crvSiPMResponsesCollection->end()) 
    {
      const CrvSiPMResponses &crvSiPMResponses = iterSiPMResponses->second;

      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<CrvSiPMResponses::CrvSingleSiPMResponse> &singleSiPMResponses = crvSiPMResponses.GetSiPMResponses(SiPM);
        for(size_t i=0; i<singleSiPMResponses.size(); i++) PEs[SiPM] += singleSiPMResponses[i]._charge;
      }
    }

    double theta=acos(simParticleCollection->getOrThrow(cet::map_vector_key(1)).startMomentum().vect().unit().y());
    double phi=atan2(simParticleCollection->getOrThrow(cet::map_vector_key(1)).startMomentum().vect().z(),simParticleCollection->getOrThrow(cet::map_vector_key(1)).startMomentum().vect().x());

    for(int SiPM=0; SiPM<4; SiPM++)
    {
      _recoPulses->Fill(SiPM,startX,startZ,recoPEs[SiPM],PEs[SiPM],theta,phi);
    }


    double firstTimeGlobal1=NAN;
    double firstTimeGlobal2=NAN;
    double lastTimeGlobal1=NAN;
    double lastTimeGlobal2=NAN;
    double trackLength = (simParticleCollection->getOrThrow(cet::map_vector_key(1)).startPosition() - simParticleCollection->getOrThrow(cet::map_vector_key(1)).endPosition()).mag();
    iterRecoPulses=crvRecoPulsesCollection->begin();
    for(; iterRecoPulses!=crvRecoPulsesCollection->end(); iterRecoPulses++) 
    {
      const CrvRecoPulses &crvRecoPulses = iterRecoPulses->second;
      double firstTimeCounter[4]={NAN,NAN,NAN,NAN};
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &singlePulses = crvRecoPulses.GetRecoPulses(SiPM);
        for(size_t i=0; i<singlePulses.size(); i++)
        {
          double time = singlePulses[i]._leadingEdge;
          if(isnan(firstTimeCounter[SiPM]) || firstTimeCounter[SiPM]>time) firstTimeCounter[SiPM]=time;
        }
      }
      double timeDifferenceCounter = fabs(firstTimeCounter[0]-firstTimeCounter[2]);
      _leadingEdgesCounter->Fill(timeDifferenceCounter);
      timeDifferenceCounter = fabs(firstTimeCounter[1]-firstTimeCounter[3]);
      _leadingEdgesCounter->Fill(timeDifferenceCounter);

      double firstTime=std::min(firstTimeCounter[0],firstTimeCounter[2]);     
      double lastTime=std::max(firstTimeCounter[0],firstTimeCounter[2]);     
      if(isnan(firstTimeGlobal1) || firstTimeGlobal1>firstTime) firstTimeGlobal1=firstTime;
      if(isnan(lastTimeGlobal1) || lastTimeGlobal1<lastTime) lastTimeGlobal1=lastTime;

      firstTime=std::min(firstTimeCounter[1],firstTimeCounter[3]);     
      lastTime=std::max(firstTimeCounter[1],firstTimeCounter[3]);     
      if(isnan(firstTimeGlobal2) || firstTimeGlobal2>firstTime) firstTimeGlobal2=firstTime;
      if(isnan(lastTimeGlobal2) || lastTimeGlobal2<lastTime) lastTimeGlobal2=lastTime;
    }
    double timeDifferenceGlobal = lastTimeGlobal1-firstTimeGlobal1;
    _leadingEdgesGlobal->Fill(timeDifferenceGlobal,trackLength);
    timeDifferenceGlobal = lastTimeGlobal2-firstTimeGlobal2;
    _leadingEdgesGlobal->Fill(timeDifferenceGlobal,trackLength);


  } // end analyze

} // end namespace mu2e

using mu2e::CRVTest;
DEFINE_ART_MODULE(CRVTest)
