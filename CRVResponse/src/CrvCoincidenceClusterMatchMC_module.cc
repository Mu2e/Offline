//
// A module to find the MC information for the clusters of coincidences of CRV pulses
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "Offline/CRVResponse/inc/CrvMCHelper.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TNtuple.h>

namespace mu2e
{
  class CrvCoincidenceClusterMatchMC : public art::EDProducer
  {
    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvCoincidenceClusterFinderModuleLabel{Name("crvCoincidenceClusterFinderModuleLabel"),
                                                                      Comment("label of CoincidenceClusterFinder module")};
      fhicl::Atom<std::string> crvWaveformsModuleLabel{Name("crvWaveformsModuleLabel"), Comment("label of CrvDigiMC module")};
      fhicl::Atom<std::string> ewmLabel{Name("eventWindowMarkerModuleLabel"), Comment("label of EventWindowMarker module")};
      fhicl::Atom<bool> doNtuples{Name("doNtuples"), Comment("produce TNtuples")};
    };
    typedef art::EDProducer::Table<Config> Parameters;

    explicit CrvCoincidenceClusterMatchMC(const Parameters& config);
    void produce(art::Event& e);
    void beginJob();
    void beginRun(art::Run &run);
    void endJob();

    private:
    std::string _crvCoincidenceClusterFinderModuleLabel;  //module label of the CrvCoincidenceClusterFinder module
                                                          //it is possible to have more than one instance of the CrvCoincidenceClusterFinder module
    std::string _crvWaveformsModuleLabel;  //module label of the CrvWaveform module.
                                           //this is optional. only needed, if MC information is required
    std::string _ewmLabel; //module label of the EventWindowMarker module

    bool        _doNtuples;
    TNtuple    *_ntuple;
    TNtuple    *_ntupleHitTimes;
  };

  CrvCoincidenceClusterMatchMC::CrvCoincidenceClusterMatchMC(const Parameters& config) :
    art::EDProducer(config),
    _crvCoincidenceClusterFinderModuleLabel(config().crvCoincidenceClusterFinderModuleLabel()),
    _crvWaveformsModuleLabel(config().crvWaveformsModuleLabel()),
    _ewmLabel(config().ewmLabel()),
    _doNtuples(config().doNtuples())
  {
    produces<CrvCoincidenceClusterMCCollection>();

    if(_doNtuples)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("CrvClusterMCmatch");
      _ntuple = tfdir.make<TNtuple>("CrvClusters", "CrvCluster", "run:subrun:event:valid:crvSector:x:y:z:xMC:yMC:zMC:timeDiff:posDiff:eDep:PEs");
      _ntupleHitTimes = tfdir.make<TNtuple>("HitTimes", "HitTime", "crvSector:sipm:hitZ:hitT:recoTime:LETime:recoPE:noRecoFlags");
    }
  }

  void CrvCoincidenceClusterMatchMC::beginJob()
  {
  }

  void CrvCoincidenceClusterMatchMC::endJob()
  {
  }

  void CrvCoincidenceClusterMatchMC::beginRun(art::Run &run)
  {
  }

  void CrvCoincidenceClusterMatchMC::produce(art::Event& event)
  {

    std::unique_ptr<CrvCoincidenceClusterMCCollection> crvCoincidenceClusterMCCollection(new CrvCoincidenceClusterMCCollection);

    art::Handle<CrvCoincidenceClusterCollection> crvCoincidenceClusterCollection;
    event.getByLabel(_crvCoincidenceClusterFinderModuleLabel,"",crvCoincidenceClusterCollection);

    if(crvCoincidenceClusterCollection.product()==NULL) return;

    art::Handle<CrvDigiMCCollection> crvDigiMCCollection;
    if(_crvWaveformsModuleLabel!="") event.getByLabel(_crvWaveformsModuleLabel,"",crvDigiMCCollection); //this is an optional part for MC information

    CrvCoincidenceClusterCollection::const_iterator iter;
    for(iter=crvCoincidenceClusterCollection->begin(); iter!=crvCoincidenceClusterCollection->end(); iter++)
    {
      bool   hasMCInfo                = (crvDigiMCCollection.isValid()?true:false); //MC
      double visibleEnergyDeposited   = 0;         //MC
      double earliestHitTime          = 0;         //MC
      double avgHitTime               = 0;         //MC
      art::Ptr<SimParticle> simParticle;           //MC
      CLHEP::Hep3Vector     earliestHitPos;        //MC
      CLHEP::Hep3Vector     avgHitPos;             //MC

      //get the event window marker
      art::Handle<mu2e::EventWindowMarker> ewmh;
      event.getByLabel(_ewmLabel, ewmh);

      //loop through all reco pulses and try to find the MC information
      std::vector<CrvCoincidenceClusterMC::PulseInfo> pulses; //collection of all pulses (sim particle, energy dep.)
      std::set<art::Ptr<CrvStep> > steps;  //collection of all crvSteps for total energy calculation
                                           //use a set to avoid double counts, e.g. for two pulses coming from same digi
      const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses = iter->GetCrvRecoPulses();
      for(size_t i=0; i<crvRecoPulses.size(); i++)
      {
        const art::Ptr<CrvRecoPulse> crvRecoPulse = crvRecoPulses[i];
        art::Ptr<SimParticle> simParticleThisPulse;
        double visibleEnergyDepositedThisPulse = 0;
        double earliestHitTimeThisPulse = 0; //not used here
        double avgHitTimeThisPulse = 0; //not used here
        CLHEP::Hep3Vector earliestHitPosThisPulse; //not used here
        CLHEP::Hep3Vector avgHitPosThisPulse; //not used here

        //get MC information, if available
        if(hasMCInfo)
        {
          //get the step points of this reco pulse and add it to the collection of all step points
          CrvMCHelper::GetStepPointsFromCrvRecoPulse(crvRecoPulse, crvDigiMCCollection, steps);

          //get the sim particle and deposited energy of this reco pulse
          CrvMCHelper::GetInfoFromCrvRecoPulse(crvRecoPulse, crvDigiMCCollection, visibleEnergyDepositedThisPulse,
                                               earliestHitTimeThisPulse, earliestHitPosThisPulse, avgHitTimeThisPulse,
                                               avgHitPosThisPulse, simParticleThisPulse, ewmh);
          if(_doNtuples) _ntupleHitTimes->Fill(iter->GetCrvSectorType(),crvRecoPulse->GetSiPMNumber(),avgHitPosThisPulse.x(),avgHitTimeThisPulse,
                                               crvRecoPulse->GetPulseTime(),crvRecoPulse->GetLEtime(),crvRecoPulse->GetPEs(),
                                               crvRecoPulse->GetRecoPulseFlags().none());
        }

        //add the MC information (sim particle, dep. energy) of this reco pulse to the collection of pulses
        //unless the reco pulse was caused by a noise hit (and doesn't have an associated simParticle)
        pulses.emplace_back(simParticleThisPulse,visibleEnergyDepositedThisPulse);
      }//loop over reco pulses

      //based on all step points, get the most likely sim particle, total energy, etc.
      CrvMCHelper::GetInfoFromStepPoints(steps, visibleEnergyDeposited, earliestHitTime, earliestHitPos,
                                         avgHitTime, avgHitPos, simParticle, ewmh);

      //insert the cluster information into the vector of the crv coincidence clusters (collection of pulses, most likely sim particle, etc.)
      //if the cluster was caused by noise, the simParticle will be null, the visibleEnergyDeposited, earliest/average hit time and hit position will be 0
      crvCoincidenceClusterMCCollection->emplace_back(hasMCInfo, pulses, simParticle, visibleEnergyDeposited,
                                                      earliestHitTime, earliestHitPos, avgHitTime, avgHitPos);
      if(_doNtuples)
      {
          _ntuple->Fill(event.run(),event.subRun(),event.event(),
                        iter->HasTwoReadoutSides(),iter->GetCrvSectorType(),iter->GetAvgHitPos().x(),iter->GetAvgHitPos().y(),iter->GetAvgHitPos().z(),
                        avgHitPos.x(),avgHitPos.y(),avgHitPos.z(),
                        iter->GetAvgHitTime()-avgHitTime,(iter->GetAvgHitPos()-avgHitPos).mag(),
                        visibleEnergyDeposited,iter->GetSidePEs()[0]+iter->GetSidePEs()[1]);
      }
    }//loop over all clusters

    event.put(std::move(crvCoincidenceClusterMCCollection));
  } // end produce

} // end namespace mu2e

using mu2e::CrvCoincidenceClusterMatchMC;
DEFINE_ART_MODULE(CrvCoincidenceClusterMatchMC)
