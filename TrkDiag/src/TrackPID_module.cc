//
// Reconstruction-level PID determination from a track.  The current implementation
// just uses the calorimeter matching information (TrkCaloHit), eventually it will
// also process dE/dx, etc.
//
// Original author: Dave Brown (LBNL)
//

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
// utilities
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
// data
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"
#include "Offline/TrkDiag/inc/TrackPID.hxx"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std;

namespace TMVA_SOFIE_TrackPID {
  class Session;
}

namespace mu2e {
  class TrackPID : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int> debug{ Name("debugLevel"),
          Comment("Debug Level"), 0};
        fhicl::Atom<float> MaxDE{ Name("MaxDE"),
          Comment("Maximum difference between calorimeter cluster EDep energy and the track energy (assuming electron mass)")};
        fhicl::Atom<float> DT{ Name("DeltaTOffset"),
          Comment("Track - Calorimeter time offset")}; // this should be a condition FIXME!
        fhicl::Atom<art::InputTag> KSC { Name("KalSeedCollection"),
          Comment("KalSeedCollection producer")};
        fhicl::Table<typename MVATools::Config> MVAConfig{fhicl::Name("MVAConfig"),
          fhicl::Comment("MVA Configuration") };
        fhicl::Atom<std::string> datFilename{Name("datFilename"), Comment("Fiename for the .dat file to use")};
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit TrackPID(const Parameters& conf);

    private:
      void produce(art::Event& event) override;
      bool _debug;
      float _maxde, _dtoffset;
      art::InputTag _kalSeedTag;
      MVATools* _tchmva;

      std::shared_ptr<TMVA_SOFIE_TrackPID::Session> mva_;
  };

  TrackPID::TrackPID(const Parameters& config ) :
    art::EDProducer(config),
    _debug(config().debug()),
    _maxde(config().MaxDE()),
    _dtoffset(config().DT()),
    _kalSeedTag(config().KSC()),
    _tchmva(new MVATools(config().MVAConfig()))
  {
    _tchmva->initMVA();
    if(_debug> 0)_tchmva->showMVA();
    produces<MVAResultCollection>();

    ConfigFileLookupPolicy configFile;
    mva_ = std::make_shared<TMVA_SOFIE_TrackPID::Session>(ConfigFile(config().datFilename()));
  }

  void TrackPID::produce(art::Event& event ) {
    mu2e::GeomHandle<mu2e::Calorimeter> calo;
    // create output
    unique_ptr<MVAResultCollection> mvacol(new MVAResultCollection());
    // get the KalSeeds
    art::Handle<KalSeedCollection> kalSeedHandle;
    event.getByLabel(_kalSeedTag, kalSeedHandle);
    const auto& kalSeeds = *kalSeedHandle;

    for (const auto& kseed : kalSeeds) {
      TrkCaloHitPID tchpid;
      tchpid.setMVAStatus(MVAStatus::unset);
      tchpid.setMVAValue(-1.0);
      std::array<float,4> features; // features used for training

      static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
      if (kseed.status().hasAllProperties(goodfit)){
        if(kseed.hasCaloCluster() &&
            kseed.caloHit()._flag.hasAllProperties(StrawHitFlag::active)){
          auto const& tchs = kseed.caloHit();
          auto const& cc = tchs.caloCluster();
          XYZVectorF trkmom;
          auto ikseg = kseed.nearestSegment(tchs._rptoca);
          if(ikseg != kseed.segments().end())
            ikseg->mom(ikseg->localFlt(tchs.trkLen()),trkmom);
          else
            throw cet::exception("RECO")<<"mu2e::TrackPID: KalSeed segment missing" << endl;
          // compute the energy difference, assuming an electron mass
          features[0] = cc.energyDep() = sqrt(trkmom.Mag2())
          //features[1] = tchs._trkdepth()
          // move into detector coordinates.  Yikes!!
          XYZVectorF cpos = XYZVectorF(calo->geomUtil().mu2eToTracker(calo->geomUtil().diskFFToMu2e( cc->diskID(), cc->cog3Vector())));
          features[1] = sqrt(cpos.Perp2());
          // compute transverse direction WRT position
          cpos.SetZ(0.0);
          trkmom.SetZ(0.0);
          features[2] = cpos.Dot(trkmom)/sqrt(cpos.Mag2()*trkmom.Mag2());
          // the following includes the (Calibrated) light-propagation time delay.  It should eventually be put in the reconstruction FIXME!
          // This velocity should come from conditions FIXME!
          features[3] = tchs.t0().t0()-tchs.time()- std::min((float)200.0,std::max((float)0.0,tchs.hitLen()))*0.005 - _dtoffset;
          // hard cut on the energy difference.  This rejects cosmic rays which hit the calo and produce an upstream-going track that is then
          // reconstructed as a downstream particle associated to this cluster
          if(features[0] < _maxde){
            // evaluate the MVA
            auto mvaout - mva->infer(features.data());
            mvacol->push_back(MVAResult(mvaout[0]));
          }
        }
      }
    }
    // put the output products into the event
    event.put(move(mvacol));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackPID)
