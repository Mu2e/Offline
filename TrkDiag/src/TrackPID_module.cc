//
// Reconstruction-level PID determination from a track.  The current implementation
// just uses the calorimeter matching information (TrkCaloHit), eventually it will
// also process dE/dx, etc.
//
// Original author: Dave Brown (LBNL)
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// utilities
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
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
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        fhicl::Atom<float> MaxDE{Name("MaxDE"), Comment("Maximum difference between calorimeter cluster EDep energy and the track energy (assuming electron mass)")};
        fhicl::Atom<float> DT{ Name("DeltaTOffset"),
          Comment("Track - Calorimeter time offset")}; // this should be a condition FIXME!
        fhicl::Atom<art::InputTag> kalSeedPtrTag{Name("KalSeedPtrCollection"), Comment("Input tag for KalSeedPtrCollection")};
        fhicl::Atom<bool> printMVA{Name("printMVA"), Comment("print the MVA used"), false};
        fhicl::Atom<std::string> datFilename{Name("datFilename"), Comment("Filename for the .dat file to use")};
        fhicl::Atom<int> debug{Name("debugLevel"), Comment("Debug printout Level"), 0};
      };

      using Parameters = art::EDProducer::Table<Config>;
      TrackPID(const Parameters& conf);

    private:
      void produce(art::Event& event) override;
      void initializeMVA(std::string xmlfilename);

      float _maxde, _dtoffset;
      art::InputTag _kalSeedPtrTag;
      bool _printMVA;
      int _debug;

      std::shared_ptr<TMVA_SOFIE_TrackPID::Session> mva_;
  };

  TrackPID::TrackPID(const Parameters& conf) :
    art::EDProducer(conf),
    _maxde(conf().MaxDE()),
    _dtoffset(conf().DT()),
    _kalSeedPtrTag(conf().kalSeedPtrTag()),
    _printMVA(conf().printMVA()),
    _debug(conf().debug())
  {
    produces<MVAResultCollection>();

    ConfigFileLookupPolicy configFile;
    mva_ = std::make_shared<TMVA_SOFIE_TrackPID::Session>(configFile(conf().datFilename()));
  }

  void TrackPID::produce(art::Event& event ) {
    mu2e::GeomHandle<mu2e::Calorimeter> calo;
    // create output
    unique_ptr<MVAResultCollection> mvacol(new MVAResultCollection());
    // get the KalSeedsPtrs
    art::Handle<KalSeedPtrCollection> kalSeedPtrHandle;
    event.getByLabel(_kalSeedPtrTag, kalSeedPtrHandle);
    const auto& kalSeedPtrs = *kalSeedPtrHandle;

    // Go through the tracks and calculate the track PID
    for (const auto& kalSeedPtr : kalSeedPtrs) {
      const auto& kalSeed = *kalSeedPtr;
      std::array<float,4> features{-9999,-9999,-9999,-9999}; // features used for training
      double mvaval = -1;

      // Fill the features
      static TrkFitFlag goodfit(TrkFitFlag::kalmanOK);
      if (kalSeed.status().hasAllProperties(goodfit)){
        if(kalSeed.hasCaloCluster() && kalSeed.caloHit()._flag.hasAllProperties(StrawHitFlag::active)){
          auto const& tchs = kalSeed.caloHit();
          auto const& cc = tchs.caloCluster();
          XYZVectorD trkmom = kalSeed.nearestSegment(tchs._rptoca)->momentum3();
          features[0] = cc->energyDep() - sqrt(trkmom.Mag2());
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
            auto mvaout = mva_->infer(features.data());
            mvaval = mvaout[0];
          }
          else if (_debug > 0) {
            printf("energy difference at %f, above threshold %f", features[0], _maxde);
          }
        }
      }
      if (_debug > 0) {
        printf("TrackPID ; input features: %f ; %f ; %f ; %f ; output: %f",
          features[0], features[1], features[2], features[3], mvaval);
      }
      mvacol->push_back(MVAResult(mvaval));
    }
    if (mvacol->size() != kalSeedPtrs.size()) {
      throw cet::exception("TrackPID") << "KalSeedPtr and MVAResult sizes are inconsistent: KalSeedPtr.size() = " << kalSeedPtrs.size() << " ; MVAResult.size() = " << mvacol->size();
    }
    // put the output products into the event
    event.put(move(mvacol));
  }
}// mu2e

DEFINE_ART_MODULE(mu2e::TrackPID)
