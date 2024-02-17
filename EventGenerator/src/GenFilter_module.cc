// Purpose: Event filter for DIO simulations
// author: S Middleton  2024
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include <iostream>
#include <string>

#include "TTree.h"
using namespace std;
namespace mu2e {

  class GenFilter : public art::EDFilter {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> SimToken{Name("StageParticleCollection"),Comment("")};
        fhicl::Atom<art::InputTag>SimTag{Name("StageParticleCollection"),Comment("SimTag")};
        fhicl::Atom<double> maxr_min{Name("maxr_min"),0};
        fhicl::Atom<double> maxr_max{Name("maxr_max"),1e7};
        fhicl::Atom<bool> makeplots{Name("makeplots"),false};
        fhicl::Atom<bool> isNull{Name("isNull"),true};
      };
      explicit GenFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;

    private:
      art::InputTag _SimToken;
      const StageParticleCollection* _SimCol;
      double maxr_min_;
      double maxr_max_;
      bool makeplots_;
      bool isNull_;
      TTree* genTree;
      Float_t _maxr;
      Float_t _momT;
      Float_t _posT;
      Float_t _cosTheta;
      Float_t _time;
  };

  GenFilter::GenFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , _SimToken(config().SimToken())
    , maxr_min_(config().maxr_min())
    , maxr_max_(config().maxr_max())
    , makeplots_{config().makeplots()}
    , isNull_{config().isNull()}
  {
    if(makeplots_){
      art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("maxr", &_maxr, "maxr/F");
      genTree->Branch("momT", &_momT, "momT/F");
      genTree->Branch("posT", &_posT, "posT/F");
      genTree->Branch("cosTheta", &_cosTheta, "cosTheta/F");
      genTree->Branch("time", &_time, "time/F");
    }
  }

  bool GenFilter::filter(art::Event& event) {
    if(isNull_) return true;
    bool passed = false;
    auto sim = event.getValidHandle<StageParticleCollection>(_SimToken);
    _SimCol = sim.product();
    for(const auto& aParticle : *_SimCol){
    //  make momentum and position vectors
      GeomHandle<DetectorSystem> det;
      ROOT::Math::XYZVectorF pos = XYZVectorF(det->toDetector(aParticle.position()));
      ROOT::Math::XYZTVector pos0(pos.x(), pos.y(), pos.z(), aParticle.time());
      ROOT::Math::PxPyPzMVector mom0(aParticle.momentum().x(), aParticle.momentum().y(), aParticle.momentum().z(), aParticle.momentum().t());

      // extract charge
      static GlobalConstantsHandle<ParticleDataList> pdt;
      auto charge = pdt->particle(aParticle.pdgId()).charge();

      // extact field
      GeomHandle<BFieldManager> bfmgr;
      mu2e::GeomHandle<mu2e::Tracker> tracker;
      auto tracker_origin = det->toMu2e(tracker->origin());
      //XYZVectorF pos3Vec = XYZVectorF(aParticle.position().x(),aParticle.position().y(),aParticle.position().z());
      ROOT::Math::XYZVector bnom(bfmgr->getBField(tracker_origin));

      // make the loophelix
      KinKal::LoopHelix lh(pos0, mom0, charge, bnom);
      // calculate rmax and add maxr to siminfo
      _maxr =sqrt(lh.cx()*lh.cx()+lh.cy()*lh.cy())+fabs(lh.rad());
      if(makeplots_){
        // fill other branches for plots
        _momT =  sqrt(mom0.x()*mom0.x() + mom0.y()*mom0.y());
        _posT = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
        _cosTheta = cos(atan2(_momT,mom0.z()));
        _time = aParticle.time();
        genTree->Fill();
      }
      if((_maxr < maxr_max_ and _maxr > maxr_min_ )){ passed = true; }
    }
    return passed;
  }
}

using mu2e::GenFilter;
DEFINE_ART_MODULE(GenFilter)
