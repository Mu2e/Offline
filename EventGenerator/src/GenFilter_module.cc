// Purpose: Filter events based on generator particles, focused on muon target stop daughters (e.g. DIO)
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
        fhicl::Atom<bool> isNull_{Name("isNull"), Comment("Flag to turn off event filtering"), false};
        fhicl::Atom<bool> maxr_cut_{Name("filterMaxR"), Comment("Filter on the maximum radius of particles in a B-field"), false};
        fhicl::Atom<double> maxr_min{Name("maxr_min"), Comment("Minimum particle radius"), -1.};
        fhicl::Atom<double> maxr_max{Name("maxr_max"), Comment("Maximum particle radius"), 1.e7};
        fhicl::Atom<bool> energy_cut_{Name("filterEnergy"), Comment("Filter on particle energies"), false};
        fhicl::Atom<double> emin{Name("emin"), Comment("Minimum particle energy"), 0.};
        fhicl::Atom<double> emax{Name("emax"), Comment("Maximum particle energy"), 1.e10};
        fhicl::Atom<bool> momentum_cut_{Name("filterMomentum"), Comment("Filter on particle momenta"), false};
        fhicl::Atom<double> pmin{Name("pmin"), Comment("Minimum particle momentum"), 0.};
        fhicl::Atom<double> pmax{Name("pmax"), Comment("Maximum particle momentum"), 1.e10};
        fhicl::Atom<bool> pT_cut_{Name("filterPT"), Comment("Filter on particle transverse momenta"), false};
        fhicl::Atom<double> ptmin{Name("ptmin"), Comment("Minimum transverse momentum"), 0.};
        fhicl::Atom<double> ptmax{Name("ptmax"), Comment("Maximum transverse momentum"), 1.e10};
        fhicl::Atom<bool> makeplots{Name("makeplots"), Comment("Fill an output ntuple of particle information"), false};
      };
      explicit GenFilter(const art::EDFilter::Table<Config>& config);
      virtual bool filter(art::Event& event) override;

    private:
      art::InputTag _SimToken;
      const StageParticleCollection* _SimCol;
      const bool isNull_; //turn off all filtering
      const bool maxr_cut_;
      const double maxr_min_;
      const double maxr_max_;
      const bool energy_cut_;
      const double emin_;
      const double emax_;
      const bool momentum_cut_;
      const double pmin_;
      const double pmax_;
      const bool pT_cut_;
      const double ptmin_;
      const double ptmax_;
      const bool makeplots_;
      TTree* genTree;
      Float_t _maxr;
      Float_t _energy;
      Float_t _mom;
      Float_t _momT;
      Float_t _posT;
      Float_t _cosTheta;
      Float_t _time;
  };

  GenFilter::GenFilter(const art::EDFilter::Table<Config>& config) :
     EDFilter{config}
    , _SimToken(config().SimToken())
    , isNull_{config().isNull_()}
    , maxr_cut_{config().maxr_cut_()}
    , maxr_min_(config().maxr_min())
    , maxr_max_(config().maxr_max())
    , energy_cut_{config().energy_cut_()}
    , emin_(config().emin())
    , emax_(config().emax())
    , momentum_cut_{config().momentum_cut_()}
    , pmin_(config().pmin())
    , pmax_(config().pmax())
    , pT_cut_{config().pT_cut_()}
    , ptmin_(config().ptmin())
    , ptmax_(config().ptmax())
    , makeplots_{config().makeplots()}
  {
    if(makeplots_){
      art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      genTree->Branch("maxr", &_maxr, "maxr/F");
      genTree->Branch("energy", &_energy, "energy/F");
      genTree->Branch("mom", &_mom, "mom/F");
      genTree->Branch("momT", &_momT, "momT/F");
      genTree->Branch("posT", &_posT, "posT/F");
      genTree->Branch("cosTheta", &_cosTheta, "cosTheta/F");
      genTree->Branch("time", &_time, "time/F");
    }
  }

  bool GenFilter::filter(art::Event& event) {
    if(isNull_) return true;
    if(!(maxr_cut_ || energy_cut_ || momentum_cut_ || pT_cut_)) return true;
    auto sim = event.getValidHandle<StageParticleCollection>(_SimToken);
    _SimCol = sim.product();
    bool pass_maxr_cut = !maxr_cut_;
    bool pass_energy_cut = !energy_cut_;
    bool pass_momentum_cut = !momentum_cut_;
    bool pass_pT_cut = !momentum_cut_;
    for(const auto& aParticle : *_SimCol){
      //  make momentum and position vectors
      GeomHandle<DetectorSystem> det;
      ROOT::Math::XYZVectorF pos = XYZVectorF(det->toDetector(aParticle.position()));
      ROOT::Math::XYZTVector pos0(pos.x(), pos.y(), pos.z(), aParticle.time());
      ROOT::Math::PxPyPzMVector mom0(aParticle.momentum().x(), aParticle.momentum().y(), aParticle.momentum().z(), aParticle.momentum().t());
      _energy = aParticle.momentum().e();
      _mom = aParticle.momentum().rho();
      _momT = sqrt(mom0.x()*mom0.x() + mom0.y()*mom0.y());

      if(maxr_cut_) { // cut on maximum R in the solenoid
        // extract charge
        static GlobalConstantsHandle<ParticleDataList> pdt;
        auto charge = pdt->particle(aParticle.pdgId()).charge();

        // check if the particle is charged, else ignore this cut
        if(fabs(charge) > 1.e-5) {
          // extact field
          GeomHandle<BFieldManager> bfmgr;
          mu2e::GeomHandle<mu2e::Tracker> tracker;
          auto tracker_origin = det->toMu2e(tracker->origin());
          ROOT::Math::XYZVector bnom(bfmgr->getBField(tracker_origin));

          // make the loophelix
          KinKal::LoopHelix lh(pos0, mom0, charge, bnom);
          // calculate rmax and add maxr to siminfo
          _maxr = sqrt(lh.cx()*lh.cx()+lh.cy()*lh.cy())+fabs(lh.rad());
          pass_maxr_cut |= _maxr < maxr_max_ && _maxr > maxr_min_;
        } else {  //not defined for neutral particles
          pass_maxr_cut = true; //default to passing for neutral particles
          _maxr = -1.;
        }
      } else {
        _maxr = -1.; //not filtered on, so not calculated
      }

      // particle energy/momentum cuts
      if(energy_cut_) {
        pass_energy_cut |= _energy > emin_ && _energy < emax_;
      }
      if(momentum_cut_) {
        pass_momentum_cut |= _mom > pmin_ && _mom < pmax_;
      }
      if(pT_cut_) {
        pass_pT_cut |= _momT > ptmin_ && _momT < ptmax_;
      }

      if(makeplots_){
        // fill other branches for plots
        _posT = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
        _cosTheta = cos(atan2(_momT,mom0.z()));
        _time = aParticle.time();
        genTree->Fill();
      }
    }
    return pass_maxr_cut && pass_energy_cut && pass_momentum_cut && pass_pT_cut;
  }
}

using mu2e::GenFilter;
DEFINE_ART_MODULE(GenFilter)
