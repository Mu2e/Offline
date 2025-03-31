// Generate antiproton events upstream of the stopping target using rough approximations
// Michael MacKenize, 2024


#include <iostream>
#include <string>
#include <cmath>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/TargetFoil.hh"

#include "TH1.h"

namespace mu2e {

  //================================================================
  class SimpleAntiprotonGun : public art::EDProducer {
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> t0{Name("t0"), Comment("Start time for the antiproton distribution"), 400.}; //flat time spectrum
      fhicl::Atom<double> tmax{Name("tmax"), Comment("Maximum time for the antiproton distribution"), 2000.};
      fhicl::Atom<double> pzmax{Name("pzmax"), Comment("Antiproton maximum generated z momentum (flat spectrum)"), 100.}; //flat momentum spectrum
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbosity level"), 0};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of the conversion kinematics"), false};
    };
    struct Stop_t {
      float x;
      float y;
      float z;
      float t;
      Stop_t(): x(0.), y(0.), z(0.), t(0.) {}
      Stop_t(float x_in, float y_in, float z_in, float t_in) : x(x_in), y(y_in), z(z_in), t(t_in) {}
    };

    PDGCode::type       pdgId_;
    GenId               genId_;
    double              mass_;

    double              t0_;
    double              tmax_;
    double              pzmax_;

    int                 verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;

    double _foilY; //foil origin
    double _foilX;
    double _foilZ;
    double _foilRIn; //foil radial dimensions
    double _foilROut;

    //-----------------------------------------------------------------------------
    // histogramming
    //-----------------------------------------------------------------------------
    bool    makeHistograms_;
    TH1*   _hX;
    TH1*   _hY;
    TH1*   _hZ;
    TH1*   _hR;
    TH1*   _hTime;
    TH1*   _hPz;

  private:

  public:
    using Parameters= art::EDProducer::Table<Config>;
    explicit SimpleAntiprotonGun(const Parameters& conf);

    virtual void beginRun(art::Run& run);
    virtual void produce(art::Event& event);
    Stop_t generateStop();
  };

  //================================================================
  SimpleAntiprotonGun::SimpleAntiprotonGun(const Parameters& conf)
    : EDProducer{conf}
    , pdgId_(PDGCodeDetail::anti_proton)
    , genId_(GenId::antiproton)
    , mass_(GlobalConstantsHandle<ParticleDataList>()->particle(pdgId_).mass())
    , t0_(conf().t0())
    , tmax_(conf().tmax())
    , pzmax_(conf().pzmax())
    , verbosity_(conf().verbosity())
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
    , makeHistograms_(conf().makeHistograms())
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosity_ > 0) {
      std::cout<<"SimpleAntiprotonGun: using gen ID " << genId_ << std::endl;
    }

    if(t0_ > tmax_) {
      throw cet::exception("BADCONFIG") << "Time range is unphysical: t0 = " << t0_ << " tmax = " << tmax_ << "\n";
    }
    if(pzmax_ < 0.) {
      throw cet::exception("BADCONFIG") << "Momentum range is unphysical: pz_max = " << pzmax_ << "\n";
    }

    if(makeHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;
      _hX      = tfs->make<TH1F>("hX"     , "X"           ,  100,     0.,  100.);
      _hY      = tfs->make<TH1F>("hY"     , "Y"           ,  100,     0.,  100.);
      _hZ      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400., 6400.);
      _hR      = tfs->make<TH1F>("hR"     , "R"           ,  100,     0.,  100.);
      _hTime   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,     0., std::max(tmax_, 2000.));
      _hPz     = tfs->make<TH1F>("hPz"    , "Pz"          ,  100,     0.,  std::max(pzmax_, 100.));
    }
  }

  //================================================================
  void SimpleAntiprotonGun::beginRun(art::Run&) {

    // Initialize the target foil position

    // Find the foil that's first in z (should be the first index)
    double z_first = 0.;
    GeomHandle<StoppingTarget> target;
    for(int ifoil = 0; ifoil < target->nFoils(); ++ifoil) {
      auto foil = target->foil(ifoil);
      auto origin = foil.centerInMu2e();
      if(ifoil == 0 || origin.z() < z_first) {
        z_first = origin.z();
        _foilX = origin.x();
        _foilY = origin.y();
        _foilZ = origin.z() - 2.*foil.halfThickness(); //ensure it's upstream of the ST with a small gap
        _foilRIn = foil.rIn();
        _foilROut = foil.rOut();
      }
      if(verbosity_ > 0) printf("SimpleAntiprotonGun::%s: Foil %2i has (z, r0, r1) = (%8.1f, %4.1f, %4.1f), lowest z = %.1f\n", __func__, ifoil, origin.z(), foil.rIn(), foil.rOut(), _foilZ);
    }
  }

  //================================================================
  SimpleAntiprotonGun::Stop_t SimpleAntiprotonGun::generateStop() {
    Stop_t stop;
    stop.z = _foilZ; //just upstream of the foil
    const double r = _foilRIn + (_foilROut - _foilRIn)*randFlat_.fire(); //flat radially
    const double phi =  CLHEP::twopi*randFlat_.fire(); //flat in phi
    stop.x = r*std::cos(phi) + _foilX;
    stop.y = r*std::sin(phi) + _foilY;

    //assume an exponential decay in time
    stop.t = t0_ + (tmax_ - t0_)*randFlat_.fire();

    return stop;
  }

  //================================================================
  void SimpleAntiprotonGun::produce(art::Event& event) {
    auto output{std::make_unique<GenParticleCollection>()};


    // Generate the starting position and time
    const auto start = generateStop();

    const CLHEP::Hep3Vector pos(start.x, start.y, start.z);

    // start at a stop
    const double energy = mass_;
    CLHEP::Hep3Vector p3(0., 0., pzmax_*randFlat_.fire()); //moving downstream only
    CLHEP::HepLorentzVector fourmom(p3, energy);
    output->emplace_back(GenParticle(pdgId_,
                                     genId_,
                                     pos,
                                     fourmom,
                                     start.t));

    event.put(std::move(output));

    if(verbosity_ > 1) {
      printf("SimpleAntiprotonGun::%s: Add antiproton with (x,y,z,t) = (%7.1f, %5.1f, %8.1f, %6.1f)\n",
             __func__, start.x, start.y, start.z, start.t);
    }
    //-----------------------------------------------------------------------------
    // if requested, fill histograms
    //-----------------------------------------------------------------------------
    if(makeHistograms_) {
      _hX->Fill(pos.x() - _foilX);
      _hY->Fill(pos.y() - _foilY);
      _hZ->Fill(pos.z());
      _hTime->Fill(start.t);
      const double r = std::sqrt(std::pow(pos.x() - _foilX, 2) + std::pow(pos.y() - _foilY, 2));
      _hR->Fill(r);
      _hPz->Fill(p3.z());
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimpleAntiprotonGun)
