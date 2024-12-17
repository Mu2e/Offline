// Generate antiproton events in the stopping target using rough approximation
// Michael MacKenize, 2024


#include <iostream>
#include <string>
#include <cmath>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
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
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
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
      fhicl::Atom<double> foilSurvivalRate{Name("foilSurvivalRate"), Comment("Rate of antiprotons surviving each foil"), 0.01};
      fhicl::Atom<double> timeArrivalRate{Name("timeArrivalRate"), Comment("Exponential rate for antiproton arrival"), 130.};
      fhicl::Atom<double> t0{Name("t0"), Comment("Start time for the antiproton distribution"), 600.}; //time distribution falls again below 600 ns
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"), Comment("Input sim particle collection")};
      fhicl::Atom<bool> ignoreInput{Name("ignoreInput"), Comment("Ignore the input sim particle information, using random stops"), false};
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
    ProcessCode         processCode_;
    double              mass_;

    double              foilSurvivalRate_;
    double              timeArrivalRate_;
    double              t0_;
    art::ProductToken<SimParticleCollection> const simsToken_;
    bool                ignoreInput_;

    int                 verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    CLHEP::RandFlat randFlat_;

    std::vector<double> _foilYs; //foil origins
    std::vector<double> _foilXs;
    std::vector<double> _foilZs;
    std::vector<double> _foilRIns; //foil radial dimensions
    std::vector<double> _foilROuts;

    //-----------------------------------------------------------------------------
    // histogramming
    //-----------------------------------------------------------------------------
    bool    makeHistograms_;
    TH1*   _hX;
    TH1*   _hY;
    TH1*   _hZ;
    TH1*   _hR;
    TH1*   _hTime;

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
    , processCode_(ProcessCode::mu2eAntiproton)
    , mass_(GlobalConstantsHandle<ParticleDataList>()->particle(pdgId_).mass())
    , foilSurvivalRate_(conf().foilSurvivalRate())
    , timeArrivalRate_(conf().timeArrivalRate())
    , t0_(conf().t0())
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , ignoreInput_{conf().ignoreInput()}
    , verbosity_(conf().verbosity())
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randFlat_{eng_}
    , makeHistograms_(conf().makeHistograms())
  {
    produces<mu2e::StageParticleCollection>();

    if(verbosity_ > 0) {
      std::cout<<"SimpleAntiprotonGun: using process code " << processCode_ << std::endl;
    }


    if(makeHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;
      _hX      = tfs->make<TH1F>("hX"     , "X"           ,  100,     0.,  100.);
      _hY      = tfs->make<TH1F>("hY"     , "Y"           ,  100,     0.,  100.);
      _hZ      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400., 6400.);
      _hR      = tfs->make<TH1F>("hR"     , "R"           ,  100,     0.,  100.);
      _hTime   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,     0., 2000.);
    }
  }

  //================================================================
  void SimpleAntiprotonGun::beginRun(art::Run&) {

    // Initialize the target foil positions
    if(_foilZs.size() == 0) { //only do this once per job
      GeomHandle<StoppingTarget> target;
      for(int ifoil = 0; ifoil < target->nFoils(); ++ifoil) {
        auto foil = target->foil(ifoil);
        auto origin = foil.centerInMu2e();
        _foilXs.push_back(origin.x());
        _foilYs.push_back(origin.y());
        _foilZs.push_back(origin.z());
        _foilRIns.push_back(foil.rIn());
        _foilROuts.push_back(foil.rOut());
      }
    }
  }

  //================================================================
  SimpleAntiprotonGun::Stop_t SimpleAntiprotonGun::generateStop() {
    // select a stopping target foil using a wrapped exponential distribution
    const int nfoils = _foilZs.size();
    const int foil = int(randExp_.fire(foilSurvivalRate_)*nfoils) % nfoils;

    Stop_t stop;
    stop.z = _foilZs[foil]; //place in the center of the foil as an approximation
    const double r = _foilRIns[foil] + (_foilROuts[foil] - _foilRIns[foil])*randFlat_.fire(); //flat radially
    const double phi =  CLHEP::twopi*randFlat_.fire(); //flat in phi
    stop.x = r*std::cos(phi) + _foilXs[foil];
    stop.y = r*std::sin(phi) + _foilYs[foil];

    //assume an exponential decay in time
    stop.t = t0_ + randExp_.fire(timeArrivalRate_);

    return stop;
  }

  //================================================================
  void SimpleAntiprotonGun::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    // Input sim particle collection
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto& sims = *simh;

    if(sims.begin() == sims.end())
      throw cet::exception("BADINPUT") << "SimpleAntiprotonGun::" << __func__ << ": No input sim particles found\n";

    // Only take the first sim particle (only one per event is expected)
    const art::Ptr<SimParticle> sim(simh, sims.begin()->first.asInt());

    // Either use the input antiproton event or ignore it and create an independent stop position
    const auto stop = (ignoreInput_) ? generateStop() : Stop_t(sim->endPosition().x(), sim->endPosition().y(),
                                                               sim->endPosition().z(), sim->endGlobalTime());

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    // start at a stop
    const double energy = mass_;
    CLHEP::Hep3Vector p3(0., 0., 0.);
    CLHEP::HepLorentzVector fourmom(p3, energy);
    output->emplace_back(sim,
                         processCode_,
                         pdgId_,
                         pos,
                         fourmom,
                         stop.t);

    event.put(std::move(output));

    //-----------------------------------------------------------------------------
    // if requested, fill histograms
    //-----------------------------------------------------------------------------
    if(makeHistograms_) {
      _hX->Fill(pos.x() - _foilXs[0]);
      _hY->Fill(pos.y() - _foilYs[0]);
      _hZ->Fill(pos.z());
      _hTime->Fill(stop.t);
      const double r = std::sqrt(std::pow(pos.x() - _foilXs[0], 2) + std::pow(pos.y() - _foilYs[0], 2));
      _hR->Fill(r);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimpleAntiprotonGun)
