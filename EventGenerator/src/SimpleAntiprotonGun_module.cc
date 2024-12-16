// Generate antiproton events in the stopping target using rough approximation
// Michael MacKenize, 2024


#include <iostream>
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
      fhicl::Atom<double> foilSurvivalRate{Name("foilSurvivalRate"), Comment("Rate of antiprotons surviving each foil"), 0.01};
      fhicl::Atom<double> timeArrivalRate{Name("timeArrivalRate"), Comment("Exponential rate for antiproton arrival"), 130.};
      fhicl::Atom<double> t0{Name("t0"), Comment("Start time for the antiproton distribution"), 600.}; //time distribution falls again below 600 ns
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbosity level"), 0};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of the conversion kinematics"), false};
    };
    struct Stop_t {
      float x;
      float y;
      float z;
      float t;
    };


    PDGCode::type       pdgId_;
    GenId               genId_;
    double              mass_;

    double              foilSurvivalRate_;
    double              timeArrivalRate_;
    double              t0_;

    int               verbosity_;

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
    , genId_(GenId::antiproton)
    , mass_(GlobalConstantsHandle<ParticleDataList>()->particle(pdgId_).mass())
    , foilSurvivalRate_(conf().foilSurvivalRate())
    , timeArrivalRate_(conf().timeArrivalRate())
    , t0_(conf().t0())
    , verbosity_(conf().verbosity())
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randFlat_{eng_}
    , makeHistograms_(conf().makeHistograms())
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosity_ > 0) {
      std::cout<<"SimpleAntiprotonGun: using GenId = " << genId_ << std::endl;

      std::cout<<"SimpleAntiprotonGun: producing particle "<< pdgId_ << ", mass = "<< mass_ << std::endl;
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
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = generateStop();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    // start at a stop
    const double energy = mass_;
    CLHEP::Hep3Vector p3(0., 0., 0.);
    CLHEP::HepLorentzVector fourmom(p3, energy);
    output->emplace_back(pdgId_,
                         genId_,
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
