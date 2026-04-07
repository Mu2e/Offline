// Generate two DIOs from a muon stop, with artificial position/stop time shifts
// Orginal author: Michael MacKenzie, 2026, based on DIOGenerator_tool

// Framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// CLHEP
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"

// Offline
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

// ROOT
#include "TH1.h"

// C++
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

namespace mu2e {
  //================================================================
  class DoubleDIOGenerator : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> sim_col       {Name("inputSimParticles")     , Comment("SimParticleCollection with muons stopping in the target")};
      fhicl::Atom<std::string>   material      {Name("stoppingTargetMaterial"), Comment("Material for muon lifetime"), "Al"};
      fhicl::DelegatedParameter  dio_tool      {Name("dioTool")               , Comment("DIO spectrum (and variables) to be generated")};
      fhicl::Atom<double>        tstop_scale   {Name("StopTimeScale")         , Comment("Scale for the stopping time shifting")};
      fhicl::Atom<double>        stop_radius   {Name("StopRadius")            , Comment("Set the radius of the stopping position distribution")};
      fhicl::Atom<double>        stop_x0       {Name("StopX0")                , Comment("Set the x center of the stopping position distribution"), -3904.};
      fhicl::Atom<double>        stop_y0       {Name("StopY0")                , Comment("Set the y center of the stopping position distribution"), 0.};
      fhicl::Atom<bool>          same_stop     {Name("SameStop")              , Comment("Use identical stop positions")};
      fhicl::Atom<bool>          do_histograms {Name("doHistograms")          , Comment("Create histograms"), false};
      fhicl::Atom<unsigned>      verbosity     {Name("verbosity")             , Comment("Verbosity output"), 0};

    };

    // Histogramming info
    struct Hist_t {
      TH1* mom;
      TH1* cz;
      TH1* total_mom;
      TH1* delta_theta;
      TH1* delta_r;
      TH1* delta_r_calo;
      TH1* delta_t;
    };
    enum {kMaxHists = 100};

    using Parameters= art::EDProducer::Table<Config>;
    explicit DoubleDIOGenerator(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    // Inputs
    double   muon_lifetime_;
    art::ProductToken<SimParticleCollection> const sims_token_;
    double   tstop_scale_;
    double   stop_radius_;
    double   stop_x0_;
    double   stop_y0_;
    bool     same_stop_;
    bool     do_histograms_;
    unsigned verbosity_;

    // Fields
    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    CLHEP::RandFlat randFlat_;
    std::unique_ptr<ParticleGeneratorTool> Generator_;
    Hist_t* hists_[kMaxHists];

    // Methods
    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double muon_time, ParticleGeneratorTool* gen);
    void bookHistograms(const int index, const char* title) {
      if(index >= kMaxHists) throw cet::exception("BADCONFIG") << "Too many histogram books! Index = " << index;
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory dir = tfs->mkdir(std::format("hist_{}", index).c_str(), title);

      hists_[index] = new Hist_t;
      auto Hist = hists_[index];

      Hist->mom          = dir.make<TH1F>("mom"         , "Electron momentum; p (MeV/c)"                 , 100,  0.,  110.);
      Hist->cz           = dir.make<TH1F>("cz"          , "Electron cos(#theta_{z}); cos(#theta_{z})"    , 100, -1.,    1.);
      Hist->total_mom    = dir.make<TH1F>("total_mom"   , "Total momentum; p (MeV/c)"                    , 100,  0.,  220.);
      Hist->delta_theta  = dir.make<TH1F>("delta_theta" , "Angle between electrons;#Delta#theta (rad)"   , 100,  0.,    1.);
      Hist->delta_r      = dir.make<TH1F>("delta_r"     , "Origin distance;#Deltar (mm)"                 , 100,  0.,  200.);
      Hist->delta_r_calo = dir.make<TH1F>("delta_r_calo", "Distance projected at calo;#Deltar(calo) (mm)", 100,  0., 1000.);
      Hist->delta_t      = dir.make<TH1F>("delta_t"     , "Origin time difference;#Deltat (ns)"          , 100,  0.,  200.);
    }

    void fillHistograms(const int index, const StageParticleCollection* gen) {
      // Check the inputs are valid
      Hist_t* Hist = hists_[index];
      if(!Hist) throw cet::exception("BADCONFIG") << "Filling a histogram book that wasn't initialized! Index = " << index;
      if(gen->size() != 2) {
        std::cout << "[DoubleDIOGenerator::" << __func__ << "] Output has " << gen->size() << " gen particles --> skipping\n";
        return;
      }

      // Get the gen particles
      const auto& gen_1 = gen->at(0);
      const auto& gen_2 = gen->at(1);

      // Get the gen info
      const auto pos_1 = gen_1.position();
      const auto pos_2 = gen_2.position();
      const auto mom_1 = gen_1.momentum().vect();
      const auto mom_2 = gen_2.momentum().vect();

      // Calculate some relative quanities
      const double p1dotp2 = mom_1.dot(mom_2);
      const double p_1 = mom_1.mag();
      const double p_2 = mom_2.mag();
      const double cz = p1dotp2/(p_1*p_2);

      const double dr = (pos_1 - pos_2).mag();
      const double dt = std::abs(gen_1.time() - gen_2.time());

      const auto dir_1 = mom_1 / p_1;
      const auto dir_2 = mom_2 / p_2;

      const float dz_to_calo = 4500.; // FIXME: Rough distance in Run 1B, for approximate delta R
      const auto pos_calo_1 = pos_1 + dz_to_calo/dir_1.z()*dir_1;
      const auto pos_calo_2 = pos_2 + dz_to_calo/dir_2.z()*dir_2;
      const float dr_calo = (pos_calo_1 - pos_calo_2).mag();

      // Fill the histograms
      Hist->mom          ->Fill(p_1);
      Hist->mom          ->Fill(p_2);
      Hist->cz           ->Fill(mom_1.cosTheta());
      Hist->cz           ->Fill(mom_2.cosTheta());
      Hist->total_mom    ->Fill(p_1 + p_2);
      Hist->delta_theta  ->Fill(cz);
      Hist->delta_r      ->Fill(dr);
      Hist->delta_r_calo ->Fill(dr_calo);
      Hist->delta_t      ->Fill(dt);
    }

  };

  //================================================================
  DoubleDIOGenerator::DoubleDIOGenerator(const Parameters& config)
    : EDProducer{config}
    , muon_lifetime_(GlobalConstantsHandle<PhysicsParams>()->getDecayTime(config().material()))
    , sims_token_(consumes<SimParticleCollection>(config().sim_col()))
    , tstop_scale_(config().tstop_scale())
    , stop_radius_(config().stop_radius())
    , stop_x0_(config().stop_x0())
    , stop_y0_(config().stop_y0())
    , same_stop_(config().same_stop())
    , do_histograms_(config().do_histograms())
    , verbosity_(config().verbosity())
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randExp_(eng_)
    , randFlat_(eng_)
  {
    produces<mu2e::StageParticleCollection>();

    //-----------------------------------------------------------------
    // Initialize the generator
    //-----------------------------------------------------------------
    const auto pset = config().dio_tool.get<fhicl::ParameterSet>();
    Generator_ = art::make_tool<ParticleGeneratorTool>(pset);
    Generator_->finishInitialization(eng_, config().material());

    //-----------------------------------------------------------------
    // Book histograms, if requested
    //-----------------------------------------------------------------
    if(do_histograms_) {
      for(int index = 0; index < kMaxHists; ++index) hists_[index] = nullptr;

      bookHistograms(0, "all events");
    }
  }

  //================================================================
  void DoubleDIOGenerator::produce(art::Event& event) {
    auto output(std::make_unique<StageParticleCollection>());

    const auto simh = event.getValidHandle<SimParticleCollection>(sims_token_);
    const auto mus  = stoppedMuMinusList(simh);

    if(!mus.empty()) {
      const size_t index = randFlat_.fireInt(mus.size());
      const auto& mustop = mus[index];
      const double time  = mustop->endGlobalTime() + randExp_.fire(muon_lifetime_);
      addParticles(output.get(), mustop, time, Generator_.get());
    }

    event.put(std::move(output));
  }

  //================================================================
  void DoubleDIOGenerator::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> mustop,
                            double muon_time,
                            ParticleGeneratorTool* gen)
  {

    bool first = true;
    const auto pos = mustop->endPosition();
    constexpr int nsamples = 2; // sample the DIO spectrum twice
    for(int isample = 0; isample < nsamples; ++isample) {
      auto daughters = gen->generate();
      for(const auto& d: daughters) {
        // For the first daughter use the stop position+time, for the rest offset set this randomly
        const double x = (same_stop_ || first) ? pos.x() : stop_x0_ + stop_radius_*2.*(randFlat_.fire() - 0.5);
        const double y = (same_stop_ || first) ? pos.y() : stop_y0_ + stop_radius_*2.*(randFlat_.fire() - 0.5);
        const double z = pos.z(); // assume the same z location
        const double t = (first) ? muon_time : muon_time + tstop_scale_*2.*(randFlat_.fire() - 0.5);
        if(verbosity_ > 2) std::cout << "  Adding DIO with x = {" << x << ", " << y << ", " << z << "}"
                                     << " t = " << t
                                     << " p = {" << d.fourmom.x() << ", " << d.fourmom.y() << ", " << d.fourmom.z() << "}"
                                     << " E = " << d.fourmom.t()
                                     << std::endl;

        output->emplace_back(mustop,
                             d.creationCode,
                             d.pdgId,
                             CLHEP::Hep3Vector(x, y, z),
                             d.fourmom,
                             t
                             );
      }
      first = false; // shift the next sample's origin/time
    }

    // Optionally fill histograms
    if(do_histograms_) {
      fillHistograms(0, output);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::DoubleDIOGenerator)
