// EdepAna: simple analyzer for primary particle energy losses and calo/tracker energy deposition
// Original author: Michael MacKenzie, 2026

// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// Offline
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TString.h"
#include "TF1.h"
#include "TFitResult.h"

// c++
#include <format>
#include <string>
#include <vector>
#include <iostream>

namespace mu2e {

  class EdepAna : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> primaryTag     { Name("primary")                 , Comment("PrimaryParticle tag") };
      fhicl::Atom<art::InputTag> caloShowerTag  { Name("CaloShowerStepCollection"), Comment("CaloShowerStep collection") };
      fhicl::Atom<art::InputTag> strawGasStepTag{ Name("StrawGasStepCollection")  , Comment("StrawGasStep collection") };
      fhicl::Atom<art::InputTag> stepPointTag   { Name("StepPointMCCollection")   , Comment("StepPointMC collection") };
      fhicl::Atom<int>           debugLevel     { Name("debugLevel")              , Comment("Debug level"), 0 };
    };

    // Histograms
    struct Hist_t {
      TH1* h_primary_energy_;
      TH1* h_primary_pdg_;
      TH1* h_primary_start_z_;
      TH1* h_primary_start_r_;
      TH1* h_total_calo_energy_;
      TH1* h_step_energy_; // individual CaloShowerStep energies
      TH2* h_step_energy_vs_time_; // step energy vs time
      TH1* h_trk_front_energy_;
      TH1* h_trk_front_energy_diff_;
      TH2* h_primary_vs_edep_;
      TH1* h_primary_edep_;
      TH1* h_primary_energy_edep_diff_;
      TH1* h_trk_front_energy_edep_diff_;
    };
    enum {kMaxHists = 100};

    // Event info
    struct Info_t {
      const SimParticle* primsim = nullptr;
      const StepPointMC* front_trk_sp = nullptr;
      double primsim_edep = 0.;
      double calo_total_edep = 0.;
      double weight_ = 1.;
    };

    // Output tree information
    enum {kMaxPrimaries = 10};
    struct Tree_t {
      int   nprim;
      float prim_start_x         [kMaxPrimaries];
      float prim_start_y         [kMaxPrimaries];
      float prim_start_z         [kMaxPrimaries];
      float prim_start_px        [kMaxPrimaries];
      float prim_start_py        [kMaxPrimaries];
      float prim_start_pz        [kMaxPrimaries];
      float prim_start_e         [kMaxPrimaries];
      float prim_start_m         [kMaxPrimaries];
      int   prim_start_pdg       [kMaxPrimaries];
      float prim_calo_edep       [kMaxPrimaries];
      float prim_trk_front_energy[kMaxPrimaries];
      float event_calo_edep;
      float event_trk_edep;
      float weight;
      int   run;
      int   subrun;
      int   event;
      Long64_t ngen; // running N(gen) count
    };
    Tree_t data_;
    TTree* tree_ = nullptr;

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit EdepAna(const Parameters& conf);
    virtual void analyze(const art::Event& event) override;
    virtual void endJob() override;
    virtual void beginSubRun(const art::SubRun& subrun);

    void bookHistograms(const int index, const char* title);
    void fillHistograms(Hist_t* Hist, const double Weight = 1.);

    bool isDescendant(const SimParticle* ancestor, const SimParticle* descendant) const {
      if(!ancestor || !descendant) return false;
      if(ancestor->id() == descendant->id()) return true; // same particle
      auto parent = descendant->parent();
      while(parent.isNonnull()) {
        if(parent->id() == ancestor->id()) return true;
        parent = parent->parent();
      }
      return false;
    }

    double edepBySim(const SimParticle* sim, const bool check_parents) const {
      if(!sim) return 0.;
      if(!shower_col_) return 0.;
      double edep = 0.;
      for(const auto& shower : *shower_col_) {
        if(shower.simParticle()->id() == sim->id() || (check_parents && isDescendant(sim, &*shower.simParticle()))) {
          edep += shower.energyDepBirks();
        }
      }
      return edep;
    }

    const StepPointMC* simTrkFrontStep(const SimParticle* sim) {
      if(!sim || !step_point_col_) return nullptr;
      const StepPointMC* front_trk_sp = nullptr;
      for(const auto& sp : *step_point_col_) {
        if(sp.simParticle()->id() == sim->id() &&
           (sp.volumeId() == VirtualDetectorId::TT_FrontHollow || sp.volumeId() == VirtualDetectorId::TT_FrontPA)) {
          if(!front_trk_sp || sp.time() < front_trk_sp->time()) front_trk_sp = &sp;
          break;
        }
      }
      return front_trk_sp;
    }

  //------------------------------------------------------------------------------------------------------------
  // Landau core with power-law tails
  static double landau(double x, double mean, double a, double b) {

    // Evaluate the function
    const double dx = x-mean;
    double val = exp(a*(b*dx-exp(b*dx)));
    return val;
  }

  //------------------------------------------------------------------------------------------------------------
  // Landau core with power-law tails
  static double landau_crystal_ball(double x, double mean, double a, double b, double alpha1, double alpha2, double n1, double n2) {
    // See: https://github.com/pavel1murat/murat/blob/main/scripts/fit_cb4.C

    // Requirements for normalization are N1 and N2 are > 1:
    if(n1 < 1. || n2 < 1.) return 0.;

    // Evaluate the function
    const double dx = x-mean;
    double val = 0.;
    if (dx < -alpha1) { // Low tail
      const double B1 = -alpha1+n1/(a*b*(1-exp(-b*alpha1)));
      const double A1 = exp(a*(-b*alpha1-exp(-b*alpha1)))*pow(B1+alpha1,n1);
      val  = A1/pow(B1-dx,n1);
    } else if (dx < alpha2) { // Landau core
      val = exp(a*(b*dx-exp(b*dx)));
    } else { // High tail
      const double B2 = -alpha2-n2/(a*b*(1. - exp(b*alpha2)));
      const double A2 = exp(a*(b*alpha2-exp(b*alpha2)))*pow(B2+alpha2,n2);
      val  = A2/pow(B2+dx,n2);
    }

    return val;
  }

  //------------------------------------------------------------------------------------------------------------
  static double landau_crystal_ball_func(double* X, double* P) {
    return P[0]*landau_crystal_ball(X[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7]);
  }
  //------------------------------------------------------------------------------------------------------------
  static double landau_func(double* X, double* P) {
    return P[0]*landau(X[0], P[1], P[2], P[3]);
  }

  //------------------------------------------------------------------------------------------------------------
  static void get_landau_seed(TH1* h, double& mean_seed, double& fwhm_seed) {
    if(!h) {
      mean_seed = 0.;
      fwhm_seed = 0.;
      return;
    }
    mean_seed = h->GetBinCenter(h->GetMaximumBin());
    const double h_height = h->GetMaximum();
    const double x1_seed = h->GetBinCenter(h->FindFirstBinAbove(h_height/2.));
    const double x2_seed = h->GetBinCenter(h->FindLastBinAbove(h_height/2.));
    fwhm_seed = x2_seed - x1_seed;
  }

  //------------------------------------------------------------------------------------------------------------
  static int fit_landau(TH1* h, double& mean, double& fwhm) {
    mean = 0.; fwhm = 0.;
    if(!h) return -1;
    if(h->GetEntries() < 100) return -1; // not enough stats for a fit
    double mean_seed, fwhm_seed;
    get_landau_seed(h, mean_seed, fwhm_seed);
    TF1* f = new TF1("landau", landau_func, mean_seed - fwhm_seed, std::min(0., mean_seed + fwhm_seed), 4);
    f->SetParNames("Norm", "#mu", "a", "b");
    f->SetParLimits(1, -50., 0.); // mean
    f->FixParameter(2, 0.5); // a
    f->SetParLimits(3, 0.01, 100.); // b
    f->SetParameters(h->GetMaximum(), mean_seed, 0.5, fwhm_seed/5.);
    auto fit_res = h->Fit(f, "SR");
    mean = f->GetParameter(1);
    const double max_val = f->GetMaximum();
    const double half_max = max_val / 2.;
    if(!fit_res->IsValid()) {
      std::cerr << "Fit failed with code " << fit_res->Status() << std::endl;
    } else {
      try {
        f->SetRange(-100., 100.); // ensure the fit function covers the full histogram range for mean/FWHM calculation
        const double x_left = f->GetX(half_max, mean - 100., mean);
        const double x_right = std::min(f->GetX(half_max, mean, mean + 100.), 0.); // don't include above 0 values
        fwhm = x_right - x_left;
      } catch(...) {
        std::cerr << "Error calculating FWHM from fit\n";
      }
    }
    delete f;
    return fit_res->Status();
  }

  private:
    art::InputTag primary_tag_;
    art::InputTag calo_shower_tag_;
    art::InputTag straw_gas_tag_;
    art::InputTag step_point_tag_;
    int debug_level_;
    std::unique_ptr<StopWatch> watch_ = std::make_unique<StopWatch>();

    Hist_t* hists_[kMaxHists];
    Info_t info_;
    const PrimaryParticle*                 primary_        = nullptr;
    const CaloShowerStepCollection*        shower_col_     = nullptr;
    const StrawGasStepCollection*          straw_step_col_ = nullptr;
    const StepPointMCCollection*           step_point_col_ = nullptr;
    double total_events_ = 0.;
    double events_above_50_mev_ = 0.;
    double total_calo_edep_ = 0.;
    double total_tracker_edep_ = 0.;
    unsigned long ngen_ = 0;
  };


  EdepAna::EdepAna(const Parameters& conf)
    : art::EDAnalyzer{conf}
    , primary_tag_(conf().primaryTag())
    , calo_shower_tag_(conf().caloShowerTag())
    , straw_gas_tag_(conf().strawGasStepTag())
    , step_point_tag_(conf().stepPointTag())
    , debug_level_(conf().debugLevel())
  {
    // register products consumed
    consumes<PrimaryParticle>(primary_tag_);
    consumes<CaloShowerStepCollection>(calo_shower_tag_);
    consumes<StrawGasStepCollection>(straw_gas_tag_);
    consumes<StepPointMCCollection>(step_point_tag_);

    for(int i = 0; i < kMaxHists; ++i) hists_[i] = nullptr;
    bookHistograms(0, "all events");
    bookHistograms(1, "edep 1 MeV");
    bookHistograms(2, "edep 10 MeV");
    bookHistograms(3, "edep 50 MeV");

    // Book the output TTree
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("tree", "Energy deposition data");
    tree_->Branch("nprimaries"              , &data_.nprim);
    tree_->Branch("primary_start_x"         , data_.prim_start_x         , "primary_start_x[nprimaries]/F");
    tree_->Branch("primary_start_y"         , data_.prim_start_y         , "primary_start_y[nprimaries]/F");
    tree_->Branch("primary_start_z"         , data_.prim_start_z         , "primary_start_z[nprimaries]/F");
    tree_->Branch("primary_start_px"        , data_.prim_start_px        , "primary_start_px[nprimaries]/F");
    tree_->Branch("primary_start_py"        , data_.prim_start_py        , "primary_start_py[nprimaries]/F");
    tree_->Branch("primary_start_pz"        , data_.prim_start_pz        , "primary_start_pz[nprimaries]/F");
    tree_->Branch("primary_start_e"         , data_.prim_start_e         , "primary_start_e[nprimaries]/F");
    tree_->Branch("primary_start_m"         , data_.prim_start_m         , "primary_start_m[nprimaries]/F");
    tree_->Branch("primary_start_pdg"       , data_.prim_start_pdg       , "primary_start_pdg[nprimaries]/I");
    tree_->Branch("primary_calo_edep"       , data_.prim_calo_edep       , "primary_calo_edep[nprimaries]/F");
    tree_->Branch("primary_trk_front_energy", data_.prim_trk_front_energy, "primary_trk_front_energy[nprimaries]/F");
    tree_->Branch("event_calo_edep"         , &data_.event_calo_edep);
    tree_->Branch("event_trk_edep"          , &data_.event_trk_edep);
    tree_->Branch("run"                     , &data_.run);
    tree_->Branch("subrun"                  , &data_.subrun);
    tree_->Branch("event"                   , &data_.event);
    tree_->Branch("weight"                  , &data_.weight);
    tree_->Branch("ngen"                    , &data_.ngen);

    if(debug_level_ > 0) watch_->Calibrate();
  }

  //--------------------------------------------------------------------------------------
  void EdepAna::beginSubRun(const art::SubRun& subrun) {
    // Get the generator counter
    auto genCounterHandle = subrun.getHandle<GenEventCount>("genCounter");
    if(genCounterHandle.isValid()) {
      ngen_ += genCounterHandle->count();
    } else {
      std::cerr << "Warning: GenEventCount not found in subrun" << std::endl;
    }
    data_.ngen = ngen_; // update the output tree data
  }

  void EdepAna::bookHistograms(const int index, const char* title) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    hists_[index] = new Hist_t;
    auto Hist = hists_[index];

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("hist_{}", index), title);

   Hist->h_primary_energy_     = dir.make<TH1F>("primary_energy"    , "Primary energy;Energy (MeV)"              , 150,    0.,  150.);
   Hist->h_primary_pdg_        = dir.make<TH1D>("primary_pdg"       , "Primary PDG ID;PDG ID"                    ,  50,  -25.,   25.);
   Hist->h_primary_start_z_    = dir.make<TH1D>("primary_start_z"   , "Primary start z;Start z (mm)"             , 500, 3000., 8000.);
   Hist->h_primary_start_r_    = dir.make<TH1D>("primary_start_r"   , "Primary start radius;Start radius (mm)"   , 100,    0.,  200.);
   Hist->h_total_calo_energy_  = dir.make<TH1F>("total_calo_energy" , "Total Calo energy from steps;Energy (MeV)", 200,    0.,  200.);
   Hist->h_step_energy_        = dir.make<TH1F>("calo_step_energy"  , "CaloShowerStep energy;E_{step} (MeV)"     , 100,    0.,  100.);
   Hist->h_trk_front_energy_   = dir.make<TH1F>("trk_front_energy"  , "Energy of StepPointMC at front of tracker;Energy (MeV)", 1500, 0., 150.);
   Hist->h_trk_front_energy_diff_ = dir.make<TH1F>("trk_front_energy_diff", "Energy difference of StepPointMC at front of tracker and primary;Energy (MeV)", 500, -100., 0.);
   Hist->h_primary_energy_edep_diff_ = dir.make<TH1F>("primary_energy_edep_diff", "Primary Edep - energy;Energy (MeV)", 300, -150., 0.);
   Hist->h_primary_edep_ = dir.make<TH1F>("primary_edep", "Primary Edep;Energy (MeV)", 300, 0., 150.);
   Hist->h_primary_vs_edep_ = dir.make<TH2F>("primary_vs_edep", "Primary energy vs Edep;Primary energy (MeV);Edep (MeV)", 150, 0., 150., 150, 0., 150.);
   Hist->h_step_energy_vs_time_ = dir.make<TH2F>("calo_step_energy_vs_time", "CaloShowerStep energy vs time;Time (ns);Energy (MeV)", 40, 0., 2000., 100, 0., 100.);
   Hist->h_trk_front_energy_edep_diff_ = dir.make<TH1F>("trk_front_energy_edep_diff", "Energy difference of StepPointMC at front of tracker and primary Edep;Energy (MeV)", 500, -50., 0.);
  }

  void EdepAna::fillHistograms(Hist_t* Hist, const double Weight) {
    if(debug_level_ > 0) watch_->SetTime(__func__);
    if(!Hist) throw std::runtime_error("Uninitialized histogram book!");

    const SimParticle* primsim = info_.primsim;
    const StepPointMC* front_trk_sp = info_.front_trk_sp;
    const double primsim_edep = info_.primsim_edep;


    if(primsim) {
      Hist->h_primary_energy_->Fill(primsim->startMomentum().e(), Weight);
      Hist->h_primary_pdg_->Fill(primsim->pdgId(), Weight);
      if(debug_level_ > 1) std::cout << "EdepAna: primary energy " << primsim->startMomentum().e() << " PDG " << primsim->pdgId() << std::endl;
      Hist->h_primary_start_z_->Fill(primsim->startPosition().z(), Weight);
      Hist->h_primary_start_r_->Fill(std::sqrt(std::pow(primsim->startPosition().x()+3904.,2) + std::pow(primsim->startPosition().y(),2)), Weight);
      const double edep = primsim_edep;
      Hist->h_primary_edep_->Fill(edep, Weight);
      Hist->h_primary_energy_edep_diff_->Fill(edep - primsim->startMomentum().e(), Weight);
      Hist->h_primary_vs_edep_->Fill(primsim->startMomentum().e(), edep);
      if(debug_level_ > 1) std::cout << "EdepAna: primary Edep " << edep << " difference " << primsim->startMomentum().e() - edep << std::endl;

      if(front_trk_sp) {
        Hist->h_trk_front_energy_->Fill(front_trk_sp->momentum().mag(), Weight);
        Hist->h_trk_front_energy_diff_->Fill(front_trk_sp->momentum().mag() - primsim->startMomentum().e(), Weight);
        Hist->h_trk_front_energy_edep_diff_->Fill(edep - front_trk_sp->momentum().mag(), Weight);
        if(debug_level_ > 1) std::cout << "EdepAna: front tracker StepPointMC energy " << front_trk_sp->momentum().mag() << " difference " << primsim->startMomentum().e() - front_trk_sp->momentum().mag() << std::endl;
      }
    }

    if(shower_col_) {
      for(const auto& css : *shower_col_) {
        const double e = css.energyDepBirks();
        Hist->h_step_energy_->Fill(e, Weight);
        Hist->h_step_energy_vs_time_->Fill(css.time(), e, Weight);
      }
      Hist->h_total_calo_energy_->Fill(info_.calo_total_edep, Weight);
    }
    if(debug_level_ > 0) watch_->StopTime(__func__);
  }


  void EdepAna::analyze(const art::Event& event) {
    if(debug_level_ > 0) watch_->SetTime(__func__);

    //-------------------------------------------------------------
    // Retrieve data products for the event
    //-------------------------------------------------------------

    if(debug_level_ > 0) watch_->SetTime("data retrieval");

    // Primary particle
    art::Handle<PrimaryParticle> primaryH;
    event.getByLabel(primary_tag_, primaryH);
    primary_ = (primaryH.isValid()) ? primaryH.product() : nullptr;

    // Calo shower steps - sum energyDepBirks across the collection
    art::Handle<CaloShowerStepCollection> cssH;
    event.getByLabel(calo_shower_tag_, cssH);
    shower_col_ = (cssH.isValid()) ? cssH.product() : nullptr;

    // Straw gas steps
    art::Handle<StrawGasStepCollection> sgsH;
    event.getByLabel(straw_gas_tag_, sgsH);
    straw_step_col_ = (sgsH.isValid()) ? sgsH.product() : nullptr;

    // StepPointMC collection
    art::Handle<StepPointMCCollection> spH;
    event.getByLabel(step_point_tag_, spH);
    step_point_col_ = (spH.isValid()) ? spH.product() : nullptr;

    if(debug_level_ > 0) watch_->StopTime("data retrieval");

    //-------------------------------------------------------------
    // Compute event-level info
    //-------------------------------------------------------------

    if(debug_level_ > 0) watch_->SetTime("event info computation");

    double totalTrackerE = 0.;
    if(straw_step_col_) {
      for(const auto& step : *straw_step_col_) {
        const double e = step.ionizingEdep();
        totalTrackerE += e;
      }
    }

    double totalCaloE = 0.;
    if(shower_col_) {
      for(const auto& css : *shower_col_) {
        const double e = css.energyDepBirks();
        totalCaloE += e;
      }
    }

    // Take the first primary as the main one for default histogramming
    info_.primsim = (primary_ && !primary_->primarySimParticles().empty()) ? &(*primary_->primarySimParticles().front()) : nullptr;
    info_.front_trk_sp = simTrkFrontStep(info_.primsim);
    info_.primsim_edep = edepBySim(info_.primsim, true);
    info_.calo_total_edep = totalCaloE;

    info_.weight_ = 1.; // can be used to apply event weights if needed

    //-------------------------------------------------------------
    // Set ouput tree data
    //-------------------------------------------------------------

    data_.weight = info_.weight_;
    data_.run = event.run();
    data_.subrun = event.subRun();
    data_.event = event.event();
    data_.event_calo_edep = totalCaloE;
    data_.event_trk_edep  = totalTrackerE;
    data_.nprim = (primary_) ? int(primary_->primarySimParticles().size()) : 0;
    if(data_.nprim >= kMaxPrimaries) throw cet::exception("Analysis") << "Too many primary particles!"
                                                                           << "N(primaries) = " << data_.nprim
                                                                           << " > " << kMaxPrimaries;

    // Store information for each primary
    for(int iprim = 0; iprim < data_.nprim; ++iprim) {
      const auto* sim = &(*primary_->primarySimParticles().at(iprim));
      data_.prim_start_x  [iprim] = sim->startPosition().x();
      data_.prim_start_y  [iprim] = sim->startPosition().y();
      data_.prim_start_z  [iprim] = sim->startPosition().z();
      data_.prim_start_px [iprim] = sim->startMomentum().x();
      data_.prim_start_py [iprim] = sim->startMomentum().y();
      data_.prim_start_pz [iprim] = sim->startMomentum().z();
      data_.prim_start_e  [iprim] = sim->startMomentum().e();
      data_.prim_start_m  [iprim] = sim->startMomentum().m();
      data_.prim_start_pdg[iprim] = sim->pdgId();
      data_.prim_calo_edep[iprim] = edepBySim(sim, true);
      const auto front_trk_sp = simTrkFrontStep(sim);
      data_.prim_trk_front_energy[iprim] = (front_trk_sp) ? front_trk_sp->momentum().mag() : 0.f;
    }
    tree_->Fill();

    //-------------------------------------------------------------
    // Increment summary data
    //-------------------------------------------------------------

    total_events_ += info_.weight_;
    total_calo_edep_ += info_.weight_ * totalCaloE;
    if(totalCaloE > 50.) events_above_50_mev_ += info_.weight_;
    total_tracker_edep_ += info_.weight_ * totalTrackerE;

    if(debug_level_ > 0) watch_->StopTime("event info computation");

    //-------------------------------------------------------------
    // Fill histograms
    //-------------------------------------------------------------

    fillHistograms(hists_[0], info_.weight_);
    if(totalCaloE >  1.) fillHistograms(hists_[1], info_.weight_);
    if(totalCaloE > 10.) fillHistograms(hists_[2], info_.weight_);
    if(totalCaloE > 50.) fillHistograms(hists_[3], info_.weight_);
    if(debug_level_ > 0) watch_->StopTime(__func__);
  }


  void EdepAna::endJob() {
    if(debug_level_ > 0) std::cout << "[EdepAna::" << __func__ << "]\n" << *watch_ << std::endl;

    //-------------------------------------------------------------
    // Fit the total gen -> calo edep response
    //-------------------------------------------------------------

    TH1* h = (hists_[2]) ? hists_[2]->h_primary_energy_edep_diff_ : nullptr;
    if(h && h->GetEntries() > 100) {
      double mean_seed, fwhm_seed;
      get_landau_seed(h, mean_seed, fwhm_seed);
      double mean, fwhm;
      const int fit_status = fit_landau(h, mean, fwhm);

      std::cout << std::format("Primary energy - Edep fit: status = {}, mean = {:.2f} MeV, FWHM = {:.2f} MeV",
                               fit_status, mean, fwhm) << std::endl;
      std::cout << std::format("Primary energy - Edep distribution: mean = {:.2f} MeV, RMS = {:.2f} MeV, MPV = {:.2f} MeV, FWHM = {:.2f} MeV",
                               h->GetMean(), h->GetRMS(), mean_seed, fwhm_seed) << std::endl;
    }

    //-------------------------------------------------------------
    // Fit the energy loss from gen -> tracker front
    //-------------------------------------------------------------

    h = (hists_[2]) ? hists_[2]->h_trk_front_energy_diff_ : nullptr;
    if(h) {
      TH1* h_ref = (hists_[0]) ? hists_[0]->h_total_calo_energy_ : nullptr;
      const double eff = h->GetEntries() > 0 && h_ref ? h->GetEntries() * 1./ h_ref->GetEntries() : 0.;
      double mpv_seed, fwhm_seed;
      get_landau_seed(h, mpv_seed, fwhm_seed);
      double mpv, fwhm;
      const int fit_status = fit_landau(h, mpv, fwhm);
      std::cout << std::format("Tracker front - primary energy fit: status = {}, mean = {:.2f} MeV, FWHM = {:.2f} MeV",
                               fit_status, mpv, fwhm) << std::endl;
      std::cout << std::format("Tracker front StepPointMC energy - primary energy distribution: MPV = {:.2f} MeV, FWHM = {:.2f} MeV, efficiency = {:.4g}", mpv_seed, fwhm_seed, eff) << std::endl;
    }

    //-------------------------------------------------------------
    // Fit the tracker front -> calo edep response
    //-------------------------------------------------------------

    h = (hists_[2]) ? hists_[2]->h_trk_front_energy_edep_diff_ : nullptr;
    if(h) {
      TH1* h_ref = (hists_[0]) ? hists_[0]->h_total_calo_energy_ : nullptr;
      const double eff = h->GetEntries() > 0 && h_ref ? h->GetEntries() * 1./ h_ref->GetEntries() : 0.;
      double mpv, fwhm;
      get_landau_seed(h, mpv, fwhm);
      std::cout << std::format("Primary Edep - tracker front StepPointMC energy distribution: MPV = {:.2f} MeV, FWHM = {:.2f} MeV, efficiency = {:.4g}", mpv, fwhm, eff) << std::endl;
    }

    //-------------------------------------------------------------
    // Report summary data
    //-------------------------------------------------------------

    const double averageCaloEdep    = (total_events_ > 0) ? total_calo_edep_    / total_events_ : 0.;
    const double averageCaloEdepGen = (ngen_         > 0) ? total_calo_edep_    / ngen_         : 0.;
    const double averageTrkEdep     = (total_events_ > 0) ? total_tracker_edep_ / total_events_ : 0.;
    const double averageTrkEdepGen  = (ngen_         > 0) ? total_tracker_edep_ / ngen_         : 0.;
    const double eventRate          = (ngen_         > 0) ? total_events_       / ngen_         : 0.;

    std::cout
      << "EdepAna summary:\n"
      << "  Saw " << total_events_ << " events" << " (" << ngen_ << " gen events) --> output rate = "
      << eventRate << " events / gen event" << std::endl
      << "  Average calo energy deposition per event: " << averageCaloEdep << " MeV" << std::endl
      << "  Average calo energy deposition per gen event: " << averageCaloEdepGen << " MeV" << std::endl
      << "  Events with calo Edep > 50 MeV: " << events_above_50_mev_ << std::endl
      << "  Average tracker energy deposition per event: " << averageTrkEdep << " MeV" << std::endl
      << "  Average tracker energy deposition per gen event: " << averageTrkEdepGen << " MeV" << std::endl;
  }

}

using mu2e::EdepAna;
DEFINE_ART_MODULE(EdepAna)
