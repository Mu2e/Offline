// EdepAna: simple analyzer for primary particle energy losses and calo/tracker energy deposition
// Original author: Michael MacKenzie, 2026

// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"

// Offline
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

// CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TF1.h"
#include "TFitResult.h"

// c++
#include <cmath>
#include <format>
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <optional>

namespace mu2e {

  class EdepAna : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> primaryTag     { Name("Primary")                 , Comment("PrimaryParticle tag") };
      fhicl::Atom<art::InputTag> caloShowerTag  { Name("CaloShowerStepCollection"), Comment("CaloShowerStep collection") };
      fhicl::Atom<art::InputTag> strawGasStepTag{ Name("StrawGasStepCollection")  , Comment("StrawGasStep collection") };
      fhicl::Atom<art::InputTag> stepPointTag   { Name("StepPointMCCollection")   , Comment("StepPointMC collection") };
      fhicl::OptionalAtom<art::InputTag> crvStepTag{ Name("CrvStepCollection")    , Comment("CrvStep collection; if set, CRV occupancy and dose histograms are filled") };
      fhicl::Atom<art::InputTag> genCountTag    { Name("GenCountTag")             , Comment("GenCount tag"), "genCounter" };
      // The Step time windows cut the raw step time in the tree and the end-of-job summary; the Hit time windows cut
      // the folded time (event-window frame, see FoldTimeOffset) in the occupancy hit-like histograms
      fhicl::Atom<double>        trkStepEMin    { Name("TrkStepEnergyThreshold")  , Comment("StrawGasStep ionizing energy threshold for counting hits (MeV)") };
      fhicl::Atom<double>        trkStepTMin    { Name("TrkStepTimeMin")          , Comment("StrawGasStep minimum time for counting hits (ns)") };
      fhicl::Atom<double>        trkStepTMax    { Name("TrkStepTimeMax")          , Comment("StrawGasStep maximum time for counting hits (ns)") };
      fhicl::Atom<double>        caloStepEMin   { Name("CaloStepEnergyThreshold") , Comment("CaloShowerStep visible energy threshold for counting hits (MeV)") };
      fhicl::Atom<double>        caloStepTMin   { Name("CaloStepTimeMin")         , Comment("CaloShowerStep minimum time for counting hits (ns)") };
      fhicl::Atom<double>        caloStepTMax   { Name("CaloStepTimeMax")         , Comment("CaloShowerStep maximum time for counting hits (ns)") };
      fhicl::Atom<double>        foldTimeOffset { Name("FoldTimeOffset")          , Comment("Offset added to step times before folding into the microbunch (ns): the nominal ProtonBunchTimeMC, "
                                                                                              "-TimeFromProtonsToDRMarker, so folded times are in the digitizers' event-window frame. "
                                                                                              "The per-event phase (within one 25 ns clock period) is not included. 0 folds the raw step time") };
      fhicl::Atom<double>        trkHitTMin     { Name("TrkHitTimeMin")           , Comment("Minimum folded time for the tracker hit-like histograms (ns, event-window frame)") };
      fhicl::Atom<double>        trkHitTMax     { Name("TrkHitTimeMax")           , Comment("Maximum folded time for the tracker hit-like histograms (ns, event-window frame)") };
      fhicl::Atom<double>        caloHitTMin    { Name("CaloHitTimeMin")          , Comment("Minimum folded time for the calo hit-like histograms (ns, event-window frame)") };
      fhicl::Atom<double>        caloHitTMax    { Name("CaloHitTimeMax")          , Comment("Maximum folded time for the calo hit-like histograms (ns, event-window frame)") };
      fhicl::Atom<int>           debugLevel     { Name("DebugLevel")              , Comment("Debug level"), 0 };
    };

    // Histograms
    struct Hist_t {
      TH1* h_primary_energy_;
      TH1* h_primary_pdg_;
      TH1* h_primary_start_z_;
      TH1* h_primary_start_r_;
      TH1* h_total_calo_energy_vis_;
      TH1* h_step_energy_; // individual CaloShowerStep energies
      TH2* h_step_energy_vs_time_; // step energy vs time
      TH1* h_trk_front_p_;
      TH1* h_trk_front_energy_;
      TH1* h_trk_front_energy_diff_;
      TH2* h_primary_vs_edep_;
      TH1* h_primary_edep_;
      TH1* h_primary_energy_edep_diff_;
      TH1* h_trk_front_energy_edep_diff_;
    };
    enum {kMaxHists = 100};

    // Occupancy (step count) and dose (deposited energy) histograms, filled once per event for all steps.
    // Each observable is a pair of histograms: [kOcc] filled with the step weight, [kEdep] filled with weight * edep.
    // Tracker: StrawGasStep::ionizingEdep, Calo: CaloShowerStep::energyDepG4 (true deposited energy, not Birks-quenched
    // visible energy), CRV: CrvStep::visibleEDep (the only energy stored in a CrvStep).
    // Positions are in the detector frame (origin at the tracker center). Raw-time axes span 10 us (later steps,
    // e.g. from muon stops, go to overflow); the folded-time histograms fold the step time plus FoldTimeOffset into the
    // microbunch, as the tracker digitizer does, giving times in the event-window frame.
    // calo_disk<d>_ncrystals_r is geometry, filled once per job: hadd-ed files hold one copy per job.
    enum {kOcc = 0, kEdep = 1};
    // Crystal-centre radius bins (mm): the crystals are 34.3 mm wide, so bins of this width hold tens of them (disk crystal centres span ~390-650 mm, so no bin is empty)
    static constexpr int    kNRCrysBins = 8;
    static constexpr double kRCrysMin   = 385.;
    static constexpr double kRCrysMax   = 665.;
    using HistPair_t = std::array<TH1*, 2>;
    struct OccHist_t {
      TH1* h_nevents_ = nullptr; // number of (weighted) events analyzed
      TH1* h_ngen_    = nullptr; // number of generated events, set in endJob
      // Tracker
      HistPair_t trk_plane_, trk_panel_, trk_straw_, trk_r_, trk_z_, trk_time_, trk_time_folded_;
      std::array<TH2*, 2> trk_rz_ = {nullptr, nullptr};
      // Calo
      HistPair_t calo_disk_, calo_crystal_, calo_r_, calo_z_, calo_time_, calo_time_folded_;
      std::array<HistPair_t, 2> calo_disk_r_; // r by disk
      // phi-averaged by disk: binned in the crystal-centre radius (no sub-crystal structure), with the crystal count per bin
      std::array<HistPair_t, 2> calo_disk_rcrys_;
      std::array<TH1*, 2> calo_disk_rcrys_hits_ = {nullptr, nullptr};
      std::array<TH1*, 2> calo_disk_ncrystals_r_ = {nullptr, nullptr};
      std::array<std::array<TH2*, 2>, 2> calo_disk_xy_; // [disk][kOcc/kEdep]
      std::array<TH2*, 2> calo_rz_ = {nullptr, nullptr};
      // Hit-like counts (steps above the energy threshold), single histograms; _hits only for folded time in the Hit time window
      TH1* trk_time_folded_ecut_ = nullptr;
      TH1* calo_time_folded_ecut_ = nullptr;
      TH1* trk_straw_hits_ = nullptr;
      TH1* calo_crystal_hits_ = nullptr; // steps above CaloStepEnergyThreshold (Birks), folded time in [CaloHitTimeMin, CaloHitTimeMax]
      // CRV (booked only if a CrvStep collection is configured)
      HistPair_t crv_sector_, crv_bar_, crv_z_, crv_time_, crv_time_folded_;
    };

    // Event info
    struct Info_t {
      const SimParticle* primsim = nullptr;
      const StepPointMC* front_trk_sp = nullptr;
      double primsim_edep_vis = 0.;
      double calo_total_edep = 0.;
      double calo_total_edep_vis = 0.;
      double weight_ = 1.;
    };

    // Output tree information
    enum {kMaxPrimaries = 10};
    struct Tree_t {
      int   nprim;
      float prim_start_x           [kMaxPrimaries];
      float prim_start_y           [kMaxPrimaries];
      float prim_start_z           [kMaxPrimaries];
      float prim_start_px          [kMaxPrimaries];
      float prim_start_py          [kMaxPrimaries];
      float prim_start_pz          [kMaxPrimaries];
      float prim_start_e           [kMaxPrimaries];
      float prim_start_m           [kMaxPrimaries];
      int   prim_start_pdg         [kMaxPrimaries];
      float prim_calo_edep         [kMaxPrimaries];
      float prim_calo_edep_vis     [kMaxPrimaries];
      float prim_trk_front_p       [kMaxPrimaries];
      float prim_trk_front_energy  [kMaxPrimaries];
      float prim_trk_front_r_center[kMaxPrimaries];
      float prim_trk_front_radius  [kMaxPrimaries];
      float prim_trk_front_time    [kMaxPrimaries];
      float prim_trk_front_cos     [kMaxPrimaries];
      float event_calo_edep;
      float event_calo_edep_vis;
      float event_trk_edep;
      int   event_calo_nsteps;           // all CaloShowerSteps
      int   event_calo_nsteps_ecut;      // above the energy threshold
      int   event_calo_nsteps_ecut_tcut; // above the energy threshold and in the time window
      int   event_trk_nsteps;            // all StrawGasSteps
      int   event_trk_nsteps_ecut;
      int   event_trk_nsteps_ecut_tcut;
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
    virtual void beginRun(const art::Run& run) override;
    virtual void beginSubRun(const art::SubRun& subrun) override;

    void bookHistograms(const int index, const char* title);
    void fillHistograms(Hist_t* Hist, const double Weight = 1.);
    void bookOccupancyHistograms();
    void fillOccupancyHistograms(const double Weight = 1.);

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

    double edepBySim(const SimParticle* sim, const bool check_parents, const bool do_vis = false) const {
      if(!sim) return 0.;
      if(!shower_col_) return 0.;
      double edep = 0.;
      for(const auto& shower : *shower_col_) {
        if(shower.simParticle()->id() == sim->id() || (check_parents && isDescendant(sim, &*shower.simParticle()))) {
          edep += (do_vis) ? shower.energyDepBirks() : shower.energyDepG4();
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
        }
      }
      return front_trk_sp;
    }

    float energyFromStepPoint(const SimParticle* sim, const StepPointMC* step) {
      if(!sim || !step) return 0.f;
      // use the sim mass at start to get the energy from the step point
      return std::sqrt(std::pow(step->momentum().mag(), 2) + std::pow(sim->startMomentum().m(), 2));
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
    TF1* f = new TF1("landau_func", landau_func, mean_seed - fwhm_seed, std::min(0., mean_seed + fwhm_seed), 4);
    f->SetParNames("Norm", "#mu", "a", "b");
    f->SetParLimits(1, -50., 0.); // mean
    f->FixParameter(2, 0.5); // a
    f->SetParLimits(3, 0.01, 100.); // b
    f->SetParameters(h->GetMaximum(), mean_seed, 0.5, fwhm_seed/5.);
    auto fit_res = h->Fit(f, "QNSR");
    if(!fit_res.Get()) { // e.g. empty fit range
      std::cerr << "Fit failed: no fit result" << std::endl;
      delete f;
      return -1;
    }
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
    art::InputTag gen_count_tag_;
    std::optional<art::InputTag> crv_step_tag_;
    double trk_step_emin_;
    double trk_step_tmin_;
    double trk_step_tmax_;
    double calo_step_emin_;
    double calo_step_tmin_;
    double calo_step_tmax_;
    double fold_time_offset_;
    double trk_hit_tmin_;
    double trk_hit_tmax_;
    double calo_hit_tmin_;
    double calo_hit_tmax_;
    int debug_level_;
    std::unique_ptr<StopWatch> watch_ = std::make_unique<StopWatch>();

    Hist_t* hists_[kMaxHists];
    std::unique_ptr<OccHist_t> occ_;
    // Fold a step time into the event-window frame [0, mbtime_) like StrawDigisFromStrawGasSteps::microbunchTime:
    // add the proton bunch time, then take it modulo the microbunch period (negative times wrap forward)
    double foldTime(const double t) const { const double f = std::fmod(t + fold_time_offset_, mbtime_); return f < 0. ? f + mbtime_ : f; }
    double mbtime_ = 0.; // microbunch period (ns), from PhysicsParams, set when booking
    Info_t info_;
    const PrimaryParticle*                 primary_        = nullptr;
    const CaloShowerStepCollection*        shower_col_     = nullptr;
    const StrawGasStepCollection*          straw_step_col_ = nullptr;
    const StepPointMCCollection*           step_point_col_ = nullptr;
    const CrvStepCollection*               crv_col_        = nullptr;
    double total_events_ = 0.;
    double events_above_50_mev_ = 0.;
    double total_calo_edep_ = 0.;
    double total_calo_edep_vis_ = 0.;
    double total_tracker_edep_ = 0.;
    double total_calo_nsteps_[3] = {0., 0., 0.}; // raw, energy cut, energy + time cut
    double total_trk_nsteps_ [3] = {0., 0., 0.};
    unsigned long ngen_ = 0;
  };


  EdepAna::EdepAna(const Parameters& conf)
    : art::EDAnalyzer{conf}
    , primary_tag_(conf().primaryTag())
    , calo_shower_tag_(conf().caloShowerTag())
    , straw_gas_tag_(conf().strawGasStepTag())
    , step_point_tag_(conf().stepPointTag())
    , gen_count_tag_(conf().genCountTag())
    , crv_step_tag_(std::nullopt)
    , trk_step_emin_(conf().trkStepEMin())
    , trk_step_tmin_(conf().trkStepTMin())
    , trk_step_tmax_(conf().trkStepTMax())
    , calo_step_emin_(conf().caloStepEMin())
    , calo_step_tmin_(conf().caloStepTMin())
    , calo_step_tmax_(conf().caloStepTMax())
    , fold_time_offset_(conf().foldTimeOffset())
    , trk_hit_tmin_(conf().trkHitTMin())
    , trk_hit_tmax_(conf().trkHitTMax())
    , calo_hit_tmin_(conf().caloHitTMin())
    , calo_hit_tmax_(conf().caloHitTMax())
    , debug_level_(conf().debugLevel())
  {
    art::InputTag crv_tag;
    if(conf().crvStepTag(crv_tag)) crv_step_tag_ = crv_tag;

    // register products consumed
    consumes<PrimaryParticle>(primary_tag_);
    consumes<CaloShowerStepCollection>(calo_shower_tag_);
    consumes<StrawGasStepCollection>(straw_gas_tag_);
    consumes<StepPointMCCollection>(step_point_tag_);
    if(crv_step_tag_) consumes<CrvStepCollection>(*crv_step_tag_);
    consumes<GenEventCount, art::InSubRun>(gen_count_tag_);

    for(int i = 0; i < kMaxHists; ++i) hists_[i] = nullptr;
    bookHistograms(0, "all events");
    bookHistograms(1, "edep 1 MeV");
    bookHistograms(2, "edep 10 MeV");
    bookHistograms(3, "edep 50 MeV");

    // Book the output TTree
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("tree", "Energy deposition data");
    tree_->Branch("nprimaries"                , &data_.nprim);
    tree_->Branch("primary_start_x"           , data_.prim_start_x           , "primary_start_x[nprimaries]/F");
    tree_->Branch("primary_start_y"           , data_.prim_start_y           , "primary_start_y[nprimaries]/F");
    tree_->Branch("primary_start_z"           , data_.prim_start_z           , "primary_start_z[nprimaries]/F");
    tree_->Branch("primary_start_px"          , data_.prim_start_px          , "primary_start_px[nprimaries]/F");
    tree_->Branch("primary_start_py"          , data_.prim_start_py          , "primary_start_py[nprimaries]/F");
    tree_->Branch("primary_start_pz"          , data_.prim_start_pz          , "primary_start_pz[nprimaries]/F");
    tree_->Branch("primary_start_e"           , data_.prim_start_e           , "primary_start_e[nprimaries]/F");
    tree_->Branch("primary_start_m"           , data_.prim_start_m           , "primary_start_m[nprimaries]/F");
    tree_->Branch("primary_start_pdg"         , data_.prim_start_pdg         , "primary_start_pdg[nprimaries]/I");
    tree_->Branch("primary_calo_edep"         , data_.prim_calo_edep         , "primary_calo_edep[nprimaries]/F");
    tree_->Branch("primary_calo_edep_vis"     , data_.prim_calo_edep_vis     , "primary_calo_edep_vis[nprimaries]/F");
    tree_->Branch("primary_trk_front_p"       , data_.prim_trk_front_p       , "primary_trk_front_p[nprimaries]/F");
    tree_->Branch("primary_trk_front_energy"  , data_.prim_trk_front_energy  , "primary_trk_front_energy[nprimaries]/F");
    tree_->Branch("primary_trk_front_r_center", data_.prim_trk_front_r_center, "primary_trk_front_r_center[nprimaries]/F");
    tree_->Branch("primary_trk_front_radius"  , data_.prim_trk_front_radius  , "primary_trk_front_radius[nprimaries]/F");
    tree_->Branch("primary_trk_front_time"    , data_.prim_trk_front_time    , "primary_trk_front_time[nprimaries]/F");
    tree_->Branch("primary_trk_front_cos"     , data_.prim_trk_front_cos     , "primary_trk_front_cos[nprimaries]/F");
    tree_->Branch("event_calo_edep"         , &data_.event_calo_edep);
    tree_->Branch("event_calo_edep_vis"     , &data_.event_calo_edep_vis);
    tree_->Branch("event_trk_edep"          , &data_.event_trk_edep);
    tree_->Branch("event_calo_nsteps"          , &data_.event_calo_nsteps);
    tree_->Branch("event_calo_nsteps_ecut"     , &data_.event_calo_nsteps_ecut);
    tree_->Branch("event_calo_nsteps_ecut_tcut", &data_.event_calo_nsteps_ecut_tcut);
    tree_->Branch("event_trk_nsteps"           , &data_.event_trk_nsteps);
    tree_->Branch("event_trk_nsteps_ecut"      , &data_.event_trk_nsteps_ecut);
    tree_->Branch("event_trk_nsteps_ecut_tcut" , &data_.event_trk_nsteps_ecut_tcut);
    tree_->Branch("run"                     , &data_.run);
    tree_->Branch("subrun"                  , &data_.subrun);
    tree_->Branch("event"                   , &data_.event);
    tree_->Branch("weight"                  , &data_.weight);
    tree_->Branch("ngen"                    , &data_.ngen);

    if(debug_level_ > 0) watch_->Calibrate();
  }

  //--------------------------------------------------------------------------------------
  void EdepAna::beginRun(const art::Run&) {
    // The occupancy histograms need the geometry (crystal, bar counts), which is only available once a run starts
    if(!occ_) bookOccupancyHistograms();
  }

  //--------------------------------------------------------------------------------------
  void EdepAna::beginSubRun(const art::SubRun& subrun) {
    // Get the generator counter
    auto genCounterHandle = subrun.getHandle<GenEventCount>(gen_count_tag_);
    if(genCounterHandle.isValid()) {
      ngen_ += genCounterHandle->count();
    } else {
      std::cerr << "Warning: GenEventCount not found in subrun" << std::endl;
    }
    data_.ngen = ngen_; // update the output tree data
  }

  void EdepAna::bookHistograms(const int index, const char* title) {
    if(index >= kMaxHists) throw cet::exception("Analysis") << "Too many histograms!";
    hists_[index] = new Hist_t;
    auto Hist = hists_[index];

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("hist_{}", index), title);

   Hist->h_primary_energy_     = dir.make<TH1F>("primary_energy"    , "Primary energy;Energy (MeV)"              , 150,    0.,  150.);
   Hist->h_primary_pdg_        = dir.make<TH1D>("primary_pdg"       , "Primary PDG ID;PDG ID"                    ,  50,  -25.,   25.);
   Hist->h_primary_start_z_    = dir.make<TH1D>("primary_start_z"   , "Primary start z;Start z (mm)"             , 500, 3000., 8000.);
   Hist->h_primary_start_r_    = dir.make<TH1D>("primary_start_r"   , "Primary start radius;Start radius (mm)"   , 100,    0.,  200.);
   Hist->h_total_calo_energy_vis_  = dir.make<TH1F>("total_calo_energy" , "Total visible Calo energy from steps;Energy (MeV)", 200,    0.,  200.);
   Hist->h_step_energy_        = dir.make<TH1F>("calo_step_energy"  , "Visible CaloShowerStep energy;E_{step} (MeV)"     , 100,    0.,  100.);
   Hist->h_trk_front_p_        = dir.make<TH1F>("trk_front_p"       , "Momentum of StepPointMC at front of tracker;p (MeV/c)", 1500, 0., 150.);
   Hist->h_trk_front_energy_   = dir.make<TH1F>("trk_front_energy"  , "Energy of StepPointMC at front of tracker;Energy (MeV)", 1500, 0., 150.);
   Hist->h_trk_front_energy_diff_ = dir.make<TH1F>("trk_front_energy_diff", "Energy difference of StepPointMC at front of tracker and primary;Energy (MeV)", 500, -100., 0.);
   Hist->h_primary_energy_edep_diff_ = dir.make<TH1F>("primary_energy_edep_diff", "Primary Edep - energy;Energy (MeV)", 300, -150., 0.);
   Hist->h_primary_edep_ = dir.make<TH1F>("primary_edep", "Primary Edep;Energy (MeV)", 300, 0., 150.);
   Hist->h_primary_vs_edep_ = dir.make<TH2F>("primary_vs_edep", "Primary energy vs Edep;Primary energy (MeV);Edep (MeV)", 150, 0., 150., 150, 0., 150.);
   Hist->h_step_energy_vs_time_ = dir.make<TH2F>("calo_step_energy_vs_time", "Visible CaloShowerStep energy vs time;Time (ns);Energy (MeV)", 40, 0., 2000., 100, 0., 100.);
   Hist->h_trk_front_energy_edep_diff_ = dir.make<TH1F>("trk_front_energy_edep_diff", "Energy difference of StepPointMC at front of tracker and primary Edep;Energy (MeV)", 500, -50., 0.);
  }

  void EdepAna::fillHistograms(Hist_t* Hist, const double Weight) {
    if(debug_level_ > 0) watch_->SetTime(__func__);
    if(!Hist) throw cet::exception("Analysis") << "Uninitialized histogram book!";

    const SimParticle* primsim = info_.primsim;
    const StepPointMC* front_trk_sp = info_.front_trk_sp;
    const double primsim_edep_vis = info_.primsim_edep_vis;


    if(primsim) {
      const auto startMomentum = primsim->startMomentum();
      Hist->h_primary_energy_->Fill(startMomentum.e(), Weight);
      Hist->h_primary_pdg_->Fill(primsim->pdgId(), Weight);
      if(debug_level_ > 1) std::cout << "EdepAna: primary energy " << startMomentum.e() << " PDG " << primsim->pdgId() << std::endl;
      Hist->h_primary_start_z_->Fill(primsim->startPosition().z(), Weight);
      GeomHandle<DetectorSystem> det;
      Hist->h_primary_start_r_->Fill(det->toDetector(primsim->startPosition()).perp(), Weight);
      const double edep = primsim_edep_vis;
      Hist->h_primary_edep_->Fill(edep, Weight);
      Hist->h_primary_energy_edep_diff_->Fill(edep - startMomentum.e(), Weight);
      Hist->h_primary_vs_edep_->Fill(startMomentum.e(), edep);
      if(debug_level_ > 1) std::cout << "EdepAna: primary Edep " << edep << " difference " << startMomentum.e() - edep << std::endl;

      if(front_trk_sp) {
        Hist->h_trk_front_p_->Fill(front_trk_sp->momentum().mag(), Weight);
        const float step_energy = energyFromStepPoint(primsim, front_trk_sp);
        Hist->h_trk_front_energy_->Fill(step_energy, Weight);
        Hist->h_trk_front_energy_diff_->Fill(step_energy - startMomentum.e(), Weight);
        Hist->h_trk_front_energy_edep_diff_->Fill(edep - step_energy, Weight);
        if(debug_level_ > 1) std::cout << "EdepAna: front tracker StepPointMC energy " << step_energy << " difference " << startMomentum.e() - step_energy << std::endl;
      }
    }

    if(shower_col_) {
      for(const auto& css : *shower_col_) {
        const double e = css.energyDepBirks();
        Hist->h_step_energy_->Fill(e, Weight);
        Hist->h_step_energy_vs_time_->Fill(css.time(), e, Weight);
      }
      Hist->h_total_calo_energy_vis_->Fill(info_.calo_total_edep_vis, Weight);
    }
    if(debug_level_ > 0) watch_->StopTime(__func__);
  }


  //--------------------------------------------------------------------------------------
  void EdepAna::bookOccupancyHistograms() {
    occ_ = std::make_unique<OccHist_t>();
    // Microbunch period, as used by the digitizers to fold step times
    mbtime_ = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("occupancy", "Step occupancy and energy deposition (dose) histograms");
    GeomHandle<Calorimeter> calo;

    const char* kind [2] = {"occ", "edep"};
    const char* ylab [2] = {"Steps", "Energy deposited (MeV)"};
    const char* ktitle[2] = {"occupancy", "energy deposition"};

    // Book a pair of 1D histograms (occupancy, edep) named <name>_occ / <name>_edep
    auto book1 = [&](HistPair_t& hp, const char* name, const char* title, const char* xlab, int nb, double lo, double hi) {
      for(int k = 0; k < 2; ++k)
      { hp[k] = dir.make<TH1D>(std::format("{}_{}", name, kind[k]).c_str(),
                               std::format("{} {};{};{}", title, ktitle[k], xlab, ylab[k]).c_str(), nb, lo, hi);
        hp[k]->Sumw2(); }
    };
    auto book2 = [&](std::array<TH2*, 2>& hp, const char* name, const char* title, const char* xlab, int nx, double xlo, double xhi,
                     const char* ylabel, int ny, double ylo, double yhi) {
      for(int k = 0; k < 2; ++k)
      { hp[k] = dir.make<TH2D>(std::format("{}_{}", name, kind[k]).c_str(),
                               std::format("{} {} ({});{};{}", title, ktitle[k], ylab[k], xlab, ylabel).c_str(), nx, xlo, xhi, ny, ylo, yhi);
        hp[k]->Sumw2(); }
    };

    occ_->h_nevents_ = dir.make<TH1D>("nevents", "Events analyzed (weighted);;N(events)", 1, 0., 1.);
    occ_->h_ngen_    = dir.make<TH1D>("ngen"   , "Generated events;;N(gen)"             , 1, 0., 1.);

    // Tracker, detector frame: r ~ 350-700 mm, |z| < ~1600 mm
    book1(occ_->trk_plane_ , "trk_plane" , "Tracker by plane" , "Plane"       , StrawId::_nplanes , 0., StrawId::_nplanes);
    book1(occ_->trk_panel_ , "trk_panel" , "Tracker by panel" , "Unique panel", StrawId::_nupanels, 0., StrawId::_nupanels);
    book1(occ_->trk_straw_ , "trk_straw" , "Tracker by straw" , "Unique straw", StrawId::_nustraws, 0., StrawId::_nustraws);
    book1(occ_->trk_r_     , "trk_r"     , "Tracker by radius", "r (mm)"      , 175, 300., 750.);
    book1(occ_->trk_z_     , "trk_z"     , "Tracker by z"     , "z (mm)"      , 320, -1600., 1600.);
    book1(occ_->trk_time_  , "trk_time"  , "Tracker by time"  , "Time (ns)"   , 500, 0., 10000.);
    book1(occ_->trk_time_folded_, "trk_time_folded", "Tracker by time folded into the microbunch", "Time mod microbunch period (ns)", 100, 0., mbtime_);
    occ_->trk_time_folded_ecut_ = dir.make<TH1D>("trk_time_folded_ecut", "Tracker steps above the energy threshold by folded time;Time mod microbunch period (ns);Steps", 100, 0., mbtime_);
    occ_->trk_straw_hits_ = dir.make<TH1D>("trk_straw_hits", "Tracker steps above the energy threshold per straw, folded time in the TrkHitTime window;Unique straw;Steps", StrawId::_nustraws, 0., StrawId::_nustraws);
    for(TH1* h : {occ_->trk_time_folded_ecut_, occ_->trk_straw_hits_}) h->Sumw2();
    book2(occ_->trk_rz_    , "trk_rz"    , "Tracker r vs z"   , "z (mm)", 320, -1600., 1600., "r (mm)", 175, 300., 750.);

    // Calorimeter, detector frame: disk radii ~370-660 mm
    const int ncrys = calo->nCrystals();
    book1(occ_->calo_disk_   , "calo_disk"   , "Calo by disk"   , "Disk"        , 2, 0., 2.);
    book1(occ_->calo_crystal_, "calo_crystal", "Calo by crystal", "Crystal ID"  , ncrys, 0., ncrys);
    book1(occ_->calo_r_      , "calo_r"      , "Calo by radius" , "r (mm)"      , 100, 300., 700.);
    book1(occ_->calo_z_      , "calo_z"      , "Calo by z"      , "z (mm)"      , 250, 1500., 3500.);
    book1(occ_->calo_time_   , "calo_time"   , "Calo by time"   , "Time (ns)"   , 500, 0., 10000.);
    book1(occ_->calo_time_folded_, "calo_time_folded", "Calo by time folded into the microbunch", "Time mod microbunch period (ns)", 100, 0., mbtime_);
    occ_->calo_time_folded_ecut_ = dir.make<TH1D>("calo_time_folded_ecut", "Calo steps above the energy threshold by folded time;Time mod microbunch period (ns);Steps", 100, 0., mbtime_);
    occ_->calo_time_folded_ecut_->Sumw2();
    occ_->calo_crystal_hits_ = dir.make<TH1D>("calo_crystal_hits", "Calo steps above the energy threshold per crystal, folded time in the CaloHitTime window;Crystal ID;Steps", ncrys, 0., ncrys);
    occ_->calo_crystal_hits_->Sumw2();
    book2(occ_->calo_rz_     , "calo_rz"     , "Calo r vs z"    , "z (mm)", 250, 1500., 3500., "r (mm)", 100, 300., 700.);
    for(int d = 0; d < 2; ++d) {
      book1(occ_->calo_disk_r_[d], std::format("calo_disk{}_r", d).c_str(), std::format("Calo disk {} by radius", d).c_str(), "r (mm)", 100, 300., 700.);
      book2(occ_->calo_disk_xy_[d], std::format("calo_disk{}_xy", d).c_str(), std::format("Calo disk {} x vs y", d).c_str(),
            "x (mm)", 140, -700., 700., "y (mm)", 140, -700., 700.);
      // Phi average: binned in the radius of the crystal centre, so a bin's step count divided by its crystal count
      // (calo_disk<d>_ncrystals_r, filled here once per crystal) is the mean per crystal at that radius
      const std::string dn = std::format("calo_disk{}_rcrys", d);
      book1(occ_->calo_disk_rcrys_[d], dn.c_str(), std::format("Calo disk {} by crystal-centre radius", d).c_str(), "Crystal-centre r (mm)", kNRCrysBins, kRCrysMin, kRCrysMax);
      occ_->calo_disk_rcrys_hits_[d] = dir.make<TH1D>((dn + "_hits").c_str(), std::format("Calo disk {} steps above the energy threshold by crystal-centre radius, folded time in the CaloHitTime window;Crystal-centre r (mm);Steps", d).c_str(), kNRCrysBins, kRCrysMin, kRCrysMax);
      occ_->calo_disk_rcrys_hits_[d]->Sumw2();
      occ_->calo_disk_ncrystals_r_[d] = dir.make<TH1D>(std::format("calo_disk{}_ncrystals_r", d).c_str(), std::format("Calo disk {} crystals by crystal-centre radius;Crystal-centre r (mm);Crystals", d).c_str(), kNRCrysBins, kRCrysMin, kRCrysMax);
    }
    // Geometry: every crystal's centre radius in the detector frame
    {
      GeomHandle<DetectorSystem> det;
      for(int cid = 0; cid < ncrys; ++cid) {
        const auto& cry = calo->crystal(cid);
        occ_->calo_disk_ncrystals_r_[cry.diskID()]->Fill(det->toDetector(cry.position()).perp());
      }
    }

    // CRV, optional; planar modules so only sector, bar, z and time
    if(crv_step_tag_) {
      GeomHandle<CosmicRayShield> crs;
      const int nsectors = crs->getCRSScintillatorShields().size();
      const int nbars    = crs->getAllCRSScintillatorBars().size();
      book1(occ_->crv_sector_, "crv_sector", "CRV by sector"   , "Sector (shield number)", nsectors, 0., nsectors);
      book1(occ_->crv_bar_   , "crv_bar"   , "CRV by bar"      , "Bar index"             , nbars   , 0., nbars);
      book1(occ_->crv_z_     , "crv_z"     , "CRV by z"        , "z (mm)"                , 300, -15000., 15000.);
      book1(occ_->crv_time_  , "crv_time"  , "CRV by time"     , "Time (ns)"             , 500, 0., 10000.);
      book1(occ_->crv_time_folded_, "crv_time_folded", "CRV by time folded into the microbunch", "Time mod microbunch period (ns)", 100, 0., mbtime_);
    }
  }

  //--------------------------------------------------------------------------------------
  // One pass over each collection; all steps, no energy or time cuts
  void EdepAna::fillOccupancyHistograms(const double Weight) {
    if(debug_level_ > 0) watch_->SetTime(__func__);
    GeomHandle<DetectorSystem> det;
    occ_->h_nevents_->Fill(0.5, Weight);

    if(straw_step_col_) {
      for(const auto& step : *straw_step_col_) {
        // StrawGasStep positions are already in the detector (tracker) frame: they are compared directly to the
        // Tracker straw positions in MakeStrawGasSteps, so no conversion is needed
        const auto& s = step.startPosition();
        const auto& e = step.endPosition();
        const CLHEP::Hep3Vector pos(0.5*(s.x() + e.x()), 0.5*(s.y() + e.y()), 0.5*(s.z() + e.z()));
        const double r = pos.perp();
        const double edep = step.ionizingEdep();
        const auto id = step.strawId();
        if(edep > trk_step_emin_) {
          const double tf = foldTime(step.time());
          occ_->trk_time_folded_ecut_->Fill(tf, Weight);
          if(tf >= trk_hit_tmin_ && tf <= trk_hit_tmax_) occ_->trk_straw_hits_->Fill(id.uniqueStraw(), Weight);
        }
        for(int k = 0; k < 2; ++k) {
          const double w = (k == kOcc) ? Weight : Weight*edep;
          occ_->trk_plane_[k]->Fill(id.plane(), w);
          occ_->trk_panel_[k]->Fill(id.uniquePanel(), w);
          occ_->trk_straw_[k]->Fill(id.uniqueStraw(), w);
          occ_->trk_r_    [k]->Fill(r, w);
          occ_->trk_z_    [k]->Fill(pos.z(), w);
          occ_->trk_time_ [k]->Fill(step.time(), w);
          occ_->trk_time_folded_[k]->Fill(foldTime(step.time()), w);
          occ_->trk_rz_   [k]->Fill(pos.z(), r, w);
        }
      }
    }

    if(shower_col_) {
      GeomHandle<Calorimeter> calo;
      for(const auto& css : *shower_col_) {
        // CaloShowerStep::position() is crystal-local (CaloShowerStepMaker uses mu2eToCrystal); volumeG4ID() is the crystal ID
        const int cid = css.volumeG4ID();
        const int disk = calo->crystal(cid).diskID();
        const auto pos = det->toDetector(calo->crystalToMu2e(cid, css.position()));
        const double r = pos.perp();
        const double edep = css.energyDepG4();
        const double rcrys = det->toDetector(calo->crystal(cid).position()).perp();
        if(css.energyDepBirks() > calo_step_emin_) {
          const double tf = foldTime(css.time());
          occ_->calo_time_folded_ecut_->Fill(tf, Weight);
          if(tf >= calo_hit_tmin_ && tf <= calo_hit_tmax_) {
            occ_->calo_crystal_hits_->Fill(cid, Weight);
            occ_->calo_disk_rcrys_hits_[disk]->Fill(rcrys, Weight);
          }
        }
        for(int k = 0; k < 2; ++k) {
          const double w = (k == kOcc) ? Weight : Weight*edep;
          occ_->calo_disk_   [k]->Fill(disk, w);
          occ_->calo_crystal_[k]->Fill(cid , w);
          occ_->calo_r_      [k]->Fill(r, w);
          occ_->calo_z_      [k]->Fill(pos.z(), w);
          occ_->calo_time_   [k]->Fill(css.time(), w);
          occ_->calo_time_folded_[k]->Fill(foldTime(css.time()), w);
          occ_->calo_rz_     [k]->Fill(pos.z(), r, w);
          occ_->calo_disk_r_ [disk][k]->Fill(r, w);
          occ_->calo_disk_rcrys_[disk][k]->Fill(rcrys, w);
          occ_->calo_disk_xy_[disk][k]->Fill(pos.x(), pos.y(), w);
        }
      }
    }

    if(crv_col_) {
      GeomHandle<CosmicRayShield> crs;
      for(const auto& step : *crv_col_) {
        const double z = det->toDetector(CLHEP::Hep3Vector(0.5*(step.startPos().x() + step.endPos().x()),
                                                           0.5*(step.startPos().y() + step.endPos().y()),
                                                           0.5*(step.startPos().z() + step.endPos().z()))).z();
        const int sector = crs->getBar(step.barIndex()).id().getShieldNumber();
        // Start time, not the start/end midpoint: some CrvSteps have a NaN end time (e.g. MDC2025ac MuStopPileupCat dts)
        const double t = step.startTime();
        for(int k = 0; k < 2; ++k) {
          const double w = (k == kOcc) ? Weight : Weight*step.visibleEDep();
          occ_->crv_sector_[k]->Fill(sector, w);
          occ_->crv_bar_   [k]->Fill(step.barIndex().asInt(), w);
          occ_->crv_z_     [k]->Fill(z, w);
          occ_->crv_time_  [k]->Fill(t, w);
          occ_->crv_time_folded_[k]->Fill(foldTime(t), w);
        }
      }
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

    // CrvStep collection (optional)
    crv_col_ = nullptr;
    if(crv_step_tag_) {
      art::Handle<CrvStepCollection> crvH;
      event.getByLabel(*crv_step_tag_, crvH);
      crv_col_ = (crvH.isValid()) ? crvH.product() : nullptr;
    }

    if(debug_level_ > 0) watch_->StopTime("data retrieval");

    //-------------------------------------------------------------
    // Compute event-level info
    //-------------------------------------------------------------

    if(debug_level_ > 0) watch_->SetTime("event info computation");

    // Approximate hit counts: count steps, steps above an energy threshold, and steps above the threshold in the time window
    int trkNSteps[3] = {0, 0, 0};
    double totalTrackerE = 0.;
    if(straw_step_col_) {
      for(const auto& step : *straw_step_col_) {
        const double e = step.ionizingEdep();
        totalTrackerE += e;
        ++trkNSteps[0];
        if(e <= trk_step_emin_) continue;
        ++trkNSteps[1];
        if(step.time() >= trk_step_tmin_ && step.time() <= trk_step_tmax_) ++trkNSteps[2];
      }
    }

    int caloNSteps[3] = {0, 0, 0};
    double totalCaloE(0.), totalCaloEVis(0.);
    if(shower_col_) {
      for(const auto& css : *shower_col_) {
        totalCaloEVis += css.energyDepBirks(); // visible energy
        totalCaloE    += css.energyDepG4(); // total energy
        ++caloNSteps[0];
        if(css.energyDepBirks() <= calo_step_emin_) continue;
        ++caloNSteps[1];
        if(css.time() >= calo_step_tmin_ && css.time() <= calo_step_tmax_) ++caloNSteps[2];
      }
    }

    // Take the first primary as the main one for default histogramming
    info_.primsim = (primary_ && !primary_->primarySimParticles().empty()) ? &(*primary_->primarySimParticles().front()) : nullptr;
    info_.front_trk_sp = simTrkFrontStep(info_.primsim);
    info_.primsim_edep_vis = edepBySim(info_.primsim, true, true);
    info_.calo_total_edep = totalCaloE;
    info_.calo_total_edep_vis = totalCaloEVis;

    info_.weight_ = 1.; // can be used to apply event weights if needed

    //-------------------------------------------------------------
    // Set ouput tree data
    //-------------------------------------------------------------

    data_.weight = info_.weight_;
    data_.run = event.run();
    data_.subrun = event.subRun();
    data_.event = event.event();
    data_.event_calo_edep = totalCaloE;
    data_.event_calo_edep_vis = totalCaloEVis;
    data_.event_trk_edep  = totalTrackerE;
    data_.event_calo_nsteps           = caloNSteps[0];
    data_.event_calo_nsteps_ecut      = caloNSteps[1];
    data_.event_calo_nsteps_ecut_tcut = caloNSteps[2];
    data_.event_trk_nsteps            = trkNSteps[0];
    data_.event_trk_nsteps_ecut       = trkNSteps[1];
    data_.event_trk_nsteps_ecut_tcut  = trkNSteps[2];
    data_.nprim = (primary_) ? int(primary_->primarySimParticles().size()) : 0;
    if(data_.nprim >= kMaxPrimaries) throw cet::exception("Analysis") << "Too many primary particles!"
                                                                           << "N(primaries) = " << data_.nprim
                                                                           << " > " << kMaxPrimaries;

    // Store information for each primary
    GeomHandle<DetectorSystem> det;
    GeomHandle<BFieldManager> bfmgr;
    GlobalConstantsHandle<ParticleDataList> pdt;
    for(int iprim = 0; iprim < data_.nprim; ++iprim) {
      const auto* sim = &(*primary_->primarySimParticles().at(iprim));
      const auto pos = sim->startPosition();
      const auto mom = sim->startMomentum();
      data_.prim_start_x      [iprim] = pos.x();
      data_.prim_start_y      [iprim] = pos.y();
      data_.prim_start_z      [iprim] = pos.z();
      data_.prim_start_px     [iprim] = mom.x();
      data_.prim_start_py     [iprim] = mom.y();
      data_.prim_start_pz     [iprim] = mom.z();
      data_.prim_start_e      [iprim] = mom.e();
      data_.prim_start_m      [iprim] = mom.m();
      data_.prim_start_pdg    [iprim] = sim->pdgId();
      data_.prim_calo_edep    [iprim] = edepBySim(sim, true);
      data_.prim_calo_edep_vis[iprim] = (iprim == 0) ? info_.primsim_edep_vis : edepBySim(sim, true, true);
      const auto front_trk_sp = (iprim == 0) ? info_.front_trk_sp : simTrkFrontStep(sim);
      data_.prim_trk_front_energy[iprim] = energyFromStepPoint(sim, front_trk_sp);
      if(front_trk_sp) {
        const auto trk_mom = front_trk_sp->momentum();
        const auto trk_pos = front_trk_sp->position();
        const float p = trk_mom.mag();
        const float cos = (p > 0.) ? trk_mom.z() / p : 0.f;
        // Transverse helix center from the local Bz and the particle charge:
        // center = pos + (py, -px) / (c q Bz), radius = pT / |c q Bz|; -1 if q*Bz = 0
        const double bz = bfmgr->getBField(trk_pos).z()*CLHEP::tesla; // BFieldManager returns Tesla
        const double q  = pdt->particle(sim->pdgId()).charge();
        float r_center = -1.f;
        float radius   = -1.f;
        if(q*bz != 0.) {
          const double scale = 1./(CLHEP::c_light*q*bz);
          const auto det_pos = det->toDetector(trk_pos);
          const double x_c = det_pos.x() + trk_mom.y()*scale;
          const double y_c = det_pos.y() - trk_mom.x()*scale;
          r_center = std::sqrt(x_c*x_c + y_c*y_c);
          radius   = trk_mom.perp()*std::fabs(scale);
        }
        data_.prim_trk_front_p       [iprim] = p;
        data_.prim_trk_front_r_center[iprim] = r_center;
        data_.prim_trk_front_radius  [iprim] = radius;
        data_.prim_trk_front_time    [iprim] = front_trk_sp->time();
        data_.prim_trk_front_cos     [iprim] = cos;
      } else {
        data_.prim_trk_front_p       [iprim] =  0.f;
        data_.prim_trk_front_r_center[iprim] = -1.f;
        data_.prim_trk_front_radius  [iprim] = -1.f;
        data_.prim_trk_front_time    [iprim] = -1.f;
        data_.prim_trk_front_cos     [iprim] =  0.f;
      }
    }
    tree_->Fill();

    //-------------------------------------------------------------
    // Increment summary data
    //-------------------------------------------------------------

    total_events_ += info_.weight_;
    total_calo_edep_ += info_.weight_ * totalCaloE;
    total_calo_edep_vis_ += info_.weight_ * totalCaloEVis;
    if(totalCaloEVis > 50.) events_above_50_mev_ += info_.weight_;
    total_tracker_edep_ += info_.weight_ * totalTrackerE;
    for(int i = 0; i < 3; ++i) {
      total_calo_nsteps_[i] += info_.weight_ * caloNSteps[i];
      total_trk_nsteps_ [i] += info_.weight_ * trkNSteps [i];
    }

    if(debug_level_ > 0) watch_->StopTime("event info computation");

    //-------------------------------------------------------------
    // Fill histograms
    //-------------------------------------------------------------

    fillHistograms(hists_[0], info_.weight_);
    if(totalCaloEVis >  1.) fillHistograms(hists_[1], info_.weight_);
    if(totalCaloEVis > 10.) fillHistograms(hists_[2], info_.weight_);
    if(totalCaloEVis > 50.) fillHistograms(hists_[3], info_.weight_);
    fillOccupancyHistograms(info_.weight_);
    if(debug_level_ > 0) watch_->StopTime(__func__);
  }


  void EdepAna::endJob() {
    if(debug_level_ > 0) std::cout << "[EdepAna::" << __func__ << "]\n" << *watch_ << std::endl;
    if(occ_ && occ_->h_ngen_) occ_->h_ngen_->SetBinContent(1, ngen_);

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
      TH1* h_ref = (hists_[0]) ? hists_[0]->h_total_calo_energy_vis_ : nullptr;
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
      TH1* h_ref = (hists_[0]) ? hists_[0]->h_total_calo_energy_vis_ : nullptr;
      const double eff = h->GetEntries() > 0 && h_ref ? h->GetEntries() * 1./ h_ref->GetEntries() : 0.;
      double mpv, fwhm;
      get_landau_seed(h, mpv, fwhm);
      std::cout << std::format("Primary Edep - tracker front StepPointMC energy distribution: MPV = {:.2f} MeV, FWHM = {:.2f} MeV, efficiency = {:.4g}", mpv, fwhm, eff) << std::endl;
    }

    //-------------------------------------------------------------
    // Report summary data
    //-------------------------------------------------------------

    const double averageCaloEdep    = (total_events_ > 0) ? total_calo_edep_    / total_events_ : -1.;
    const double averageCaloEdepGen = (ngen_         > 0) ? total_calo_edep_    / ngen_         : -1.;
    const double averageTrkEdep     = (total_events_ > 0) ? total_tracker_edep_ / total_events_ : -1.;
    const double averageTrkEdepGen  = (ngen_         > 0) ? total_tracker_edep_ / ngen_         : -1.;
    const double eventRate          = (ngen_         > 0) ? total_events_       / ngen_         : -1.;

    std::cout
      << "EdepAna summary:\n"
      << "  Saw " << total_events_ << " events" << " (" << ngen_ << " gen events) --> output rate = "
      << eventRate << " events / gen event" << std::endl
      << "  Average calo energy deposition per event: " << averageCaloEdep << " MeV" << std::endl
      << "  Average calo energy deposition per gen event: " << averageCaloEdepGen << " MeV" << std::endl
      << "  Events with calo Edep > 50 MeV: " << events_above_50_mev_ << std::endl
      << "  Average tracker energy deposition per event: " << averageTrkEdep << " MeV" << std::endl
      << "  Average tracker energy deposition per gen event: " << averageTrkEdepGen << " MeV" << std::endl;

    // Approximate hit counts from step counts
    auto printCounts = [&](const char* det, const double* counts, const double emin, const double tmin, const double tmax) {
      std::cout << std::format("  {} step counts (E > {:g} MeV, {:g} <= t <= {:g} ns):", det, emin, tmin, tmax) << std::endl;
      const char* labels[3] = {"all", "E cut", "E + t cut"};
      for(int i = 0; i < 3; ++i) {
        const double perEvent = (total_events_ > 0) ? counts[i] / total_events_ : -1.;
        const double perGen   = (ngen_         > 0) ? counts[i] / ngen_         : -1.;
        std::cout << std::format("    {:<10}: total = {:<12g} per event = {:<12.5g} per gen event = {:.5g}",
                                 labels[i], counts[i], perEvent, perGen) << std::endl;
      }
    };
    printCounts("Calo"   , total_calo_nsteps_, calo_step_emin_, calo_step_tmin_, calo_step_tmax_);
    printCounts("Tracker", total_trk_nsteps_ , trk_step_emin_ , trk_step_tmin_ , trk_step_tmax_ );
  }

}

using mu2e::EdepAna;
DEFINE_ART_MODULE(EdepAna)
