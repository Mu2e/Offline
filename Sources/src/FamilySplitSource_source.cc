// FamilySplitSource: one pass over cosmic stage-2 dts files made with Geant4 secondary splitting
// (Mu2eG4 physics.secondaryBiasing), writing for every input event
//   j = 0      the event with every split ("forced") family removed (the "stripped" event), and
//   j = 1..k   one event per forced family, holding the rest of the event plus that family only.
// Output event number = input event number * stride + j.  Configuration:
// Production/JobConfig/cosmic/S2FamilySplit.fcl.
//
// Why: splitting makes N alternative secondaries at one interaction, all in one event. Downstream
// stages count events, so two alternatives that both reach the tracker would be counted once and
// the event would hold two signal-like electrons at once. One family per event removes both problems
// and lets N be large.
//
// A forced particle is a direct daughter of a muon, created by one of `processes`, born inside one of
// `regions` (cylinders along z, Mu2e coordinates); its family is itself and all its descendants. Geant4
// secondary splitting made `factor` alternatives of it, each of weight 1/factor.
//
// Every output event carries an EventWeight (label `weightLabel`), times the base weight (1, or the
// forced-decay weight below):
//   family event:   1 / factor of the region its forced particle was born in;
//   stripped event: 1 - (sum over the file's forced families of 1/factor) / (events kept from the file).
// The weights add up to the number of events kept, so the dataset keeps the livetime of its input
// (CosmicLivetime is copied unchanged). Anything a family carries (the particles it makes) keeps its
// expectation exactly. What the rest of an input event holds (the muon, unsplit particles) is written in
// the stripped event and again in each of its family events, with total weight 1 - S/E + s, where s is
// the sum of that event's family weights: right on average over the file, biased for quantities that
// correlate with s. (The exact per-event stripped weight, 1 - s, can be negative.) In the cosmic
// configuration S/E is about 0.02 and the muon count agrees with an unsplit sample to 0.1%; a soft-kill
// exemption that let soft split copies survive in one region (S/E 0.028) made the count of muons
// crossing the tracker 0.991 of it.
//
// Muon decays in flight (optional, muonDecayWindow = T > 0, proper time in ns). They are rare, so they
// can come from a second stage-2 stream in which Mu2eG4 physics.muonMaxPreAssignedDecayProperTime = T makes
// every muon decay within T (a stage-2 muon starts at proper time 0). Only stage-2 primary muons that
// decay in flight count (end momentum > 0; for a decay at rest the SimParticle records the time the muon
// stopped, not the decay time, and its electron is below 53 MeV anyway):
//   main stream   (forcedDecayStream: false): events whose primary muon decays in flight before T are
//                                              dropped;
//   forced stream (forcedDecayStream: true):  only events whose primary muon decays in flight are kept, and
//                                              every weight is multiplied by 1 - exp(-T / muonLifetime).
// The two streams then add up to the unbiased sample, each normalized to its own livetime. This is exact
// for one primary muon per event whatever happens to it (decay in flight or at rest, capture, leaving the
// world); events with several primary muons are approximate and are counted in the summary.
//
// SubRuns are all written, also those left without events, so no livetime is lost.
//
// Products are re-created with their stage-2 names, so the job's process_name must be the stage-2 one
// (`inputProcessName`): art derives a ProductID from the product name, so the art::Ptr in the copied
// products still point to the right collections without remapping. Detector step collections
// (StrawGas, CaloShower, Crv, Surface, virtual-detector StepPointMC) and MC trajectories are filtered;
// all other products are copied. The input TriggerResults is not copied (art makes its own). Any other
// branch in the input stops the job.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"

#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/FileFormatVersion.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Provenance/SubRunAuxiliary.h"
#include "canvas/Persistency/Provenance/canonicalProductName.h"
#include "canvas/Utilities/TypeID.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableFragment.h"
#include "art/Framework/Core/FileBlock.h"

#include "Offline/DataProducts/inc/StageNormalization.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"

namespace mu2e {

  namespace {
    // art branch names are canonical product names ending in '.'
    std::string branchName(std::string const& canonical) {
      return (!canonical.empty() && canonical.back() == '.') ? canonical : canonical + ".";
    }

    // One art product branch of the input Events or SubRuns tree.
    template <typename T>
    struct Branch {
      std::string label, instance, name;
      TBranch* branch = nullptr;
      art::Wrapper<T>* wrapper = nullptr;
      Branch(std::string l, std::string i) : label(std::move(l)), instance(std::move(i)) {}
      void setName(std::string const& process) {
        name = branchName(art::canonicalProductName(art::TypeID(typeid(T)).friendlyClassName(), label, instance, process));
      }
      bool tryAttach(TTree* tree) {
        branch = tree->GetBranch(name.c_str());
        wrapper = nullptr;
        if (branch != nullptr) { branch->SetAddress(&wrapper); }
        return branch != nullptr;
      }
      void attach(TTree* tree) {
        branch = tree->GetBranch(name.c_str());
        if (branch == nullptr) {
          throw cet::exception("FamilySplitSource") << "no branch " << name << " in tree " << tree->GetName() << "\n";
        }
        wrapper = nullptr;
        branch->SetAddress(&wrapper);
      }
      // for optional products: the branch may exist with the product marked not present
      std::optional<T> readIfPresent(Long64_t entry) {
        if (branch->GetEntry(entry) <= 0 || wrapper == nullptr) {
          throw cet::exception("FamilySplitSource") << "cannot read " << name << " entry " << entry << "\n";
        }
        if (!wrapper->isPresent()) { return std::nullopt; }
        return *wrapper->product();
      }
      T const& read(Long64_t entry) {
        if (branch->GetEntry(entry) <= 0 || wrapper == nullptr || !wrapper->isPresent()) {
          throw cet::exception("FamilySplitSource") << "cannot read " << name << " entry " << entry << "\n";
        }
        return *wrapper->product();
      }
    };
  } // namespace

  class FamilySplitSourceDetail {
  public:
    struct Region {
      fhicl::Atom<double> xc{fhicl::Name("xc"), fhicl::Comment("axis x [mm], Mu2e coordinates")};
      fhicl::Atom<double> yc{fhicl::Name("yc"), 0.};
      fhicl::Atom<double> rmin{fhicl::Name("rmin"), 0.};
      fhicl::Atom<double> rmax{fhicl::Name("rmax")};
      fhicl::Atom<double> zmin{fhicl::Name("zmin"), fhicl::Comment("[mm], Mu2e coordinates")};
      fhicl::Atom<double> zmax{fhicl::Name("zmax")};
      fhicl::Atom<double> factor{fhicl::Name("factor"), fhicl::Comment("splitting factor used in this region")};
    };
    struct Config {
      fhicl::Sequence<std::string> processes{fhicl::Name("processes"),
          fhicl::Comment("creation codes of forced daughters of muons")};
      fhicl::Sequence<fhicl::Table<Region>> regions{fhicl::Name("regions")};
      fhicl::Atom<unsigned> stride{fhicl::Name("stride"), 1000};
      fhicl::Atom<std::string> weightLabel{fhicl::Name("weightLabel"), "FamilySplit"};
      fhicl::Atom<std::string> inputProcessName{fhicl::Name("inputProcessName"),
          fhicl::Comment("process name of the stage-2 products; must equal this job's process_name"), "CosmicS2"};
      fhicl::Atom<double> muonDecayWindow{fhicl::Name("muonDecayWindow"),
          fhicl::Comment("[ns] proper time; 0: muon decays not treated. Otherwise see forcedDecayStream"), 0.};
      fhicl::Atom<bool> forcedDecayStream{fhicl::Name("forcedDecayStream"),
          fhicl::Comment("false: drop events whose stage-2 primary muon decays in flight before muonDecayWindow.\n"
                         "true: the input was made with muonMaxPreAssignedDecayProperTime = muonDecayWindow; keep only\n"
                         "events whose primary muon decays in flight, weight times 1 - exp(-muonDecayWindow / muonLifetime)"),
          false};
      fhicl::Atom<double> muonLifetime{fhicl::Name("muonLifetime"), fhicl::Comment("[ns], as in Geant4"), 2196.98};
      fhicl::Atom<int> primaryStage{fhicl::Name("primaryStage"),
          fhicl::Comment("simStage of the stage-2 primaries (S2Inputs.simStageOverride)"), 1};
      fhicl::Atom<int> verbosity{fhicl::Name("verbosity"), 0};
    };
    using Parameters = fhicl::TableFragment<Config>;   // art::Source adds fileNames, maxEvents, maxSubRuns

    FamilySplitSourceDetail(Parameters const& conf, art::ProductRegistryHelper& rh, art::SourceHelper const& sh);
    FamilySplitSourceDetail(FamilySplitSourceDetail const&) = delete;
    FamilySplitSourceDetail& operator=(FamilySplitSourceDetail const&) = delete;
    ~FamilySplitSourceDetail();

    void readFile(std::string const& filename, art::FileBlock*& fb);
    bool readNext(art::RunPrincipal* const& inR, art::SubRunPrincipal* const& inSR,
                  art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR, art::EventPrincipal*& outE);
    void closeCurrentFile();

  private:
    struct Cyl { double xc, yc, rmin, rmax, zmin, zmax, factor; };
    art::SourceHelper const& sh_;
    std::vector<std::string> processes_;
    std::vector<Cyl> regions_;
    unsigned stride_;
    std::string weightLabel_, process_;
    double window_, lifetime_, baseWeight_;
    bool forcedStream_;
    int primaryStage_;
    int verbosity_;
    art::ProductID simsPID_;

    // input
    std::unique_ptr<TFile> file_;
    TTree* events_ = nullptr;
    art::EventAuxiliary* eventAux_ = nullptr;
    Branch<SimParticleCollection> sims_{"compressDetStepMCs", ""};
    Branch<StrawGasStepCollection> sgs_{"compressDetStepMCs", ""};
    Branch<CaloShowerStepCollection> css_{"compressDetStepMCs", ""};
    Branch<CrvStepCollection> crv_{"compressDetStepMCs", ""};
    Branch<SurfaceStepCollection> surf_{"compressDetStepMCs", ""};
    Branch<StepPointMCCollection> vd_{"compressDetStepMCs", "virtualdetector"};
    Branch<MCTrajectoryCollection> traj_{"compressDetStepMCs", ""};
    Branch<GenParticleCollection> gen_{"compressDetStepMCs", ""};
    Branch<StatusG4> status_{"g4run", ""};
    Branch<StepPointMCCollection> crvStage1_{"CosmicResampler", "crvStage1"};
    Branch<std::vector<art::EventID>> eids_{"CosmicResampler", ""};
    Branch<CosmicLivetime> evLive_{"CosmicResampler", ""};
    Branch<PrimaryParticle> primary_{"FindMCPrimary", ""};
    Branch<CosmicLivetime> srLive_{"CosmicResampler", "mixed"};
    Branch<GenEventCount> srCount_{"genCounter", ""};                // optional
    Branch<StageNormalization> srNorm_{"s2Counter", ""};             // optional
    std::string triggerResultsName_;
    struct SubRunProducts {
      art::Timestamp begin;
      CosmicLivetime livetime;
      std::optional<GenEventCount> count;
      std::optional<StageNormalization> norm;
      bool written = false;
    };
    std::map<art::SubRunID, SubRunProducts> subRunProducts_;

    // state
    Long64_t nEntries_ = 0, entry_ = -1;
    std::vector<char> keep_;                                  // per input entry, after the muon-decay rule
    std::vector<std::pair<size_t, double>> forced_;          // key, factor of the current input event
    std::unordered_set<size_t> roots_;                        // the forced keys of the current input event
    std::unordered_map<size_t, long> family_;                 // particle key -> forced key or -1
    size_t j_ = 1;                                            // next output index of the current input event
    double stripWeight_ = 1.;
    art::RunNumber_t lastRun_ = 0;
    art::SubRunID lastSubRun_;
    bool checkedProcess_ = false;

    // totals
    unsigned long nIn_ = 0, nOut_ = 0, nFamilies_ = 0;
    unsigned long nDecayDropped_ = 0, nMultiMuon_ = 0, nLateForced_ = 0, nEmptySubRuns_ = 0;
    std::map<size_t, unsigned long> multiplicity_;

    bool keepByDecay(SimParticleCollection const& sims);
    void writeSubRun(art::SubRunID const& id, art::Timestamp const& ts,
                     art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR);
    double regionFactor(SimParticle const& p) const;
    void findForced(SimParticleCollection const& sims, std::vector<std::pair<size_t, double>>& out) const;
    long familyOf(art::Ptr<SimParticle> const& ptr, SimParticleCollection const& sims);
    void checkBranches();
  };

  FamilySplitSourceDetail::FamilySplitSourceDetail(Parameters const& conf, art::ProductRegistryHelper& rh,
                                                   art::SourceHelper const& sh)
    : sh_(sh)
    , processes_(conf().processes())
    , stride_(conf().stride())
    , weightLabel_(conf().weightLabel())
    , process_(conf().inputProcessName())
    , window_(conf().muonDecayWindow())
    , lifetime_(conf().muonLifetime())
    , baseWeight_(1.)
    , forcedStream_(conf().forcedDecayStream())
    , primaryStage_(conf().primaryStage())
    , verbosity_(conf().verbosity()) {
    if (window_ < 0. || !(lifetime_ > 0.) || (forcedStream_ && window_ == 0.)) {
      throw cet::exception("CONFIG") << "FamilySplitSource: need muonDecayWindow >= 0, muonLifetime > 0, and "
                                     << "muonDecayWindow > 0 with forcedDecayStream\n";
    }
    if (forcedStream_) { baseWeight_ = -std::expm1(-window_ / lifetime_); }
    for (auto const& r : conf().regions()) {
      if (!(r.factor() >= 1.)) {
        throw cet::exception("CONFIG") << "FamilySplitSource: region factor must be >= 1\n";
      }
      regions_.push_back({r.xc(), r.yc(), r.rmin(), r.rmax(), r.zmin(), r.zmax(), r.factor()});
    }
    if (regions_.empty() || processes_.empty() || stride_ < 2) {
      throw cet::exception("CONFIG") << "FamilySplitSource: need regions, processes and stride >= 2\n";
    }
    simsPID_ = art::ProductID{art::canonicalProductName(art::TypeID(typeid(SimParticleCollection)).friendlyClassName(),
                                                        "compressDetStepMCs", "", process_)};
    sims_.setName(process_); sgs_.setName(process_); css_.setName(process_); crv_.setName(process_);
    surf_.setName(process_); vd_.setName(process_); traj_.setName(process_); gen_.setName(process_);
    status_.setName(process_); crvStage1_.setName(process_); eids_.setName(process_); evLive_.setName(process_);
    primary_.setName(process_); srLive_.setName(process_); srCount_.setName(process_);
    srNorm_.setName(process_);
    triggerResultsName_ = branchName(art::canonicalProductName("art::TriggerResults", "TriggerResults", "", process_));

    rh.reconstitutes<SimParticleCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<StrawGasStepCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<CaloShowerStepCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<CrvStepCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<SurfaceStepCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<StepPointMCCollection, art::InEvent>("compressDetStepMCs", "virtualdetector");
    rh.reconstitutes<MCTrajectoryCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<GenParticleCollection, art::InEvent>("compressDetStepMCs");
    rh.reconstitutes<StatusG4, art::InEvent>("g4run");
    rh.reconstitutes<StepPointMCCollection, art::InEvent>("CosmicResampler", "crvStage1");
    rh.reconstitutes<std::vector<art::EventID>, art::InEvent>("CosmicResampler");
    rh.reconstitutes<CosmicLivetime, art::InEvent>("CosmicResampler");
    rh.reconstitutes<PrimaryParticle, art::InEvent>("FindMCPrimary");
    rh.reconstitutes<EventWeight, art::InEvent>(weightLabel_);
    rh.reconstitutes<CosmicLivetime, art::InSubRun>("CosmicResampler", "mixed");
    rh.reconstitutes<GenEventCount, art::InSubRun>("genCounter");
    rh.reconstitutes<StageNormalization, art::InSubRun>("s2Counter");
  }

  FamilySplitSourceDetail::~FamilySplitSourceDetail() {
    std::cout << "FamilySplitSource: input events " << nIn_ << " output events " << nOut_
              << " family events " << nFamilies_ << " | families per input event";
    for (auto const& [k, n] : multiplicity_) { std::cout << " " << k << ":" << n; }
    std::cout << std::endl;
    if (window_ > 0.) {
      std::cout << "FamilySplitSource: " << (forcedStream_ ? "forced-decay" : "main") << " stream, muon decay window "
                << window_ << " ns: " << nDecayDropped_ << " input events dropped, base weight " << baseWeight_
                << ", events with several primary muons (approximate) " << nMultiMuon_;
      if (forcedStream_) { std::cout << ", primary muon decays later than the window " << nLateForced_; }
      std::cout << std::endl;
    }
    std::cout << "FamilySplitSource: subruns written without events " << nEmptySubRuns_ << std::endl;
  }

  bool FamilySplitSourceDetail::keepByDecay(SimParticleCollection const& sims) {
    if (window_ <= 0.) { return true; }
    unsigned nMu = 0;
    double first = std::numeric_limits<double>::max();     // earliest decay of a primary muon, proper time
    for (auto const& kv : sims) {
      SimParticle const& p = kv.second;
      if (std::abs(p.pdgId()) != 13 || static_cast<int>(p.simStage()) != primaryStage_ ||
          p.creationCode().name() != "mu2ePrimary") { continue; }
      ++nMu;
      const std::string stop = p.stoppingCode().name();
      if ((stop == "Decay" || stop == "DecayWithSpin") && p.endMomentum().vect().mag() > 0.) {
        first = std::min(first, static_cast<double>(p.endProperTime() - p.startProperTime()));
      }
    }
    if (nMu > 1) { ++nMultiMuon_; }
    bool keep;
    if (!forcedStream_) {
      keep = !(first < window_);
    } else {
      keep = first <= window_ * (1. + 1.e-6);
      if (!keep && first < std::numeric_limits<double>::max()) { ++nLateForced_; }
    }
    if (!keep) { ++nDecayDropped_; }
    return keep;
  }

  void FamilySplitSourceDetail::writeSubRun(art::SubRunID const& id, art::Timestamp const& ts,
                                            art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR) {
    auto it = subRunProducts_.find(id);
    if (it == subRunProducts_.end()) {
      throw cet::exception("FamilySplitSource") << "no SubRuns entry for " << id << "\n";
    }
    if (it->second.written) {
      throw cet::exception("FamilySplitSource") << "events of " << id
                                                << " come back after another subrun; not handled\n";
    }
    if (id.run() != lastRun_) {
      outR = sh_.makeRunPrincipal(id.run(), ts);
      lastRun_ = id.run();
    }
    outSR = sh_.makeSubRunPrincipal(id.run(), id.subRun(), ts);
    // copied unchanged: the livetime is that of the input, so the EventWeights normalize to it. The
    // StageNormalization still counts input events (nPassed), not the events written here.
    art::put_product_in_principal(std::make_unique<CosmicLivetime>(it->second.livetime), *outSR, "CosmicResampler", "mixed");
    if (it->second.count) {
      art::put_product_in_principal(std::make_unique<GenEventCount>(*it->second.count), *outSR, "genCounter");
    }
    if (it->second.norm) {
      art::put_product_in_principal(std::make_unique<StageNormalization>(*it->second.norm), *outSR, "s2Counter");
    }
    it->second.written = true;
    lastSubRun_ = id;
  }

  double FamilySplitSourceDetail::regionFactor(SimParticle const& p) const {
    auto const x = p.startPosition();
    for (auto const& c : regions_) {
      const double dx = x.x() - c.xc, dy = x.y() - c.yc;
      const double r = std::sqrt(dx * dx + dy * dy);
      if (r >= c.rmin && r < c.rmax && x.z() >= c.zmin && x.z() < c.zmax) { return c.factor; }
    }
    return 0.;
  }

  void FamilySplitSourceDetail::findForced(SimParticleCollection const& sims,
                                           std::vector<std::pair<size_t, double>>& out) const {
    out.clear();
    for (auto const& kv : sims) {
      SimParticle const& p = kv.second;
      if (!p.hasParent() || p.parent().id() != simsPID_) { continue; }
      SimParticle const* par = sims.getOrNull(cet::map_vector_key(p.parent().key()));
      if (par == nullptr || std::abs(par->pdgId()) != 13) { continue; }
      const std::string code = p.creationCode().name();
      bool ok = false;
      for (auto const& c : processes_) { if (code == c) { ok = true; break; } }
      if (!ok) { continue; }
      const double f = regionFactor(p);
      if (f > 0.) { out.emplace_back(kv.first.asUint(), f); }
    }
  }

  long FamilySplitSourceDetail::familyOf(art::Ptr<SimParticle> const& ptr, SimParticleCollection const& sims) {
    if (ptr.id() != simsPID_) {
      throw cet::exception("FamilySplitSource") << "a step points to product " << ptr.id()
                                                << ", not to the SimParticles " << simsPID_ << "\n";
    }
    std::vector<size_t> path;
    long result = -1;
    size_t k = ptr.key();
    while (true) {
      auto it = family_.find(k);
      if (it != family_.end()) { result = it->second; break; }
      path.push_back(k);
      SimParticle const* p = sims.getOrNull(cet::map_vector_key(k));
      if (p == nullptr) {
        throw cet::exception("FamilySplitSource") << "SimParticle key " << k << " not in the collection\n";
      }
      if (std::abs(p->pdgId()) == 13) { result = -1; break; }
      if (roots_.count(k) != 0) { result = static_cast<long>(k); break; }
      if (!p->hasParent() || p->parent().id() != simsPID_) { result = -1; break; }
      k = p->parent().key();
    }
    for (size_t q : path) { family_[q] = result; }
    return result;
  }

  void FamilySplitSourceDetail::checkBranches() {
    std::set<std::string> known{"EventAuxiliary", triggerResultsName_, sims_.name, sgs_.name, css_.name, crv_.name,
                                surf_.name, vd_.name, traj_.name, gen_.name, status_.name, crvStage1_.name,
                                eids_.name, evLive_.name, primary_.name};
    for (auto* obj : *events_->GetListOfBranches()) {
      if (known.count(obj->GetName()) == 0) {
        throw cet::exception("FamilySplitSource") << "input Events branch " << obj->GetName()
                                                  << " is not handled; it would be lost\n";
      }
    }
  }

  void FamilySplitSourceDetail::readFile(std::string const& filename, art::FileBlock*& fb) {
    file_.reset(TFile::Open(filename.c_str()));
    if (!file_ || file_->IsZombie()) {
      throw cet::exception("FamilySplitSource") << "cannot open " << filename << "\n";
    }
    events_ = file_->Get<TTree>("Events");
    TTree* subRuns = file_->Get<TTree>("SubRuns");
    TTree* runs = file_->Get<TTree>("Runs");
    if (!events_ || !subRuns || !runs) {
      throw cet::exception("FamilySplitSource") << filename << ": missing Events, SubRuns or Runs tree\n";
    }
    checkBranches();
    if (runs->GetListOfBranches()->GetEntries() != 1) {
      throw cet::exception("FamilySplitSource") << filename << ": Runs tree has products, not handled\n";
    }
    {
      std::set<std::string> known{"SubRunAuxiliary", srLive_.name, srCount_.name, srNorm_.name};
      for (auto* obj : *subRuns->GetListOfBranches()) {
        if (known.count(obj->GetName()) == 0) {
          throw cet::exception("FamilySplitSource") << "input SubRuns branch " << obj->GetName()
                                                    << " is not handled; it would be lost\n";
        }
      }
    }
    eventAux_ = nullptr;
    events_->GetBranch("EventAuxiliary")->SetAddress(&eventAux_);
    sims_.attach(events_); sgs_.attach(events_); css_.attach(events_); crv_.attach(events_);
    surf_.attach(events_); vd_.attach(events_); traj_.attach(events_); gen_.attach(events_);
    status_.attach(events_); crvStage1_.attach(events_); eids_.attach(events_); evLive_.attach(events_);
    primary_.attach(events_);

    // subrun products
    subRunProducts_.clear();
    art::SubRunAuxiliary* srAux = nullptr;
    subRuns->GetBranch("SubRunAuxiliary")->SetAddress(&srAux);
    srLive_.attach(subRuns);
    const bool hasCount = srCount_.tryAttach(subRuns);
    const bool hasNorm = srNorm_.tryAttach(subRuns);
    for (Long64_t i = 0; i < subRuns->GetEntries(); ++i) {
      subRuns->GetBranch("SubRunAuxiliary")->GetEntry(i);
      SubRunProducts sp{srAux->beginTime(), srLive_.read(i), std::nullopt, std::nullopt};
      if (hasCount) { sp.count = srCount_.readIfPresent(i); }
      if (hasNorm) { sp.norm = srNorm_.readIfPresent(i); }
      subRunProducts_.emplace(srAux->id(), sp);
    }
    subRuns->ResetBranchAddresses();
    if (subRunProducts_.count(lastSubRun_) != 0) {
      throw cet::exception("FamilySplitSource") << filename << ": " << lastSubRun_
                                                << " continues from the previous file; not handled\n";
    }

    // first pass over the SimParticles only: the muon-decay rule, and the weight of the stripped events of
    // this file (over the events kept)
    nEntries_ = events_->GetEntries();
    keep_.assign(nEntries_, 0);
    double sumInv = 0.;
    Long64_t nKept = 0;
    std::vector<std::pair<size_t, double>> f;
    for (Long64_t i = 0; i < nEntries_; ++i) {
      SimParticleCollection const& sims = sims_.read(i);
      if (!keepByDecay(sims)) { continue; }
      keep_[i] = 1;
      ++nKept;
      findForced(sims, f);
      for (auto const& x : f) { sumInv += 1. / x.second; }
    }
    stripWeight_ = nKept > 0 ? 1. - sumInv / static_cast<double>(nKept) : 1.;
    entry_ = -1;
    j_ = 1;
    forced_.clear();
    std::cout << "FamilySplitSource: " << filename << ": " << nEntries_ << " events, " << nKept
              << " kept, sum of family weights " << sumInv << ", stripped-event weight " << stripWeight_
              << ", base weight " << baseWeight_ << std::endl;
    fb = new art::FileBlock(art::FileFormatVersion(1, "FamilySplitSource"), filename);
  }

  void FamilySplitSourceDetail::closeCurrentFile() {
    if (events_) { events_->ResetBranchAddresses(); }
    events_ = nullptr;
    file_.reset();
  }

  bool FamilySplitSourceDetail::readNext(art::RunPrincipal* const&, art::SubRunPrincipal* const&,
                                         art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR,
                                         art::EventPrincipal*& outE) {
    if (j_ > forced_.size()) {
      do { ++entry_; } while (entry_ < nEntries_ && !keep_[entry_]);
      if (entry_ >= nEntries_) {
        // subruns left without events are written alone, so that their livetime is kept
        for (auto const& [id, sp] : subRunProducts_) {
          if (!sp.written) {
            writeSubRun(id, sp.begin, outR, outSR);
            ++nEmptySubRuns_;
            return true;
          }
        }
        return false;
      }
      if (events_->GetBranch("EventAuxiliary")->GetEntry(entry_) <= 0 || eventAux_ == nullptr) {
        throw cet::exception("FamilySplitSource") << "cannot read EventAuxiliary entry " << entry_ << "\n";
      }
      // read every product of this input event once; the wrappers stay valid until the next entry
      sims_.read(entry_); sgs_.read(entry_); css_.read(entry_); crv_.read(entry_); surf_.read(entry_);
      vd_.read(entry_); traj_.read(entry_); gen_.read(entry_); status_.read(entry_); crvStage1_.read(entry_);
      eids_.read(entry_); evLive_.read(entry_); primary_.read(entry_);
      findForced(*sims_.wrapper->product(), forced_);
      if (forced_.size() + 1 > stride_) {
        throw cet::exception("FamilySplitSource") << forced_.size() << " families in one event, stride "
                                                  << stride_ << " too small\n";
      }
      roots_.clear();
      for (auto const& x : forced_) { roots_.insert(x.first); }
      family_.clear();
      j_ = 0;
      ++nIn_;
      ++multiplicity_[forced_.size()];
    }
    art::EventAuxiliary const& aux = *eventAux_;
    const art::Timestamp ts = aux.time();
    if (aux.subRunID() != lastSubRun_) { writeSubRun(aux.subRunID(), ts, outR, outSR); }
    const unsigned long ev = static_cast<unsigned long>(aux.event()) * stride_ + j_;
    if (ev > 0xffffffffUL) {
      throw cet::exception("FamilySplitSource") << "event number " << aux.event() << " * stride overflows\n";
    }
    outE = sh_.makeEventPrincipal(aux.run(), aux.subRun(), static_cast<art::EventNumber_t>(ev), ts, false,
                                  art::EventAuxiliary::Mc);
    if (!checkedProcess_) {
      if (outE->processConfiguration().processName() != process_) {
        throw cet::exception("CONFIG") << "FamilySplitSource: process_name must be " << process_
                                       << " so that copied art::Ptr stay valid, not "
                                       << outE->processConfiguration().processName() << "\n";
      }
      checkedProcess_ = true;
    }

    SimParticleCollection const& sims = *sims_.wrapper->product();
    const long kept = j_ == 0 ? -2 : static_cast<long>(forced_[j_ - 1].first);
    auto wanted = [&](art::Ptr<SimParticle> const& p) {
      const long f = familyOf(p, sims);
      return f == -1 || f == kept;
    };
    auto sgs = std::make_unique<StrawGasStepCollection>();
    for (auto const& s : *sgs_.wrapper->product()) { if (wanted(s.simParticle())) sgs->push_back(s); }
    auto css = std::make_unique<CaloShowerStepCollection>();
    for (auto const& s : *css_.wrapper->product()) { if (wanted(s.simParticle())) css->push_back(s); }
    auto crv = std::make_unique<CrvStepCollection>();
    for (auto const& s : *crv_.wrapper->product()) { if (wanted(s.simParticle())) crv->push_back(s); }
    auto surf = std::make_unique<SurfaceStepCollection>();
    for (auto const& s : *surf_.wrapper->product()) { if (wanted(s.simParticle())) surf->push_back(s); }
    auto vd = std::make_unique<StepPointMCCollection>();
    for (auto const& s : *vd_.wrapper->product()) { if (wanted(s.simParticle())) vd->push_back(s); }
    auto traj = std::make_unique<MCTrajectoryCollection>();
    for (auto const& kv : *traj_.wrapper->product()) { if (wanted(kv.first)) traj->insert(kv); }

    art::put_product_in_principal(std::move(sgs), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::move(css), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::move(crv), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::move(surf), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::move(vd), *outE, "compressDetStepMCs", "virtualdetector");
    art::put_product_in_principal(std::move(traj), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::make_unique<SimParticleCollection>(sims), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::make_unique<GenParticleCollection>(*gen_.wrapper->product()), *outE, "compressDetStepMCs");
    art::put_product_in_principal(std::make_unique<StatusG4>(*status_.wrapper->product()), *outE, "g4run");
    art::put_product_in_principal(std::make_unique<StepPointMCCollection>(*crvStage1_.wrapper->product()), *outE,
                                  "CosmicResampler", "crvStage1");
    art::put_product_in_principal(std::make_unique<std::vector<art::EventID>>(*eids_.wrapper->product()), *outE,
                                  "CosmicResampler");
    art::put_product_in_principal(std::make_unique<CosmicLivetime>(*evLive_.wrapper->product()), *outE,
                                  "CosmicResampler");
    art::put_product_in_principal(std::make_unique<PrimaryParticle>(*primary_.wrapper->product()), *outE,
                                  "FindMCPrimary");
    const double w = baseWeight_ * (j_ == 0 ? stripWeight_ : 1. / forced_[j_ - 1].second);
    art::put_product_in_principal(std::make_unique<EventWeight>(w), *outE, weightLabel_);

    if (verbosity_ > 1) {
      std::cout << "FamilySplitSource: input " << aux.id() << " output " << ev << " j " << j_ << " of "
                << forced_.size() << " weight " << w << std::endl;
    }
    ++nOut_;
    if (j_ > 0) { ++nFamilies_; }
    ++j_;
    return true;
  }

} // namespace mu2e

DEFINE_ART_INPUT_SOURCE(art::Source<mu2e::FamilySplitSourceDetail>)
