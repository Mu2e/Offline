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
// A forced particle is one of which Geant4 made `factor` alternatives, each of weight 1/factor; its family is
// itself and all its descendants. Three kinds:
//   region:  a direct daughter of a muon, created by one of `processes`, born inside one of `regions`
//            (cylinders along z, Mu2e coordinates), factor of that region. Made by Mu2eG4
//            physics.secondaryBiasing (muIoni, muBrems) or physics.processCopies (muPairProd) in the
//            matching G4Regions.
//   decay:   a direct daughter of a muon that decayed in flight (creation code Decay or DecayWithSpin, the
//            muon's end momentum > 0), anywhere, factor `decayFactor`. Made by Mu2eG4 physics.processCopies
//            for the Decay process of mu+ and mu- everywhere; decays at rest are not copied.
//   primary: a stage-2 primary with one of `copiedPrimaries.pdgIds` and start momentum above
//            `copiedPrimaries.minMomentum`, factor `copiedPrimaries.factor`. Made by StepPointMCCopier, which
//            writes such stage-1 steps N times before Mu2eG4 turns them into primaries. An event with too many
//            such steps is not copied; StepPointMCCopier's EventWeight (`copiedPrimaries.copierLabel`, 1/N or 1)
//            says which, and in an event that was not copied these primaries are not forced.
// The factors here must be those used in stage 2; the job checks only that no muon decay and no stage-1 step
// has more copies than the factor.
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
// Counts from stage 2 (optional, countsLabel = the FamilySplitCounter module label of the stage-2 job). S and E
// above then come from the FamilySplitSummary that FamilySplitCounter wrote for each SubRun of the file, summed,
// instead of from the events found in the file. FamilySplitCounter runs in the stage-2 job before an event
// selection (S2ResamplerSplitSignal.fcl), so stage 2 can drop the events that no output event of this source
// would pass, and the stripped events keep the weight they have without the selection. The source checks that
// the events in the file give at most those counts (equal when stage 2 dropped nothing). Without countsLabel, a
// FamilySplitSummary in the file must equal the file's own counts: a file that stage 2 thinned stops the job.
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
#include <tuple>
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
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/TableFragment.h"
#include "art/Framework/Core/FileBlock.h"

#include "Offline/DataProducts/inc/StageNormalization.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/MCDataProducts/inc/FamilySplitSummary.hh"
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/Sources/inc/FamilySplitRules.hh"

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
    struct Config {
      fhicl::TableFragment<FamilySplitRules::Config> rules;
      fhicl::Atom<unsigned> stride{fhicl::Name("stride"), 1000};
      fhicl::Atom<std::string> weightLabel{fhicl::Name("weightLabel"), "FamilySplit"};
      fhicl::Atom<std::string> inputProcessName{fhicl::Name("inputProcessName"),
          fhicl::Comment("process name of the stage-2 products; must equal this job's process_name"), "CosmicS2"};
      fhicl::Atom<std::string> countsLabel{fhicl::Name("countsLabel"),
          fhicl::Comment("module label of the stage-2 FamilySplitCounter whose FamilySplitSummary (in every SubRun)\n"
                         "gives the stripped-event weight; empty: the events of the file give it"), ""};
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
    using Forced = FamilySplitRules::Forced;
    art::SourceHelper const& sh_;
    FamilySplitRules rules_;
    unsigned stride_;
    std::string weightLabel_, process_, countsLabel_;
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
    Branch<EventWeight> copierW_{"", ""};                         // with copiedPrimaries only
    Branch<CosmicLivetime> srLive_{"CosmicResampler", "mixed"};
    Branch<GenEventCount> srCount_{"genCounter", ""};                // optional
    Branch<StageNormalization> srNorm_{"s2Counter", ""};             // optional
    Branch<FamilySplitSummary> srSummary_{"", ""};                   // stage-2 counts: read, not passed on
    std::string triggerResultsName_;
    struct SubRunProducts {
      art::Timestamp begin;
      CosmicLivetime livetime;
      std::optional<GenEventCount> count;
      std::optional<StageNormalization> norm;
      std::optional<FamilySplitSummary> summary;
      bool written = false;
    };
    std::map<art::SubRunID, SubRunProducts> subRunProducts_;

    // state
    Long64_t nEntries_ = 0, entry_ = -1;
    std::vector<char> keep_;                                  // per input entry, after the muon-decay rule
    std::vector<Forced> forced_;                              // the forced particles of the current input event
    unsigned eventCopyFactor_ = 0;                            // copy factor of the current input event (1: not copied)
    std::unordered_set<size_t> roots_;                        // the forced keys of the current input event
    std::unordered_map<size_t, long> family_;                 // particle key -> forced key or -1
    size_t j_ = 1;                                            // next output index of the current input event
    double stripWeight_ = 1.;
    art::RunNumber_t lastRun_ = 0;
    art::SubRunID lastSubRun_;
    bool checkedProcess_ = false;

    // totals
    unsigned long nIn_ = 0, nOut_ = 0, nFamilies_ = 0;
    unsigned long nEmptySubRuns_ = 0;
    unsigned long nNotCopied_ = 0;                            // input events StepPointMCCopier did not copy
    std::map<size_t, unsigned long> multiplicity_;
    unsigned long nFamiliesByKind_[FamilySplitRules::nKinds] = {0, 0, 0};

    void writeSubRun(art::SubRunID const& id, art::Timestamp const& ts,
                     art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR);
    unsigned copyFactorOf(Long64_t entry);
    long familyOf(art::Ptr<SimParticle> const& ptr, SimParticleCollection const& sims);
    void checkBranches();
  };

  FamilySplitSourceDetail::FamilySplitSourceDetail(Parameters const& conf, art::ProductRegistryHelper& rh,
                                                   art::SourceHelper const& sh)
    : sh_(sh)
    , rules_(conf().rules())
    , stride_(conf().stride())
    , weightLabel_(conf().weightLabel())
    , process_(conf().inputProcessName())
    , countsLabel_(conf().countsLabel())
    , verbosity_(conf().verbosity()) {
    if (rules_.copyFactor() > 0) { copierW_ = Branch<EventWeight>(rules_.copierLabel(), ""); }
    if (stride_ < 2) {
      throw cet::exception("CONFIG") << "FamilySplitSource: need stride >= 2\n";
    }
    if (!countsLabel_.empty()) {
      srSummary_ = Branch<FamilySplitSummary>(countsLabel_, "");
      srSummary_.setName(process_);
    }
    simsPID_ = art::ProductID{art::canonicalProductName(art::TypeID(typeid(SimParticleCollection)).friendlyClassName(),
                                                        "compressDetStepMCs", "", process_)};
    sims_.setName(process_); sgs_.setName(process_); css_.setName(process_); crv_.setName(process_);
    surf_.setName(process_); vd_.setName(process_); traj_.setName(process_); gen_.setName(process_);
    status_.setName(process_); crvStage1_.setName(process_); eids_.setName(process_); evLive_.setName(process_);
    primary_.setName(process_); srLive_.setName(process_); srCount_.setName(process_); copierW_.setName(process_);
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
    std::cout << "FamilySplitSource: family events by kind: region " << nFamiliesByKind_[FamilySplitRules::kRegion]
              << " (factors by region, processes";
    for (auto const& p : rules_.processes()) { std::cout << " " << p; }
    std::cout << "), decay " << nFamiliesByKind_[FamilySplitRules::kDecay] << " (factor " << rules_.decayFactor()
              << "), stage-1 copies " << nFamiliesByKind_[FamilySplitRules::kPrimary] << " (factor "
              << rules_.copyFactor() << "; input events not copied " << nNotCopied_ << ")" << std::endl;
    if (rules_.window() > 0.) {
      std::cout << "FamilySplitSource: " << (rules_.forcedStream() ? "forced-decay" : "main")
                << " stream, muon decay window " << rules_.window() << " ns: " << rules_.nDecayDropped()
                << " input events dropped, base weight " << rules_.baseWeight()
                << ", events with several primary muons (approximate) " << rules_.nMultiMuon();
      if (rules_.forcedStream()) { std::cout << ", primary muon decays later than the window " << rules_.nLateForced(); }
      std::cout << std::endl;
    }
    std::cout << "FamilySplitSource: subruns written without events " << nEmptySubRuns_ << std::endl;
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

  // copy factor StepPointMCCopier used for this input entry: copyFactor, or 1 if it did not copy the event
  unsigned FamilySplitSourceDetail::copyFactorOf(Long64_t entry) {
    if (rules_.copyFactor() == 0) { return 0; }
    return rules_.copyFactorFromWeight(copierW_.read(entry).weight(), copierW_.name + " entry " + std::to_string(entry));
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
    if (rules_.copyFactor() > 0) { known.insert(copierW_.name); }   // read here, not passed on
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
    // stage-2 counts (FamilySplitSummary, read here, not passed on): with countsLabel the named one, required;
    // without it, any one the file has, which must then agree with the file's own counts (stage 2 dropped nothing)
    if (countsLabel_.empty()) {
      srSummary_.name.clear();
      const std::string prefix = "mu2e::FamilySplitSummary_", suffix = "__" + process_ + ".";
      for (auto* obj : *subRuns->GetListOfBranches()) {
        const std::string n = obj->GetName();
        if (n.rfind(prefix, 0) == 0 && n.size() > suffix.size() &&
            n.compare(n.size() - suffix.size(), suffix.size(), suffix) == 0) {
          if (!srSummary_.name.empty()) {
            throw cet::exception("FamilySplitSource") << filename << ": several FamilySplitSummary branches ("
                                                      << srSummary_.name << ", " << n << "); set countsLabel\n";
          }
          srSummary_.name = n;
        }
      }
    }
    const bool hasSummary = !srSummary_.name.empty();
    {
      std::set<std::string> known{"SubRunAuxiliary", srLive_.name, srCount_.name, srNorm_.name};
      if (hasSummary) { known.insert(srSummary_.name); }
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
    if (rules_.copyFactor() > 0) { copierW_.attach(events_); }

    // subrun products
    subRunProducts_.clear();
    art::SubRunAuxiliary* srAux = nullptr;
    subRuns->GetBranch("SubRunAuxiliary")->SetAddress(&srAux);
    srLive_.attach(subRuns);
    const bool hasCount = srCount_.tryAttach(subRuns);
    const bool hasNorm = srNorm_.tryAttach(subRuns);
    if (hasSummary) { srSummary_.attach(subRuns); }                  // then required in every SubRun
    for (Long64_t i = 0; i < subRuns->GetEntries(); ++i) {
      subRuns->GetBranch("SubRunAuxiliary")->GetEntry(i);
      SubRunProducts sp{srAux->beginTime(), srLive_.read(i), std::nullopt, std::nullopt};
      if (hasCount) { sp.count = srCount_.readIfPresent(i); }
      if (hasNorm) { sp.norm = srNorm_.readIfPresent(i); }
      if (hasSummary) {
        sp.summary = srSummary_.readIfPresent(i);
        if (!sp.summary) {
          throw cet::exception("FamilySplitSource") << filename << ": " << srAux->id() << " has no " << srSummary_.name
                                                    << "\n";
        }
      }
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
    std::vector<Forced> f;
    for (Long64_t i = 0; i < nEntries_; ++i) {
      SimParticleCollection const& sims = sims_.read(i);
      if (!rules_.keepByDecay(sims)) { continue; }
      keep_[i] = 1;
      ++nKept;
      rules_.findForced(sims, simsPID_, copyFactorOf(i), f);
      for (auto const& x : f) { sumInv += 1. / x.factor; }
    }
    stripWeight_ = nKept > 0 ? 1. - sumInv / static_cast<double>(nKept) : 1.;
    std::cout << "FamilySplitSource: " << filename << ": " << nEntries_ << " events, " << nKept
              << " kept, sum of family weights " << sumInv << ", stripped-event weight " << stripWeight_
              << ", base weight " << rules_.baseWeight() << std::endl;
    if (hasSummary) {
      // the stage-2 counts, taken before its event selection: the file's events can only give less
      unsigned long nKeptS2 = 0, nEventsS2 = 0;
      double sumInvS2 = 0.;
      for (auto const& [id, sp] : subRunProducts_) {
        nEventsS2 += sp.summary->nEvents();
        nKeptS2 += sp.summary->nKept();
        sumInvS2 += sp.summary->sumWeights();
      }
      if (static_cast<unsigned long>(nEntries_) > nEventsS2 || static_cast<unsigned long>(nKept) > nKeptS2 ||
          sumInv > sumInvS2 + 1.e-9 * std::max(1., sumInvS2)) {
        throw cet::exception("FamilySplitSource") << filename << ": the file has more than its stage-2 counts "
                                                  << srSummary_.name << ": events " << nEntries_ << " > " << nEventsS2
                                                  << " or kept " << nKept << " > " << nKeptS2 << " or family weights "
                                                  << sumInv << " > " << sumInvS2 << "\n";
      }
      const bool same = static_cast<unsigned long>(nEntries_) == nEventsS2 &&
                        static_cast<unsigned long>(nKept) == nKeptS2 &&
                        std::abs(sumInv - sumInvS2) <= 1.e-9 * std::max(1., sumInvS2);
      if (countsLabel_.empty() && !same) {
        throw cet::exception("FamilySplitSource") << filename << ": stage 2 dropped events after counting them ("
                                                  << srSummary_.name << ": " << nEventsS2 << " events, " << nKeptS2
                                                  << " kept, family weights " << sumInvS2 << "); set countsLabel\n";
      }
      if (!countsLabel_.empty()) { stripWeight_ = nKeptS2 > 0 ? 1. - sumInvS2 / static_cast<double>(nKeptS2) : 1.; }
      std::cout << "FamilySplitSource: " << filename << ": stage-2 counts (" << srSummary_.name << ") " << nEventsS2
                << " events, " << nKeptS2 << " kept, sum of family weights " << sumInvS2
                << (same ? ", the same as the file's" : "") << ", stripped-event weight " << stripWeight_ << std::endl;
    }
    entry_ = -1;
    j_ = 1;
    forced_.clear();
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
      eventCopyFactor_ = copyFactorOf(entry_);
      if (eventCopyFactor_ == 1) { ++nNotCopied_; }
      rules_.findForced(*sims_.wrapper->product(), simsPID_, eventCopyFactor_, forced_);
      if (forced_.size() + 1 > stride_) {
        throw cet::exception("FamilySplitSource") << forced_.size() << " families in one event, stride "
                                                  << stride_ << " too small\n";
      }
      roots_.clear();
      for (auto const& x : forced_) { roots_.insert(x.key); }
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
    const long kept = j_ == 0 ? -2 : static_cast<long>(forced_[j_ - 1].key);
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
    const double w = rules_.baseWeight() * (j_ == 0 ? stripWeight_ : 1. / forced_[j_ - 1].factor);
    art::put_product_in_principal(std::make_unique<EventWeight>(w), *outE, weightLabel_);

    if (verbosity_ > 1) {
      std::cout << "FamilySplitSource: input " << aux.id() << " output " << ev << " j " << j_ << " of "
                << forced_.size() << " weight " << w << std::endl;
    }
    ++nOut_;
    if (j_ > 0) { ++nFamilies_; ++nFamiliesByKind_[forced_[j_ - 1].kind]; }
    ++j_;
    return true;
  }

} // namespace mu2e

DEFINE_ART_INPUT_SOURCE(art::Source<mu2e::FamilySplitSourceDetail>)
