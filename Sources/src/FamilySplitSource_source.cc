// FamilySplitSource: one pass over cosmic stage-2 dts files made with Geant4 secondary splitting
// (Mu2eG4 physics.processCopies), writing for every input event
//   j = 0      the event with every split ("forced") family removed (the "stripped" event), and
//   j = 1..k   one event per forced family, holding the rest of the event plus that family only.
// Output event number = input event number * stride + j.  Configuration:
// Production/JobConfig/cosmic/S2FamilySplitMain.fcl.
//
// Why: splitting makes N alternative secondaries at one interaction, all in one event. Downstream
// stages count events, so two alternatives that both reach the tracker would be counted once and
// the event would hold two signal-like electrons at once. One family per event removes both problems
// and lets N be large.
//
// Which particles are forced (split), with which factor, and the muon-decay rule (muonDecayWindow,
// forcedDecayStream): FamilySplitRules.hh. Split particles are found by their Geant4 track weight, so the job stops
// on a file written before SimParticle had a weight.
//
// Every output event carries an EventWeight (label `weightLabel`), times the base weight (1, or the
// forced-decay weight below):
//   family event:   the weight of its forced particle (1/N);
//   stripped event: 1 - (sum over the file's forced families of their weights) / (events kept from the file).
// The weights add up to the number of events kept, so the dataset keeps the livetime of its input
// (CosmicLivetime is copied unchanged). Anything a family carries (the particles it makes) keeps its
// expectation exactly. What the rest of an input event holds (the muon, unsplit particles) is written in
// the stripped event and again in each of its family events, with total weight 1 - S/E + s, where s is
// the sum of that event's family weights: right on average over the file, biased for quantities that
// correlate with s. (The exact per-event stripped weight, 1 - s, can be negative.)
//
// SubRuns are all written, also those left without events, so no livetime is lost.
//
// Counts from stage 2 (optional, countsLabel = the FamilySplitCounter module label of the stage-2 job). S and E
// above then come from the FamilySplitSummary that FamilySplitCounter wrote for each SubRun of the file, summed,
// instead of from the events found in the file. FamilySplitCounter runs in the stage-2 job before an event
// selection (S2ResamplerSplitMain.fcl), so stage 2 can drop the events that no output event of this source
// would pass, and the stripped events keep the weight they have without the selection. The source checks that
// the events in the file give at most those counts (equal when stage 2 dropped nothing). Without countsLabel, a
// FamilySplitSummary in the file must equal the file's own counts: a file that stage 2 thinned stops the job.
//
// Products are re-created with their stage-2 names, so the job's process_name must be the stage-2 one
// (`inputProcessName`): art derives a ProductID from the product name, so the art::Ptr in the copied
// products still point to the right collections. The source gives the Ptr into the SimParticles and GenParticles
// of the event their getter, so that modules of this job can follow them. Detector step collections
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
#include <utility>
#include <vector>

#include "TBranch.h"
#include "TFile.h"
#include "TList.h"
#include "TObjArray.h"
#include "TStreamerInfo.h"
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
    // One art product branch of the input Events or SubRuns tree, by module label and instance name.
    template <typename T>
    struct Branch {
      std::string label, instance, name;
      TBranch* branch = nullptr;
      art::Wrapper<T>* wrapper = nullptr;
      Branch(std::string l, std::string i) : label(std::move(l)), instance(std::move(i)) {}
      // the branch name is the canonical product name (it ends in '.')
      void setName(std::string const& process) {
        name = art::canonicalProductName(art::TypeID(typeid(T)).friendlyClassName(), label, instance, process);
      }
      bool tryAttach(TTree* tree) {
        branch = tree->GetBranch(name.c_str());
        wrapper = nullptr;
        if (branch != nullptr) { branch->SetAddress(&wrapper); }
        return branch != nullptr;
      }
      void attach(TTree* tree) {
        if (!tryAttach(tree)) {
          throw cet::exception("FamilySplitSource") << "no branch " << name << " in tree " << tree->GetName() << "\n";
        }
      }
      // reads an entry; false if the product is marked not present
      bool load(Long64_t entry) {
        if (branch->GetEntry(entry) <= 0 || wrapper == nullptr) {
          throw cet::exception("FamilySplitSource") << "cannot read " << name << " entry " << entry << "\n";
        }
        return wrapper->isPresent();
      }
      T const& read(Long64_t entry) {
        if (!load(entry)) {
          throw cet::exception("FamilySplitSource") << "cannot read " << name << " entry " << entry << "\n";
        }
        return *wrapper->product();
      }
      // for optional products: the branch may exist with the product marked not present
      std::optional<T> readIfPresent(Long64_t entry) {
        if (!load(entry)) { return std::nullopt; }
        return *wrapper->product();
      }
      T const& product() const { return *wrapper->product(); }              // of the entry read last
      void reconstitute(art::ProductRegistryHelper& rh) const { rh.reconstitutes<T, art::InEvent>(label, instance); }
      void put(std::unique_ptr<T> p, art::EventPrincipal& e) const {
        art::put_product_in_principal(std::move(p), e, label, instance);
      }
      void putCopy(art::EventPrincipal& e) const { put(std::make_unique<T>(product()), e); }
    };

    // calls f on each branch of a tuple of branches
    template <typename Tuple, typename F>
    void forEach(Tuple&& branches, F&& f) {
      std::apply([&](auto&... b) { (f(b), ...); }, std::forward<Tuple>(branches));
    }

    // a step's SimParticle Ptr (CaloShowerStep has a setter only)
    void setSimParticle(CaloShowerStep& s, art::Ptr<SimParticle> const& p) { s.setSimParticle(p); }
    template <typename S>
    void setSimParticle(S& s, art::Ptr<SimParticle> const& p) { s.simParticle() = p; }

    // the steps of `in` whose SimParticle `wanted` accepts, their SimParticle Ptr passed through `repoint`
    template <typename C, typename F, typename R>
    std::unique_ptr<C> keepIf(C const& in, F const& wanted, R const& repoint) {
      auto out = std::make_unique<C>();
      for (auto const& s : in) {
        if (!wanted(s.simParticle())) { continue; }
        out->push_back(s);
        setSimParticle(out->back(), repoint(s.simParticle()));
      }
      return out;
    }

    // stops the job on a branch of `tree` that is not in `known`: it would be lost
    void requireKnownBranches(TTree* tree, std::set<std::string> const& known) {
      for (auto* obj : *tree->GetListOfBranches()) {
        if (known.count(obj->GetName()) == 0) {
          throw cet::exception("FamilySplitSource") << "input " << tree->GetName() << " branch " << obj->GetName()
                                                    << " is not handled; it would be lost\n";
        }
      }
    }
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
    art::ProductID simsPID_, genPID_;

    // input
    std::unique_ptr<TFile> file_;
    TTree* events_ = nullptr;
    TBranch* eventAuxBranch_ = nullptr;
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
    // every Events product passed on, in the order they are registered
    auto eventBranches() {
      return std::tie(sims_, sgs_, css_, crv_, surf_, vd_, traj_, gen_, status_, crvStage1_, eids_, evLive_, primary_);
    }
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
    std::map<size_t, unsigned long> multiplicity_;
    unsigned long nFamiliesByKind_[FamilySplitRules::nKinds] = {0, 0, 0};

    void writeSubRun(art::SubRunID const& id, art::Timestamp const& ts,
                     art::RunPrincipal*& outR, art::SubRunPrincipal*& outSR);
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
    if (stride_ < 2) {
      throw cet::exception("CONFIG") << "FamilySplitSource: need stride >= 2\n";
    }
    if (!countsLabel_.empty()) {
      srSummary_ = Branch<FamilySplitSummary>(countsLabel_, "");
      srSummary_.setName(process_);
    }
    forEach(eventBranches(), [&](auto& b) { b.setName(process_); });
    srLive_.setName(process_); srCount_.setName(process_); srNorm_.setName(process_);
    simsPID_ = art::ProductID{sims_.name};
    genPID_ = art::ProductID{gen_.name};
    triggerResultsName_ = art::canonicalProductName("art::TriggerResults", "TriggerResults", "", process_);

    forEach(eventBranches(), [&](auto const& b) { b.reconstitute(rh); });
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
              << ", decay " << nFamiliesByKind_[FamilySplitRules::kDecay] << ", stage-1 copies "
              << nFamiliesByKind_[FamilySplitRules::kPrimary] << std::endl;
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
    std::set<std::string> known{"EventAuxiliary", triggerResultsName_};
    forEach(eventBranches(), [&](auto const& b) { known.insert(b.name); });
    requireKnownBranches(events_, known);
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
    // a file written before SimParticle had a weight reads weight 1 everywhere, so no family would be found
    bool hasWeight = false;
    if (TList const* infos = file_->GetStreamerInfoCache()) {
      for (TObject const* obj : *infos) {
        auto const* info = dynamic_cast<TStreamerInfo const*>(obj);
        if (info != nullptr && std::string(info->GetName()) == "mu2e::SimParticle" &&
            info->GetElements()->FindObject("_weight") != nullptr) {
          hasWeight = true;
          break;
        }
      }
    }
    if (!hasWeight) {
      throw cet::exception("FamilySplitSource") << filename << ": its SimParticles have no weight (written before "
                                                << "SimParticle::weight), so the split particles cannot be found\n";
    }
    if (runs->GetListOfBranches()->GetEntries() != 1) {
      throw cet::exception("FamilySplitSource") << filename << ": Runs tree has products, not handled\n";
    }
    // stage-2 counts (FamilySplitSummary, read here, not passed on): with countsLabel the named one, required;
    // without it, any one the file has, which must then agree with the file's own counts (stage 2 dropped nothing)
    if (countsLabel_.empty()) {
      srSummary_.name.clear();
      const std::string prefix = art::TypeID(typeid(FamilySplitSummary)).friendlyClassName() + "_";
      const std::string suffix = "__" + process_ + ".";
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
    std::set<std::string> knownSubRun{"SubRunAuxiliary", srLive_.name, srCount_.name, srNorm_.name};
    if (hasSummary) { knownSubRun.insert(srSummary_.name); }
    requireKnownBranches(subRuns, knownSubRun);
    eventAux_ = nullptr;
    eventAuxBranch_ = events_->GetBranch("EventAuxiliary");
    eventAuxBranch_->SetAddress(&eventAux_);
    forEach(eventBranches(), [&](auto& b) { b.attach(events_); });

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
    FamilySplitSummary file;
    for (Long64_t i = 0; i < nEntries_; ++i) {
      keep_[i] = rules_.count(sims_.read(i), simsPID_, file);
    }
    auto stripped = [](FamilySplitSummary const& c) {
      return c.nKept() > 0 ? 1. - c.sumWeights() / static_cast<double>(c.nKept()) : 1.;
    };
    stripWeight_ = stripped(file);
    std::cout << "FamilySplitSource: " << filename << ": " << file.nEvents() << " events, " << file.nKept()
              << " kept, sum of family weights " << file.sumWeights() << ", stripped-event weight " << stripWeight_
              << ", base weight " << rules_.baseWeight() << std::endl;
    if (hasSummary) {
      // the stage-2 counts, taken before its event selection: the file's events can only give less
      FamilySplitSummary s2;
      for (auto const& [id, sp] : subRunProducts_) { s2.aggregate(*sp.summary); }
      const double tol = 1.e-9 * std::max(1., s2.sumWeights());
      if (file.nEvents() > s2.nEvents() || file.nKept() > s2.nKept() || file.sumWeights() > s2.sumWeights() + tol) {
        throw cet::exception("FamilySplitSource") << filename << ": the file has more than its stage-2 counts "
                                                  << srSummary_.name << ": events " << file.nEvents() << " > "
                                                  << s2.nEvents() << " or kept " << file.nKept() << " > " << s2.nKept()
                                                  << " or family weights " << file.sumWeights() << " > "
                                                  << s2.sumWeights() << "\n";
      }
      const bool same = file.nEvents() == s2.nEvents() && file.nKept() == s2.nKept() &&
                        std::abs(file.sumWeights() - s2.sumWeights()) <= tol;
      if (countsLabel_.empty() && !same) {
        throw cet::exception("FamilySplitSource") << filename << ": stage 2 dropped events after counting them ("
                                                  << srSummary_.name << ": " << s2.nEvents() << " events, "
                                                  << s2.nKept() << " kept, family weights " << s2.sumWeights()
                                                  << "); set countsLabel\n";
      }
      if (!countsLabel_.empty()) { stripWeight_ = stripped(s2); }
      std::cout << "FamilySplitSource: " << filename << ": stage-2 counts (" << srSummary_.name << ") "
                << s2.nEvents() << " events, " << s2.nKept() << " kept, sum of family weights " << s2.sumWeights()
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
      if (eventAuxBranch_->GetEntry(entry_) <= 0 || eventAux_ == nullptr) {
        throw cet::exception("FamilySplitSource") << "cannot read EventAuxiliary entry " << entry_ << "\n";
      }
      // read every product of this input event once; the wrappers stay valid until the next entry
      forEach(eventBranches(), [&](auto& b) { b.read(entry_); });
      rules_.findForced(sims_.product(), simsPID_, forced_);
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

    SimParticleCollection const& sims = sims_.product();
    const long kept = j_ == 0 ? -2 : static_cast<long>(forced_[j_ - 1].key);
    auto wanted = [&](art::Ptr<SimParticle> const& p) {
      const long f = familyOf(p, sims);
      return f == -1 || f == kept;
    };
    // the copied Ptr into this event's SimParticles and GenParticles get this event as their getter (Ptr into the
    // stage-1 products keep theirs: those products are not in the event)
    art::EDProductGetter const* simsGetter = outE->productGetter(simsPID_);
    art::EDProductGetter const* genGetter = outE->productGetter(genPID_);
    auto repoint = [&](art::Ptr<SimParticle> const& p) {
      return p.id() == simsPID_ ? art::Ptr<SimParticle>(simsPID_, p.key(), simsGetter) : p;
    };
    // detector steps and trajectories are filtered, everything else is copied
    forEach(std::tie(sgs_, css_, crv_, surf_, vd_),
            [&](auto const& b) { b.put(keepIf(b.product(), wanted, repoint), *outE); });
    auto traj = std::make_unique<MCTrajectoryCollection>();
    for (auto const& [sim, t] : traj_.product()) {
      if (!wanted(sim)) { continue; }
      MCTrajectory copy(t);
      copy.sim() = repoint(t.sim());
      traj->emplace_hint(traj->end(), repoint(sim), std::move(copy));
    }
    traj_.put(std::move(traj), *outE);
    auto simsOut = std::make_unique<SimParticleCollection>(sims);
    for (auto& [key, p] : *simsOut) {
      p.parent() = repoint(p.parent());
      for (auto& d : p.daughters()) { d = repoint(d); }
      if (p.genParticle().id() == genPID_) { p.genParticle() = art::Ptr<GenParticle>(genPID_, p.genParticle().key(), genGetter); }
    }
    sims_.put(std::move(simsOut), *outE);
    auto primary = std::make_unique<PrimaryParticle>(primary_.product());
    for (auto& p : primary->modifySimParticles()) { p = repoint(p); }
    primary_.put(std::move(primary), *outE);
    forEach(std::tie(gen_, status_, crvStage1_, eids_, evLive_), [&](auto const& b) { b.putCopy(*outE); });
    const double w = rules_.baseWeight() * (j_ == 0 ? stripWeight_ : forced_[j_ - 1].weight);
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
