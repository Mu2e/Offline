// FamilySplitRules: which particles of a cosmic stage-2 event made with Geant4 secondary splitting are forced
// (split) and with which factor, and which events the muon-decay rule keeps. Shared by FamilySplitSource, which
// writes one event per forced family, and FamilySplitCounter, which counts the same families in the stage-2 job
// itself, before an event selection; see FamilySplitSource_source.cc for the rules.
// Split muon daughters are found either by where they were born (`regions`, cylinders that repeat the Geant4
// regions of the stage-2 configuration, and `processes`) or, with useTrackWeights, by their Geant4 track weight
// (SimParticle::weight), which splitting divides by the factor and descendants inherit.
#ifndef Sources_inc_FamilySplitRules_hh
#define Sources_inc_FamilySplitRules_hh

#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "canvas/Persistency/Provenance/ProductID.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e {

  class FamilySplitRules {
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
    struct CopiedPrimaries {
      fhicl::Sequence<int> pdgIds{fhicl::Name("pdgIds"), fhicl::Comment("as in StepPointMCCopier")};
      fhicl::Atom<double> minMomentum{fhicl::Name("minMomentum"), fhicl::Comment("[MeV/c], as in StepPointMCCopier")};
      fhicl::Atom<unsigned> factor{fhicl::Name("factor"), fhicl::Comment("StepPointMCCopier copies, >= 2")};
      fhicl::Atom<std::string> copierLabel{fhicl::Name("copierLabel"),
          fhicl::Comment("module label of StepPointMCCopier: its EventWeight is 1/factor (copied) or 1 (not)")};
    };
    struct Config {
      fhicl::Atom<bool> useTrackWeights{fhicl::Name("useTrackWeights"),
          fhicl::Comment("true: a stage-2 particle whose Geant4 track weight (SimParticle::weight) is below its parent's is\n"
                         "forced, factor = parent weight / weight (an integer >= 2; the parent must have weight 1): a decay\n"
                         "family if it is a decay product of a muon that decayed in flight, a region family otherwise.\n"
                         "processes and regions must then be left out. false: by processes and regions (files written\n"
                         "before SimParticle had a weight)"), false};
      fhicl::OptionalSequence<std::string> processes{fhicl::Name("processes"),
          fhicl::Comment("without useTrackWeights: creation codes of forced daughters of muons, split by region")};
      fhicl::OptionalSequence<fhicl::Table<Region>> regions{fhicl::Name("regions"),
          fhicl::Comment("without useTrackWeights: where the daughters of muons were split, and the factors")};
      fhicl::Atom<unsigned> decayFactor{fhicl::Name("decayFactor"),
          fhicl::Comment("0: muon decays in flight are not copied. N >= 2: Mu2eG4 processCopies made every muon decay\n"
                         "in flight N times; their direct daughters are forced, factor N, wherever the muon decays.\n"
                         "With useTrackWeights the weights give the factor; N >= 2 is then checked against them"), 0};
      fhicl::OptionalTable<CopiedPrimaries> copiedPrimaries{fhicl::Name("copiedPrimaries"),
          fhicl::Comment("stage-2 primaries made from StepPointMCCopier copies; absent: none")};
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
    };

    enum Kind { kRegion = 0, kDecay = 1, kPrimary = 2, nKinds = 3 };
    struct Forced { size_t key; double factor; Kind kind; };

    explicit FamilySplitRules(Config const& conf)
      : useWeights_(conf.useTrackWeights())
      , decayFactor_(conf.decayFactor())
      , window_(conf.muonDecayWindow())
      , lifetime_(conf.muonLifetime())
      , forcedStream_(conf.forcedDecayStream())
      , primaryStage_(conf.primaryStage()) {
      if (window_ < 0. || !(lifetime_ > 0.) || (forcedStream_ && window_ == 0.)) {
        throw cet::exception("CONFIG") << "FamilySplitRules: need muonDecayWindow >= 0, muonLifetime > 0, and "
                                       << "muonDecayWindow > 0 with forcedDecayStream\n";
      }
      if (forcedStream_) { baseWeight_ = -std::expm1(-window_ / lifetime_); }
      conf.processes(processes_);
      std::vector<Region> regions;
      conf.regions(regions);
      for (auto const& r : regions) {
        if (!(r.factor() >= 1.)) {
          throw cet::exception("CONFIG") << "FamilySplitRules: region factor must be >= 1\n";
        }
        regions_.push_back({r.xc(), r.yc(), r.rmin(), r.rmax(), r.zmin(), r.zmax(), r.factor()});
      }
      CopiedPrimaries cp;
      if (conf.copiedPrimaries(cp)) {
        for (int id : cp.pdgIds()) { copyPdgIds_.insert(id); }
        copyMinMomentum_ = cp.minMomentum();
        copyFactor_ = cp.factor();
        copierLabel_ = cp.copierLabel();
        if (copyPdgIds_.empty() || copyFactor_ < 2 || copyMinMomentum_ < 0. || copierLabel_.empty()) {
          throw cet::exception("CONFIG") << "FamilySplitRules: copiedPrimaries needs pdgIds, factor >= 2, minMomentum >= 0"
                                         << " and copierLabel\n";
        }
      }
      if (decayFactor_ == 1) {
        throw cet::exception("CONFIG") << "FamilySplitRules: decayFactor must be 0 (off) or >= 2\n";
      }
      if (useWeights_) {
        if (!regions_.empty() || !processes_.empty()) {
          throw cet::exception("CONFIG") << "FamilySplitRules: useTrackWeights finds the split particles by their "
                                         << "weights; leave out regions and processes\n";
        }
      } else {
        if (regions_.empty() != processes_.empty()) {
          throw cet::exception("CONFIG") << "FamilySplitRules: need both regions and processes (or neither)\n";
        }
        if (processes_.empty() && decayFactor_ == 0 && copyFactor_ == 0) {
          throw cet::exception("CONFIG") << "FamilySplitRules: nothing is forced\n";
        }
      }
    }

    bool useTrackWeights() const { return useWeights_; }
    std::vector<std::string> const& processes() const { return processes_; }
    unsigned decayFactor() const { return decayFactor_; }
    unsigned copyFactor() const { return copyFactor_; }          // 0: no copied primaries
    std::string const& copierLabel() const { return copierLabel_; }
    double window() const { return window_; }
    bool forcedStream() const { return forcedStream_; }
    double baseWeight() const { return baseWeight_; }

    // statistics of keepByDecay
    unsigned long nDecayDropped() const { return nDecayDropped_; }
    unsigned long nMultiMuon() const { return nMultiMuon_; }
    unsigned long nLateForced() const { return nLateForced_; }

    // the muon-decay rule (forcedDecayStream): false if the event is dropped
    bool keepByDecay(SimParticleCollection const& sims) {
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

    // copy factor StepPointMCCopier used for an event, from its EventWeight w: copyFactor(), or 1 if it did not
    // copy the event (0 without copied primaries). `what` names the weight in error messages.
    unsigned copyFactorFromWeight(double w, std::string const& what) const {
      if (copyFactor_ == 0) { return 0; }
      if (!(w > 0.)) {
        throw cet::exception("FamilySplitRules") << what << ": weight " << w << "\n";
      }
      const double f = 1. / w;
      if (std::abs(f - 1.) < 1.e-9) { return 1; }
      if (std::abs(f - static_cast<double>(copyFactor_)) < 1.e-6 * copyFactor_) { return copyFactor_; }
      throw cet::exception("FamilySplitRules") << what << ": copy factor " << f
                                               << ", expected 1 or copiedPrimaries.factor = " << copyFactor_ << "\n";
    }

    // the forced particles of an event; simsPID is the ProductID of `sims`
    void findForced(SimParticleCollection const& sims, art::ProductID const& simsPID, unsigned eventCopyFactor,
                    std::vector<Forced>& out) const {
      out.clear();
      std::map<size_t, unsigned> decayElectrons;                                   // muon key -> e+- copies
      std::map<std::tuple<size_t, double, double, double>, unsigned> copies;       // stage-1 step -> copies
      for (auto const& kv : sims) {
        SimParticle const& p = kv.second;
        const size_t key = kv.first.asUint();
        if (static_cast<int>(p.simStage()) == primaryStage_ && p.creationCode().id() == ProcessCode::mu2ePrimary) {
          if (copyFactor_ == 0 || copyPdgIds_.count(p.pdgId()) == 0) { continue; }
          const CLHEP::Hep3Vector mom = p.startMomentum().vect();
          const double pm = mom.mag();
          if (std::abs(pm - copyMinMomentum_) < 1.e-6) {
            throw cet::exception("FamilySplitRules") << "stage-2 primary " << key << " momentum " << pm
                                                     << " MeV/c within 1e-6 of copiedPrimaries.minMomentum\n";
          }
          if (pm <= copyMinMomentum_) { continue; }
          const size_t parentKey = p.hasParent() ? p.parent().key() : std::numeric_limits<size_t>::max();
          if (++copies[{parentKey, mom.x(), mom.y(), mom.z()}] > eventCopyFactor) {
            throw cet::exception("FamilySplitRules") << "a stage-1 step has more than " << eventCopyFactor
                                                     << " copies (the event's copy factor)\n";
          }
          if (eventCopyFactor > 1) { out.push_back({key, static_cast<double>(eventCopyFactor), kPrimary}); }
          continue;
        }
        if (!p.hasParent() || p.parent().id() != simsPID) { continue; }
        SimParticle const* par = sims.getOrNull(cet::map_vector_key(p.parent().key()));
        if (useWeights_) {
          if (static_cast<int>(p.simStage()) != primaryStage_) { continue; }
          Forced x{};
          if (byWeight(p, key, par, decayElectrons, x)) { out.push_back(x); }
          continue;
        }
        if (par == nullptr || std::abs(par->pdgId()) != 13) { continue; }
        const std::string code = p.creationCode().name();
        if (code == "Decay" || code == "DecayWithSpin") {
          if (decayFactor_ == 0 || !(par->endMomentum().vect().mag() > 0.)) { continue; }   // at rest: not copied
          if (std::abs(p.pdgId()) == 11 && ++decayElectrons[p.parent().key()] > decayFactor_) {
            throw cet::exception("FamilySplitRules") << "a muon decay in flight has more than decayFactor = "
                                                     << decayFactor_ << " electrons\n";
          }
          out.push_back({key, static_cast<double>(decayFactor_), kDecay});
          continue;
        }
        bool ok = false;
        for (auto const& c : processes_) { if (code == c) { ok = true; break; } }
        if (!ok) { continue; }
        const double f = regionFactor(p);
        if (f > 0.) { out.push_back({key, f, kRegion}); }
      }
    }

  private:
    struct Cyl { double xc, yc, rmin, rmax, zmin, zmax, factor; };

    // useTrackWeights: is p (key `key`, parent `par` in the same collection, or nullptr if the parent is not in it)
    // the root of a split family? Geant4 splitting (physics.secondaryBiasing, physics.processCopies) gives each of
    // N alternative secondaries 1/N of the parent's weight; their descendants inherit it, so the weight drops only
    // at the roots.
    bool byWeight(SimParticle const& p, size_t key, SimParticle const* par, std::map<size_t, unsigned>& decayElectrons,
                  Forced& out) const {
      constexpr double tol = 1.e-6;
      const double w = p.weight();
      if (!(w > 0.)) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << " has weight " << w << "\n";
      }
      if (par == nullptr) {
        if (std::abs(w - 1.) > tol) {
          throw cet::exception("FamilySplitRules") << "SimParticle " << key << " has weight " << w
                                                   << " but its parent is not in the collection\n";
        }
        return false;
      }
      const double wp = par->weight();
      if (!(wp > 0.) || w > wp * (1. + tol)) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << " has weight " << w << ", its parent "
                                                 << p.parent().key() << " " << wp
                                                 << ": a weight above the parent's (Russian roulette) is not handled\n";
      }
      if (w >= wp * (1. - tol)) { return false; }
      if (std::abs(wp - 1.) > tol) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << " (weight " << w << ") is split from "
                                                 << p.parent().key() << ", whose weight is " << wp
                                                 << ", not 1: splitting inside a split family is not handled\n";
      }
      const double f = wp / w;
      const double n = std::round(f);
      if (n < 2. || std::abs(f - n) > tol * n) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << ": weight " << w << ", parent weight " << wp
                                                 << ", split factor " << f << " is not an integer >= 2\n";
      }
      const std::string code = p.creationCode().name();
      Kind kind = kRegion;
      if ((code == "Decay" || code == "DecayWithSpin") && std::abs(par->pdgId()) == 13 &&
          par->endMomentum().vect().mag() > 0.) {
        kind = kDecay;
        if (decayFactor_ > 0 && n != static_cast<double>(decayFactor_)) {
          throw cet::exception("FamilySplitRules") << "SimParticle " << key << ": a muon decay in flight copied "
                                                   << n << " times, decayFactor = " << decayFactor_ << "\n";
        }
        if (std::abs(p.pdgId()) == 11 && ++decayElectrons[p.parent().key()] > n) {
          throw cet::exception("FamilySplitRules") << "a muon decay in flight has more than its " << n
                                                   << " copies of electrons\n";
        }
      }
      out = {key, n, kind};
      return true;
    }

    double regionFactor(SimParticle const& p) const {
      auto const x = p.startPosition();
      for (auto const& c : regions_) {
        const double dx = x.x() - c.xc, dy = x.y() - c.yc;
        const double r = std::sqrt(dx * dx + dy * dy);
        if (r >= c.rmin && r < c.rmax && x.z() >= c.zmin && x.z() < c.zmax) { return c.factor; }
      }
      return 0.;
    }

    bool useWeights_ = false;
    std::vector<std::string> processes_;
    std::vector<Cyl> regions_;
    unsigned decayFactor_ = 0;
    std::set<int> copyPdgIds_;
    double copyMinMomentum_ = 0.;
    unsigned copyFactor_ = 0;
    std::string copierLabel_;
    double window_ = 0., lifetime_ = 0., baseWeight_ = 1.;
    bool forcedStream_ = false;
    int primaryStage_ = 1;
    unsigned long nDecayDropped_ = 0, nMultiMuon_ = 0, nLateForced_ = 0;
  };

} // namespace mu2e

#endif /* Sources_inc_FamilySplitRules_hh */
