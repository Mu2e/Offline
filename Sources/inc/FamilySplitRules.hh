// FamilySplitRules: which particles of a cosmic stage-2 event made with Geant4 secondary splitting are forced
// (split) and with which factor, and which events the muon-decay rule keeps. Shared by FamilySplitSource, which
// writes one event per forced family, and FamilySplitCounter, which counts the same families in the stage-2 job
// itself, before an event selection.
//
// A forced particle is one of which Geant4 made N alternatives, each with 1/N of its parent's weight; its family
// is itself and all its descendants, which inherit its Geant4 track weight (SimParticle::weight). So the weight
// drops only at a forced particle: one whose weight is below its parent's, by an integer factor N >= 2. The
// family's weight is the forced particle's own weight (1/N, or less for a family inside another one). Kinds:
//   region:  made by Mu2eG4 physics.processCopies (muIoni, muBrems, muPairProd of muons) in the split G4Regions;
//   decay:   a direct daughter of a muon that decayed in flight (creation code Decay or DecayWithSpin, the muon's
//            end momentum > 0), made by physics.processCopies of the muon Decay process; decays at rest are not
//            copied;
//   primary: a stage-2 primary whose weight is below 1: StepPointMCCopier writes a stage-1 e+-/gamma step N
//            times with weight 1/N, which Mu2eG4 gives the primary (inputs.primaryWeights).
// The regions and factors are set only in the stage-2 configuration. Files written before SimParticle had a
// weight read weight 1 everywhere: no family is found in them (FamilySplitSource refuses them).
//
// Muon decays in flight (optional, muonDecayWindow = T > 0, proper time in ns). They are rare, so they can come
// from a second stage-2 stream in which Mu2eG4 physics.muonMaxPreAssignedDecayProperTime = T makes every muon
// decay within T (a stage-2 muon starts at proper time 0). Only stage-2 primary muons that decay in flight count
// (end momentum > 0; for a decay at rest the SimParticle records the time the muon stopped, not the decay time,
// and its electron is below 53 MeV anyway):
//   main stream   (forcedDecayStream: false): events whose primary muon decays in flight before T are dropped;
//   forced stream (forcedDecayStream: true):  only events whose primary muon decays in flight are kept, and
//                                              every weight is multiplied by 1 - exp(-T / muonLifetime).
// The two streams then add up to the unbiased sample, each normalized to its own livetime. This is exact for one
// primary muon per event whatever happens to it (decay in flight or at rest, capture, leaving the world); events
// with several primary muons are approximate and are counted.
#ifndef Sources_inc_FamilySplitRules_hh
#define Sources_inc_FamilySplitRules_hh

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "canvas/Persistency/Provenance/ProductID.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/FamilySplitSummary.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e {

  class FamilySplitRules {
  public:
    struct Config {
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
    struct Forced { size_t key; double weight; Kind kind; };    // weight: the forced particle's track weight

    explicit FamilySplitRules(Config const& conf)
      : window_(conf.muonDecayWindow())
      , lifetime_(conf.muonLifetime())
      , forcedStream_(conf.forcedDecayStream())
      , primaryStage_(conf.primaryStage()) {
      if (window_ < 0. || !(lifetime_ > 0.) || (forcedStream_ && window_ == 0.)) {
        throw cet::exception("CONFIG") << "FamilySplitRules: need muonDecayWindow >= 0, muonLifetime > 0, and "
                                       << "muonDecayWindow > 0 with forcedDecayStream\n";
      }
      if (forcedStream_) { baseWeight_ = -std::expm1(-window_ / lifetime_); }
    }

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
        if (std::abs(p.pdgId()) != PDGCode::mu_minus || !isStage2Primary(p)) { continue; }
        ++nMu;
        if (isDecay(p.stoppingCode()) && inFlight(p)) {
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

    // adds one event to `counts`: one event, kept or not by the muon-decay rule, and if kept the weight of each of
    // its forced families. Returns whether the event is kept.
    bool count(SimParticleCollection const& sims, art::ProductID const& simsPID, FamilySplitSummary& counts) {
      const bool kept = keepByDecay(sims);
      counts.aggregate({1, kept ? 1UL : 0UL, 0.});
      if (!kept) { return false; }
      findForced(sims, simsPID, scratch_);
      for (auto const& x : scratch_) { counts.aggregate({0, 0, x.weight}); }
      return true;
    }

    // the forced particles of an event; simsPID is the ProductID of `sims`
    void findForced(SimParticleCollection const& sims, art::ProductID const& simsPID, std::vector<Forced>& out) const {
      out.clear();
      std::map<size_t, unsigned> decayElectrons;                                   // muon key -> e+- copies
      for (auto const& kv : sims) {
        SimParticle const& p = kv.second;
        if (static_cast<int>(p.simStage()) != primaryStage_) { continue; }
        // a stage-2 primary's parent is a stage-1 particle, in another collection
        SimParticle const* par = (p.hasParent() && p.parent().id() == simsPID) ?
                                 sims.getOrNull(cet::map_vector_key(p.parent().key())) : nullptr;
        Forced x{};
        if (isForced(p, kv.first.asUint(), par, decayElectrons, x)) { out.push_back(x); }
      }
    }

  private:
    static bool isDecay(ProcessCode const& c) { return c == ProcessCode::Decay || c == ProcessCode::DecayWithSpin; }
    static bool inFlight(SimParticle const& muon) { return muon.endMomentum().vect().mag() > 0.; }
    bool isStage2Primary(SimParticle const& p) const {
      return static_cast<int>(p.simStage()) == primaryStage_ && p.creationCode() == ProcessCode::mu2ePrimary;
    }
    // is p (key `key`, parent `par` in the same collection, or nullptr if the parent is not in it) a forced
    // particle (see the top of this file)?
    bool isForced(SimParticle const& p, size_t key, SimParticle const* par, std::map<size_t, unsigned>& decayElectrons,
                  Forced& out) const {
      constexpr double tol = 1.e-6;
      const double w = p.weight();
      if (!(w > 0.)) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << " has weight " << w << "\n";
      }
      const bool primary = isStage2Primary(p);
      if (par == nullptr && !primary) {
        if (std::abs(w - 1.) > tol) {
          throw cet::exception("FamilySplitRules") << "SimParticle " << key << " has weight " << w
                                                   << " but its parent is not in the collection\n";
        }
        return false;
      }
      const double wp = primary ? 1. : par->weight();                // a stage-2 primary against weight 1
      if (!(wp > 0.) || w > wp * (1. + tol)) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << " has weight " << w << ", its parent "
                                                 << p.parent().key() << " " << wp
                                                 << ": a weight above the parent's (Russian roulette) is not handled\n";
      }
      if (w >= wp * (1. - tol)) { return false; }
      const double f = wp / w;
      const double n = std::round(f);
      if (n < 2. || std::abs(f - n) > tol * n) {
        throw cet::exception("FamilySplitRules") << "SimParticle " << key << ": weight " << w << ", parent weight " << wp
                                                 << ", split factor " << f << " is not an integer >= 2\n";
      }
      Kind kind = primary ? kPrimary : kRegion;
      if (!primary && isDecay(p.creationCode()) && std::abs(par->pdgId()) == PDGCode::mu_minus && inFlight(*par)) {
        kind = kDecay;
        if (std::abs(p.pdgId()) == PDGCode::e_minus && ++decayElectrons[p.parent().key()] > n) {
          throw cet::exception("FamilySplitRules") << "a muon decay in flight copied " << n
                                                   << " times has more electrons\n";
        }
      }
      out = {key, w, kind};
      return true;
    }

    double window_ = 0., lifetime_ = 0., baseWeight_ = 1.;
    bool forcedStream_ = false;
    int primaryStage_ = 1;
    unsigned long nDecayDropped_ = 0, nMultiMuon_ = 0, nLateForced_ = 0;
    std::vector<Forced> scratch_;                               // count()'s forced particles
  };

} // namespace mu2e

#endif /* Sources_inc_FamilySplitRules_hh */
