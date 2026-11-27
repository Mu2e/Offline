// STMDetProvenance_module.cc
//
// Test module for the STMDet genealogy. One row per StepPointMC in STMDet, recording
// whether the depositing particle is a primary, truncated (parent dropped by the
// compression) or has a parent, together with two groupings of it: the top parent
// HPGeTree assigns (walking up only while the parent itself deposited in STMDet) and the
// root of its full genealogy. Nothing is aggregated, so either grouping can be rebuilt
// offline. The endJob summary prints the totals, the number of distinct groups under each
// grouping, and the birth volumes of the depositing particles.
//
// Volumes are stored as raw startVolumeIndex and resolved to names only for the endJob
// summary, and only if PhysicalVolumesTag is given.
//
// Yongyi Wu, Oct. 2026

#include <map>
#include <set>
#include <string>
#include <vector>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "TTree.h"

#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/STMMC/inc/STMGenealogy.hh"

namespace mu2e {

  class STMDetProvenance : public art::EDAnalyzer {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      struct Config {
        fhicl::Atom<art::InputTag> stepPointMCsTag{
          Name("StepPointMCsTag"),
          Comment("STMDet StepPointMCs, e.g. \"compressSTMDet:STMDet\"")};
        fhicl::Atom<art::InputTag> simParticlemvTag{
          Name("SimParticlemvTag"),
          Comment("The matching SimParticle collection, e.g. \"compressSTMDet\"")};
        // Optional: affects only how the endJob summary PRINTS volumes. Every number the
        // module produces is computed from the index, so omitting this costs no information.
        fhicl::OptionalAtom<art::InputTag> physicalVolumesTag{
          Name("PhysicalVolumesTag"),
          Comment("PhysicalVolumeInfoMultiCollection (SubRun). Used only to print volume "
                  "NAMES in the endJob summary; omit it and the summary lists the raw "
                  "indices instead. The tree always stores indices either way.")};
        fhicl::Atom<double> xBeamCentre{
          Name("xBeamCentre"),
          Comment("The x split HPGeTree uses to separate HPGe (x<) from LaBr (x>). Recorded "
                  "per step so the cut can be audited rather than applied here.")};
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit STMDetProvenance(const Parameters& conf);
      void beginSubRun(const art::SubRun& sr) override;
      void analyze(const art::Event& event) override;
      void endJob() override;

    private:
      // "<name> (idx N)" when the table resolves the index, otherwise "idx N". Used only by
      // endJob; nothing in the hot path calls it.
      std::string volumeLabel(unsigned index) const;

      art::ProductToken<StepPointMCCollection> stepsToken_;
      art::ProductToken<SimParticleCollection> simsToken_;
      art::InputTag volumesTag_;   // empty label = not supplied
      double xBeamCentre_ = -3904.0;

      // Flattened volume table, index -> name; stays empty when no tag was supplied.
      std::vector<std::string> volumeNames_;

      TTree* ttree = nullptr;

      // ---- identity of the step ----
      int    runId_ = 0, subRunId_ = 0, eventId_ = 0;  // so rows group by event offline
      double edep_ = 0.0;         // ionizingEdep of THIS step [MeV]
      double stepTime_ = 0.0;     // [ns]
      double x_ = 0.0, y_ = 0.0, z_ = 0.0;
      bool   isHPGeSide_ = false; // x < xBeamCentre, i.e. HPGeTree would keep it for HPGe

      // ---- the depositing particle ----
      int      trackId_ = 0;      // SimParticle key of the depositing particle
      int      pdgId_ = 0;
      int      creationCode_ = 0; // ProcessCode that made it; mu2ePrimary for a primary
      unsigned startVolumeIndex_ = 0;  // where it was BORN
      double   startGlobalTime_ = 0.0;

      // ---- the three classifications HPGeTree collapses into one ----
      bool   isPrimary_ = false;    // genParticle set AND creationCode == mu2ePrimary
      bool   isTruncated_ = false;  // parent Ptr null but NOT a primary: parent was dropped
      bool   hasParent_ = false;    // parent Ptr resolvable

      // ---- what HPGeTree would decide, recomputed here for direct comparison ----
      int    hpgeTopParentId_ = 0;  // the id HPGeTree's topParent() would return
      int    hpgeTopParentPdg_ = 0;
      bool   hpgeStoppedAtShell_ = false;  // a resolvable parent existed but had no STMDet
                                           // step, so the walk stopped early

      // ---- the full genealogy, walked WITHOUT the STMDet restriction ----
      // This is the grouping that ties a shell secondary back to the particle that entered
      // the detector.
      int      rootId_ = 0;       // topmost ancestor reachable in the collection
      int      rootPdg_ = 0;
      double   rootStartTime_ = 0.0;
      unsigned rootStartVolumeIndex_ = 0;
      int      genDepth_ = 0;     // generations from this particle up to rootId_
      bool     rootIsPrimary_ = false;
      bool     rootIsTruncated_ = false;

      // ---- endJob tallies ----
      unsigned long nSteps_ = 0, nPrimary_ = 0, nTruncated_ = 0, nWithParent_ = 0;
      unsigned long nStoppedAtShell_ = 0;
      std::map<int, unsigned long> rootPdgCounts_;
      // Where the depositing particles were born, and where the ones that made HPGeTree
      // stop early were born.
      std::map<unsigned, unsigned long> startVolumeCounts_, shellStopVolumeCounts_;
      // Distinct HPGeTree top parents and genealogy roots, summed over events.
      unsigned long nEventsWithSteps_ = 0, nHPGeTops_ = 0, nRoots_ = 0;
  };

  STMDetProvenance::STMDetProvenance(const Parameters& conf) :
    art::EDAnalyzer(conf),
    stepsToken_(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
    simsToken_ (consumes<SimParticleCollection>(conf().simParticlemvTag())),
    xBeamCentre_(conf().xBeamCentre())
  {
    // Absent or blank leaves volumesTag_ with an empty label, which is the "not supplied"
    // state beginSubRun and volumeLabel both test.
    art::InputTag vt;
    if (conf().physicalVolumesTag(vt)) volumesTag_ = vt;

    art::ServiceHandle<art::TFileService> tfs;
    ttree = tfs->make<TTree>("ttree", "STMDet step provenance");

    ttree->Branch("run",    &runId_,    "run/I");
    ttree->Branch("subRun", &subRunId_, "subRun/I");
    ttree->Branch("event",  &eventId_,  "event/I");

    ttree->Branch("edep",     &edep_,     "edep/D");
    ttree->Branch("stepTime", &stepTime_, "stepTime/D");
    ttree->Branch("x", &x_, "x/D");
    ttree->Branch("y", &y_, "y/D");
    ttree->Branch("z", &z_, "z/D");
    ttree->Branch("isHPGeSide", &isHPGeSide_, "isHPGeSide/O");

    ttree->Branch("trackId",          &trackId_,          "trackId/I");
    ttree->Branch("pdgId",            &pdgId_,            "pdgId/I");
    ttree->Branch("creationCode",     &creationCode_,     "creationCode/I");
    ttree->Branch("startVolumeIndex", &startVolumeIndex_, "startVolumeIndex/i");
    ttree->Branch("startGlobalTime",  &startGlobalTime_,  "startGlobalTime/D");

    ttree->Branch("isPrimary",   &isPrimary_,   "isPrimary/O");
    ttree->Branch("isTruncated", &isTruncated_, "isTruncated/O");
    ttree->Branch("hasParent",   &hasParent_,   "hasParent/O");

    ttree->Branch("hpgeTopParentId",    &hpgeTopParentId_,    "hpgeTopParentId/I");
    ttree->Branch("hpgeTopParentPdg",   &hpgeTopParentPdg_,   "hpgeTopParentPdg/I");
    ttree->Branch("hpgeStoppedAtShell", &hpgeStoppedAtShell_, "hpgeStoppedAtShell/O");

    ttree->Branch("rootId",                &rootId_,                "rootId/I");
    ttree->Branch("rootPdg",               &rootPdg_,               "rootPdg/I");
    ttree->Branch("rootStartTime",         &rootStartTime_,         "rootStartTime/D");
    ttree->Branch("rootStartVolumeIndex",  &rootStartVolumeIndex_,  "rootStartVolumeIndex/i");
    ttree->Branch("genDepth",              &genDepth_,              "genDepth/I");
    ttree->Branch("rootIsPrimary",         &rootIsPrimary_,         "rootIsPrimary/O");
    ttree->Branch("rootIsTruncated",       &rootIsTruncated_,       "rootIsTruncated/O");
  }

  void STMDetProvenance::beginSubRun(const art::SubRun& sr) {
    volumeNames_.clear();
    if (volumesTag_.label().empty()) return;

    art::Handle<PhysicalVolumeInfoMultiCollection> volsHandle;
    if (!sr.getByLabel(volumesTag_, volsHandle)) {
      mf::LogWarning("STMDetProvenance")
        << "PhysicalVolumesTag " << volumesTag_ << " not found in this SubRun; the endJob "
        << "summary will list raw volume indices. No other output is affected.";
      return;
    }

    // The multi-collection is a vector of per-stage tables, indexed by SimParticle::
    // simStage(); each table is a cet::map_vector keyed by volume index. This job is
    // single-stage, so every stage is flattened into one index->name map rather than
    // carrying simStage through the tree. An index that two stages name differently is
    // genuinely ambiguous under that flattening, so it is reported instead of being
    // silently resolved to whichever stage came last.
    for (const PhysicalVolumeInfoSingleStage& stage : *volsHandle) {
      for (const auto& entry : stage) {
        const unsigned idx = entry.first.asUint();
        if (volumeNames_.size() <= idx) volumeNames_.resize(idx + 1);
        const std::string& name = entry.second.name();
        if (!volumeNames_[idx].empty() && volumeNames_[idx] != name)
          mf::LogWarning("STMDetProvenance")
            << "Volume index " << idx << " maps to both \"" << volumeNames_[idx]
            << "\" and \"" << name << "\" across simulation stages; keeping the first. "
            << "Treat the summary's volume names as indicative for this file.";
        else
          volumeNames_[idx] = name;
      }
    }
  }

  std::string STMDetProvenance::volumeLabel(unsigned index) const {
    const std::string idx = "idx " + std::to_string(index);
    if (index < volumeNames_.size() && !volumeNames_[index].empty())
      return volumeNames_[index] + " (" + idx + ")";
    return idx;
  }

  void STMDetProvenance::analyze(const art::Event& event) {
    auto const& steps = event.getProduct(stepsToken_);
    auto const& sims  = event.getProduct(simsToken_);
    if (steps.empty() || sims.empty()) return;

    runId_    = event.run();
    subRunId_ = event.subRun();
    eventId_  = event.event();
    std::set<int> hpgeTops, roots;   // distinct ids in this event

    // The set HPGeTree builds: particles that have at least one STMDet step. Reproduced
    // exactly (including that it is keyed on the SimParticle id) so hpgeTopParentId_ below
    // is what HPGeTree would actually compute, not an approximation of it.
    std::set<cet::map_vector_key> stmDetIds;
    for (const StepPointMC& step : steps)
      stmDetIds.insert(sims.at(step.trackId()).id());

    for (const StepPointMC& step : steps) {
      const SimParticle& p = sims.at(step.trackId());

      edep_     = step.ionizingEdep();
      stepTime_ = step.time();
      x_ = step.position().x();
      y_ = step.position().y();
      z_ = step.position().z();
      // HPGeTree keeps x < xBeamCentre for HPGe; recorded rather than applied, so both
      // detectors' steps are in one tree and the cut's own effect is measurable.
      isHPGeSide_ = (x_ < xBeamCentre_);

      trackId_          = static_cast<int>(p.id().asUint());
      pdgId_            = p.pdgId();
      // id() returns ProcessCode::enum_type; stored as an int so the branch is a plain
      // integer the summary's name table can be joined against offline.
      creationCode_     = static_cast<int>(p.creationCode().id());
      startVolumeIndex_ = p.startVolumeIndex();
      startGlobalTime_  = p.startGlobalTime();

      // HPGeTree treats primaries and truncated particles alike. A truncated particle has a
      // null parent but a non-mu2ePrimary creation code, i.e. its parent was dropped, which
      // should not happen with noCompression + keepNGenerations:-1.
      isPrimary_   = p.isPrimary();
      isTruncated_ = p.isTruncated();
      hasParent_   = p.hasParent();

      // ---- reproduce HPGeTree's topParent() decision ----
      // Walk up only while the parent is itself in stmDetIds, exactly as HPGeTree does.
      // stoppedAtShell: a resolvable parent exists but deposited nothing in a crystal, so
      // this particle becomes its own top parent.
      {
        const STMGenealogy::Walk top = STMGenealogy::topParentInDetector(p, stmDetIds);
        hpgeTopParentId_    = static_cast<int>(top.particle->id().asUint());
        hpgeTopParentPdg_   = top.particle->pdgId();
        hpgeStoppedAtShell_ = top.stoppedAtShell;
      }

      // ---- walk the FULL genealogy, with no STMDet restriction ----
      {
        const STMGenealogy::Walk root = STMGenealogy::genealogyRoot(p);
        const SimParticle* cur = root.particle;
        genDepth_             = root.depth;
        rootId_               = static_cast<int>(cur->id().asUint());
        rootPdg_              = cur->pdgId();
        rootStartTime_        = cur->startGlobalTime();
        rootStartVolumeIndex_ = cur->startVolumeIndex();
        rootIsPrimary_        = cur->isPrimary();
        rootIsTruncated_      = cur->isTruncated();
      }

      ttree->Fill();

      ++nSteps_;
      if (isPrimary_)   ++nPrimary_;
      if (isTruncated_) ++nTruncated_;
      if (hasParent_)   ++nWithParent_;
      ++startVolumeCounts_[startVolumeIndex_];
      if (hpgeStoppedAtShell_) {
        ++nStoppedAtShell_;
        ++shellStopVolumeCounts_[startVolumeIndex_];
      }
      ++rootPdgCounts_[rootPdg_];
      hpgeTops.insert(hpgeTopParentId_);
      roots.insert(rootId_);
    }

    ++nEventsWithSteps_;
    nHPGeTops_ += hpgeTops.size();
    nRoots_    += roots.size();
  }

  void STMDetProvenance::endJob() {
    mf::LogInfo log("STMDetProvenance");
    log << "========== STMDet provenance summary ==========\n"
        << "steps total                         : " << nSteps_ << "\n"
        << "  depositing particle is a PRIMARY  : " << nPrimary_ << "\n"
        << "  depositing particle is TRUNCATED  : " << nTruncated_ << "   <- expected 0\n"
        << "  depositing particle HAS a parent  : " << nWithParent_ << "\n"
        << "  HPGeTree walk stopped at a shell  : " << nStoppedAtShell_ << "\n"
        << "events with >=1 STMDet step         : " << nEventsWithSteps_ << "\n"
        << "distinct HPGeTree top parents       : " << nHPGeTops_
        << "   <- HPGeTree's entry count\n"
        << "distinct genealogy roots            : " << nRoots_
        << "   <- entry count if grouped by full ancestry\n";

    if (volumeNames_.empty())
      log << "\n(No volume table resolved - volumes are listed as raw indices. Supply "
          << "PhysicalVolumesTag to see names.)\n";

    log << "\nroot PDG counts (steps):\n";
    for (const auto& [pdg, n] : rootPdgCounts_)
      log << "  " << pdg << ": " << n << "\n";

    log << "\ndepositing particles were born in:\n";
    for (const auto& [idx, n] : startVolumeCounts_)
      log << "  " << volumeLabel(idx) << ": " << n << "\n";

    if (!shellStopVolumeCounts_.empty()) {
      log << "\nof those, the ones that made HPGeTree stop early were born in:\n";
      for (const auto& [idx, n] : shellStopVolumeCounts_)
        log << "  " << volumeLabel(idx) << ": " << n << "\n";
    }
    log << "===============================================\n";
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::STMDetProvenance)
