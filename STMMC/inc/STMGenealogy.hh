#ifndef STMMC_inc_STMGenealogy_hh
#define STMMC_inc_STMGenealogy_hh
//
// SimParticle genealogy walks shared by STMDetTree and STMDetProvenance.
//
// topParentInDetector reproduces HPGeTree::topParent (HPGeTree_module.cc). HPGeTree and
// LaBrTree keep their own copies so their output does not change; if their walk changes,
// change this one with it.
//
// Yongyi Wu, Nov. 2026

#include <set>

#include "cetlib/map_vector.h"
#include "cetlib_except/exception.h"

#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e {
  namespace STMGenealogy {

    // Guard against cyclic parent links, which would otherwise hang the job.
    constexpr int kMaxDepth = 10000;

    // A parent that can be followed: set, and present in the collection.
    inline bool hasUsableParent(const SimParticle& p) {
      return p.hasParent() && p.parent().isAvailable();
    }

    inline void checkDepth(int depth, const SimParticle& start) {
      if (depth > kMaxDepth)
        throw cet::exception("STMGenealogy")
          << "Genealogy walk exceeded " << kMaxDepth << " generations starting from "
          << "SimParticle " << start.id().asUint() << "; the parent links are cyclic.";
    }

    struct Walk {
      const SimParticle* particle = nullptr;  // where the walk ended
      int  depth = 0;                         // generations walked up to it
      bool stoppedAtShell = false;            // topParentInDetector only: the parent is
                                              // usable but has no STMDet step
    };

    // HPGeTree's walk: up only while the parent itself has a step in STMDet. `stmDetIds`
    // holds the ids of every particle with at least one STMDet step in the event.
    inline Walk topParentInDetector(const SimParticle& p,
                                    const std::set<cet::map_vector_key>& stmDetIds) {
      Walk w;
      w.particle = &p;
      if (!hasUsableParent(p)) return w;
      if (!stmDetIds.count(p.parent()->id())) {
        w.stoppedAtShell = true;
        return w;
      }
      w.particle = p.parent().get();
      w.depth = 1;
      while (hasUsableParent(*w.particle) && stmDetIds.count(w.particle->parent()->id())) {
        w.particle = w.particle->parent().get();
        checkDepth(++w.depth, p);
      }
      return w;
    }

    // To the root of the full genealogy, with no STMDet restriction. Walked by hand rather
    // than via SimParticle::originParticle(), which follows selfParent() across
    // SIMULATION-STAGE boundaries rather than physics-genealogy ones.
    inline Walk genealogyRoot(const SimParticle& p) {
      Walk w;
      w.particle = &p;
      while (hasUsableParent(*w.particle)) {
        w.particle = w.particle->parent().get();
        checkDepth(++w.depth, p);
      }
      return w;
    }

  } // namespace STMGenealogy
} // namespace mu2e

#endif
