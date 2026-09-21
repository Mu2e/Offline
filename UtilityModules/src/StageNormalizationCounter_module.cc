//
// Count the events one output stream actually keeps, and record with them how
// many events at the ORIGIN of the simulation chain those events represent,
// as a StageNormalization SubRun product.
//
// Place one instance at the END of each trigger path, so it sees exactly the
// events that reach that path's output:
//
//   targetStopPath : [ ..., TargetStopFilter, compressPVTargetStops, targetStopCounter ]
//
// It is a producer, not an analyzer, for that reason alone: analyzers run only
// in end_paths and so cannot see what a trigger path passed. It never rejects
// anything and produces no event-level data.
//
// Seeding the normalization
// -------------------------
// For a stage whose source reads events 1:1 (a generator or a filter job) the
// generated-equivalent count is the job's own GenEventCount, which is what
// genCountTag names. That is the default. Whether that count reaches the
// ORIGIN of the chain is genCountIsOrigin, and it is recorded with it.
//
// For a RESAMPLING stage it is not: the job's GenEventCount records the number
// of draws, and the resampled input arrives as a mixing secondary whose own
// GenEventCount never reaches the output. There the resampling mixer computes
// the normalization (it alone knows both the draw count and the pool it drew
// from) and writes it as a SubRun product; set upstreamTag to that product and
// this module reports it against its own passed count rather than re-deriving
// it.
//
// One instance per stream, and the label matters
// ----------------------------------------------
// SubRun products ignore an output's SelectEvents, so EVERY output stream of a
// multi-stream job carries EVERY counter's product -- exactly as every stream
// already carries every RandomPrescaleFilter's PrescaleFilterFraction. A
// consumer must therefore select by module label, and the label is the only
// thing that says which stream a count belongs to. Name them after the stream
// (targetStopCounter, polyStopCounter, ...), and never assume the sole product
// in a file is the relevant one.
//
// Original author: Michael MacKenzie, 2026
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "Offline/DataProducts/inc/StageNormalization.hh"
#include "Offline/MCDataProducts/inc/GenEventCount.hh"

#include <string>

namespace mu2e {

  class StageNormalizationCounter : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> genCountTag { Name("genCountTag"),
          Comment("SubRun GenEventCount of THIS job, used to seed the normalization "
                  "for a stage that reads its input 1:1."), "genCounter" };
        fhicl::OptionalAtom<art::InputTag> upstreamTag { Name("upstreamTag"),
          Comment("SubRun StageNormalization carrying the generated-equivalent count "
                  "already computed for this stage -- set this for a RESAMPLING "
                  "stage, where the mixer computes it and genCountTag would count "
                  "draws instead.") };
        fhicl::Atom<bool> countAsGenerated { Name("countAsGenerated"),
          Comment("This module IS the generated count: record nGenEquivalent = "
                  "the events it sees, rather than reading either tag. Put it "
                  "immediately after the generator, for a stage whose "
                  "GenEventCount cannot serve -- the cosmic S1 paths keep theirs "
                  "at the END of the path so CosmicLivetime scales correctly, and "
                  "a second GenEventCounter is refused by GenEventCounter itself. "
                  "A later counter then names THIS module in upstreamTag."), false };
        fhicl::Atom<bool> genCountIsOrigin { Name("genCountIsOrigin"),
          Comment("Whether the GenEventCount named by genCountTag counts events "
                  "at the ORIGIN of the chain. True for a generator, and for a "
                  "1:1 stage whose chain has never resampled, since the count "
                  "propagates. Set FALSE for a 1:1 stage over a resampled file "
                  "-- there GenEventCount is the resampling stage's draw count, "
                  "not the origin's, and calling it origin-referenced overstates "
                  "the chain by the upstream efficiency."), true };
        fhicl::Atom<int> diagLevel { Name("diagLevel"), Comment("Printout level"), 0 };
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit StageNormalizationCounter(const Parameters& conf);

      void produce(art::Event& event) override;
      void beginSubRun(art::SubRun& sr) override;
      void endSubRun(art::SubRun& sr) override;

    private:
      art::InputTag genCountTag_;
      art::InputTag upstreamTag_;
      bool resampled_;      // upstreamTag was supplied
      bool countAsGenerated_;
      bool genCountIsOrigin_;
      int  diagLevel_;
      uint64_t nPassed_;
  };

  StageNormalizationCounter::StageNormalizationCounter(const Parameters& conf)
    : art::EDProducer{conf}
    , genCountTag_{conf().genCountTag()}
    , resampled_{conf().upstreamTag(upstreamTag_)}
    , countAsGenerated_{conf().countAsGenerated()}
    , genCountIsOrigin_{conf().genCountIsOrigin()}
    , diagLevel_{conf().diagLevel()}
    , nPassed_{0}
  {
    if(countAsGenerated_ && resampled_) {
      throw cet::exception("BADCONFIG")
        << "StageNormalizationCounter: countAsGenerated says this module IS the "
        << "generated count, and upstreamTag says it inherits one. Give one or "
        << "the other.\n";
    }
    produces<StageNormalization, art::InSubRun>();
  }

  void StageNormalizationCounter::beginSubRun(art::SubRun&) {
    nPassed_ = 0;
  }

  void StageNormalizationCounter::produce(art::Event&) {
    // Reached only by events this path passed, which is the whole point of
    // sitting at the end of a trigger path.
    ++nPassed_;
  }

  void StageNormalizationCounter::endSubRun(art::SubRun& sr) {
    StageNormalization norm;

    if(countAsGenerated_) {
      // Placed right after the generator, so every event of the job reaches
      // it: what it counted IS the generated count, and the chain starts here.
      // This module IS the generator, so its count is the origin's.
      norm = StageNormalization(double(nPassed_), nPassed_, 1, true);
    } else if(resampled_) {
      // A resampling stage: the mixer already worked out the count, because
      // only it knows the pool the draws came from.
      auto h = sr.getHandle<StageNormalization>(upstreamTag_);
      if(h.isValid()) {
        // Depth and origin-reference both come from the mixer, which knows
        // what the pool's normalization actually reached.
        norm = StageNormalization(h->nGenEquivalent(), nPassed_, h->nStages(),
                                  h->fromOrigin());
      } else {
        // Deliberately left unseeded (nStages 0) rather than defaulted to
        // anything: an efficiency of 1, or of 0, would be indistinguishable
        // from a real measurement downstream. valid() says so explicitly.
        mf::LogWarning("StageNormalizationCounter")
          << "no StageNormalization '" << upstreamTag_ << "' in this SubRun; "
          << "writing an unseeded normalization (nStages=0). The resampling "
          << "mixer did not run, or is not configured to produce it.";
        norm = StageNormalization(0., nPassed_, 0, false);
      }
    } else {
      // A 1:1 stage: this job's own generated count, propagated from the
      // input file by RootOutput. genCountIsOrigin says whether that count
      // reaches the origin -- it does not for a file that was resampled.
      auto h = sr.getHandle<GenEventCount>(genCountTag_);
      if(!h.isValid()) {
        throw cet::exception("BADCONFIG")
          << "StageNormalizationCounter: no GenEventCount '" << genCountTag_
          << "' in this SubRun. Set genCountTag to the counter this job runs, "
          << "upstreamTag for a resampling stage, or countAsGenerated where "
          << "the generated count has to be taken here.\n";
      }
      norm = StageNormalization(double(h->count()), nPassed_, 1,
                                genCountIsOrigin_);
    }

    if(diagLevel_ > 0) {
      mf::LogInfo("StageNormalizationCounter")
        << "passed " << norm.nPassed() << " events representing "
        << norm.nGenEquivalent() << " generated events over "
        << norm.nStages() << " stage(s)"
        << (norm.fromOrigin() ? " back to the origin" : " (NOT reaching the origin)")
        << "; efficiency " << norm.efficiency();
    }

    sr.put(std::make_unique<StageNormalization>(norm), art::fullSubRun());
  }

}

DEFINE_ART_MODULE(mu2e::StageNormalizationCounter)
