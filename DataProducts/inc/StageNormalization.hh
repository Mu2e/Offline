//
// art product recording how many events at the ORIGIN of a simulation chain
// the events of one output stream represent, so that a rate per generated
// event (or per POT) can be recovered from the file alone, with no recourse
// to SAM metadata or a hand-supplied constant.
//
// The quantity carried is a pair:
//
//   nGenEquivalent : generated events that this stream's events collectively
//                    stand for -- at the ORIGIN of the chain when
//                    fromOrigin() is true, otherwise at the stage nStages
//                    levels up
//   nPassed        : events actually written to this stream
//
// so the cumulative efficiency of everything upstream is nPassed/nGenEquivalent,
// and the pair composes along a chain of stages by simple bookkeeping:
//
//   generator stage   nGenEquivalent = GenEventCount of the job
//   1:1 filter stage  propagates unchanged (SubRun products are copied from
//                     the input file by RootOutput)
//   resampling stage  nGenEquivalent = nDraws * upstream.perEvent()
//
// The last line is why this product exists. A resampling job's source is
// EmptyEvent and its resampled input arrives as a mixing secondary, so the
// input's own GenEventCount does not reach the output: its GenEventCount
// records the number of DRAWS. Each draw stands for one input event, and one
// input event stands for upstream.perEvent() = nGenEquivalent/nPassed
// generated events, which is the factor the resampler applies.
//
// Why not FilterFraction: its chain() requires upstream.nPassed() == nSeen(),
// a strict 1:1-consumption assumption that resampling violates by
// construction (nDraws is unrelated to the size of the pool drawn from, and
// events are reused). Its counts are also integers, while a resampled
// equivalent generally is not. The API below deliberately mirrors
// FilterFraction so the two can be unified if that constraint is ever
// relaxed.
//
// Two caveats worth carrying wherever this number is used.
//
// nDraws * perEvent() is correct in EXPECTATION. When nDraws exceeds the size
// of the pool the draws reuse events, so the statistical error on anything
// derived from it is not naive-Poisson.
//
// And the forwarded form -- accumulating each drawn subrun's perEvent() -- is
// blind to pool subruns that kept NO events: a mix op is only ever called for
// a subrun something was drawn from, while the generated events of an empty
// subrun still belong in the denominator. The result reads LOW by about
// exp(-mu), mu being the mean events kept per pool subrun. It is nothing for a
// dense pool (no empty file in 5000 for Run1B target stops) and 2.2% for
// Run1B IPA stops, where 111 files of 5000 are empty. Stating the pool's
// totals instead is unbiased, because they cover every subrun whether or not
// it was drawn from -- so prefer that for a sparse pool, and note that the
// small undeclared test productions this product is aimed at are exactly
// where mu is smallest.
//
// Original author: Michael MacKenzie, 2026
//
#ifndef DataProducts_StageNormalization_hh
#define DataProducts_StageNormalization_hh
#include <cstdint>

namespace mu2e {
  class StageNormalization {
    public:
      // No defaults: a normalization that does not say how deep it goes, or
      // whether it reaches the origin, is what let a forwarded stage quietly
      // re-assert the origin after a bootstrap had denied it.
      StageNormalization(double nGenEquivalent, uint64_t nPassed,
                         unsigned nStages, bool fromOrigin) :
        nGenEquivalent_(nGenEquivalent), nPassed_(nPassed), nStages_(nStages),
        fromOrigin_(fromOrigin) {}
      StageNormalization(){}

      // accessors
      double   nGenEquivalent() const { return nGenEquivalent_; }
      uint64_t nPassed() const { return nPassed_; }
      // how many stages deep the chain this describes is; 0 means the
      // normalization was never seeded, which is NOT the same as an
      // efficiency of zero and must never be silently treated as one
      unsigned nStages() const { return nStages_; }
      // Whether nGenEquivalent counts events at the ORIGIN of the chain, or
      // only at some intermediate stage. False means the number is referenced
      // to whatever the stage nStages levels up generated, which is a real
      // measurement but NOT a rate per POT.
      //
      // It exists because the obvious source of a bootstrap's totals does not
      // reach the origin: SAM's dh.gencount is reset by a resampling stage to
      // that stage's draw count, so a stops sample records draws from the beam
      // sample, not protons.
      //
      // Whether a given pool is of that kind CAN be decided outside this
      // product, from SAM parentage and per FILE: a 1:1 stage -- splitter,
      // selector, concatenation -- has a file whose dh.gencount equals the sum
      // over its art parents, and a resampling stage does not. It cannot be
      // decided per DATASET, where a campaign may size a resampler to the same
      // round number as the stage above it (Run1Ban reads 2e9 at every level,
      // while per file TargetStops is 400000 against parents summing to
      // 10000000). prodtools does that walk when it fills the totals; nothing
      // inside a job can, which is why the answer is carried here.
      bool fromOrigin() const { return fromOrigin_; }
      // nPassed matters as much as the rest: a normalization recording no
      // events cannot say what one event represents, and perEvent() would
      // return 0, which a downstream resampler would otherwise accumulate as
      // if it were a measurement.
      bool valid() const {
        return nStages_ > 0 && nGenEquivalent_ > 0. && nPassed_ > 0;
      }

      // Efficiency of the nStages this normalization spans. It is a rate per
      // generated event at the origin only when fromOrigin() is true.
      double efficiency() const {
        return (nGenEquivalent_ > 0.) ? double(nPassed_)/nGenEquivalent_ : 0.;
      }
      // generated events each event of this stream stands for, at this
      // normalization's own depth: the factor a downstream resampler
      // multiplies its draw count by
      double perEvent() const {
        return (nPassed_ > 0) ? nGenEquivalent_/double(nPassed_) : 0.;
      }

      // merge subruns (or files) of the SAME stream: both counts add
      StageNormalization& operator +=(StageNormalization const& other);
      StageNormalization  operator + (StageNormalization const& other) const;

      // compose with the stream this one resampled from: nDraws draws, each
      // standing for upstream.perEvent() generated events
      static StageNormalization resample(uint64_t nDraws, StageNormalization const& upstream,
                                         uint64_t nPassed);

    private:
      double   nGenEquivalent_ = 0.; // generated events represented, at the
                                     // depth fromOrigin_/nStages_ describe
      uint64_t nPassed_ = 0;         // events written to this stream
      unsigned nStages_ = 0;         // stages spanned; 0 = never seeded
      bool     fromOrigin_ = false;  // does nGenEquivalent reach the origin?
  };
}
#endif
