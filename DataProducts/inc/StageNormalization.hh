//
// art product recording how many events at the ORIGIN of a simulation chain
// the events of one output stream represent, so that a rate per generated
// event (or per POT) can be recovered from the file alone, with no recourse
// to SAM metadata or a hand-supplied constant.
//
// The quantity carried is a pair:
//
//   nGenEquivalent : generated events at the origin of the chain that this
//                    stream's events collectively stand for
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
// input event stands for upstream.perEvent() = nGenEquivalent/nPassed origin
// events, which is the factor the resampler applies.
//
// Why not FilterFraction: its chain() requires upstream.nPassed() == nSeen(),
// a strict 1:1-consumption assumption that resampling violates by
// construction (nDraws is unrelated to the size of the pool drawn from, and
// events are reused). Its counts are also integers, while a resampled
// equivalent generally is not. The API below deliberately mirrors
// FilterFraction so the two can be unified if that constraint is ever
// relaxed.
//
// Caveat worth carrying wherever this number is used: nDraws * perEvent() is
// correct in EXPECTATION. When nDraws exceeds the size of the pool the draws
// reuse events, so the statistical error on anything derived from it is not
// naive-Poisson.
//
// Original author: Michael MacKenzie, 2026
//
#ifndef DataProducts_StageNormalization_hh
#define DataProducts_StageNormalization_hh
#include <cstdint>

namespace mu2e {
  class StageNormalization {
    public:
      StageNormalization(double nGenEquivalent, uint64_t nPassed, unsigned nStages = 1) :
        nGenEquivalent_(nGenEquivalent), nPassed_(nPassed), nStages_(nStages) {}
      StageNormalization(){}

      // accessors
      double   nGenEquivalent() const { return nGenEquivalent_; }
      uint64_t nPassed() const { return nPassed_; }
      // how many stages deep the chain this describes is; 0 means the
      // normalization was never seeded, which is NOT the same as an
      // efficiency of zero and must never be silently treated as one
      unsigned nStages() const { return nStages_; }
      // nPassed matters as much as the rest: a normalization recording no
      // events cannot say what one event represents, and perEvent() would
      // return 0, which a downstream resampler would otherwise accumulate as
      // if it were a measurement.
      bool valid() const {
        return nStages_ > 0 && nGenEquivalent_ > 0. && nPassed_ > 0;
      }

      // cumulative efficiency of every stage upstream of (and including) this one
      double efficiency() const {
        return (nGenEquivalent_ > 0.) ? double(nPassed_)/nGenEquivalent_ : 0.;
      }
      // origin-generated events each event of this stream stands for: the
      // factor a downstream resampler multiplies its draw count by
      double perEvent() const {
        return (nPassed_ > 0) ? nGenEquivalent_/double(nPassed_) : 0.;
      }

      // merge subruns (or files) of the SAME stream: both counts add
      StageNormalization& operator +=(StageNormalization const& other);
      StageNormalization  operator + (StageNormalization const& other) const;

      // compose with the stream this one resampled from: nDraws draws, each
      // standing for upstream.perEvent() origin events
      static StageNormalization resample(uint64_t nDraws, StageNormalization const& upstream,
                                         uint64_t nPassed);

    private:
      double   nGenEquivalent_ = 0.; // origin-stage generated events represented
      uint64_t nPassed_ = 0;         // events written to this stream
      unsigned nStages_ = 0;         // chain depth; 0 = never seeded
  };
}
#endif
