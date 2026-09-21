#include "Offline/DataProducts/inc/StageNormalization.hh"
#include <algorithm>
#include <stdexcept>

namespace mu2e {

  // Merging subruns of the same stream: both counts are extensive, so both
  // add. The chain depth must agree -- adding a two-stage normalization to a
  // three-stage one would produce a number that is not an efficiency of
  // anything, and silently so.
  StageNormalization& StageNormalization::operator +=(StageNormalization const& other) {
    if(nStages_ == 0) {          // never seeded: adopt the other wholesale
      *this = other;
      return *this;
    }
    if(other.nStages() == 0) return *this;
    if(other.nStages() != nStages_)
      throw std::runtime_error("StageNormalization: chain depths conflict");
    nGenEquivalent_ += other.nGenEquivalent();
    nPassed_ += other.nPassed();
    return *this;
  }

  StageNormalization StageNormalization::operator + (StageNormalization const& other) const {
    auto retval = *this;
    retval += other;
    return retval;
  }

  // One more stage, entered by resampling rather than by reading events 1:1.
  // Each of the nDraws draws stands for upstream.perEvent() origin events,
  // whatever the size of the pool drawn from -- that is exactly the step
  // FilterFraction::chain() cannot express.
  StageNormalization StageNormalization::resample(uint64_t nDraws,
                                                  StageNormalization const& upstream,
                                                  uint64_t nPassed) {
    if(!upstream.valid())
      throw std::runtime_error("StageNormalization: upstream normalization is unset; "
                               "the resampled input carries no StageNormalization and "
                               "no GenEventCount fallback was available");
    return StageNormalization(double(nDraws) * upstream.perEvent(), nPassed,
                              upstream.nStages() + 1);
  }

}
