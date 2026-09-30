#include "Offline/DataProducts/inc/StageNormalization.hh"
#include "cetlib_except/exception.h"
#include <algorithm>

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
      throw cet::exception("BADINPUT")
        << "StageNormalization: chain depths conflict\n";
    if(other.fromOrigin() != fromOrigin_)
      throw cet::exception("BADINPUT")
        << "StageNormalization: one of these reaches the origin and the other "
        << "does not, so their generated counts are in different units\n";
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
      throw cet::exception("BADINPUT")
        << "StageNormalization: the upstream normalization is unset, so what "
        << "one drawn event represents is unknown. Either the pool carries no "
        << "StageNormalization, or the stated poolGenCount/poolEventCount were "
        << "not usable.\n";
    // fromOrigin is carried, not assumed: a chain composed from a pool whose
    // totals stop at an intermediate stage stays honest all the way down.
    return StageNormalization(double(nDraws) * upstream.perEvent(), nPassed,
                              upstream.nStages() + 1, upstream.fromOrigin());
  }

}
