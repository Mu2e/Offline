// Info on the number of StraDigiMCs produced by a SimParticle.
//
// Andrei Gaponenko, 2014

#ifndef MCDataProducts_inc_SimParticleTrackerHitSummaryCollection_hh
#define MCDataProducts_inc_SimParticleTrackerHitSummaryCollection_hh

#include <vector>

namespace mu2e {

  class SimParticleTrackerHitSummary {
  public:
    SimParticleTrackerHitSummary(unsigned nPrincipal, unsigned nAll)
      : nPrincipal_(nPrincipal), nAll_(nAll)
    {}

    // Number of hits that were "made" by the particle
    unsigned nPrincipal()   const { return nPrincipal_; }

    // Number of hits to which the particle contributed
    unsigned nAll()    const { return nAll_; }

    // for persistency
    SimParticleTrackerHitSummary() : nPrincipal_(-1u), nAll_(-1u) {}

  private:
    unsigned nPrincipal_;
    unsigned nAll_;
  };

  typedef std::vector<SimParticleTrackerHitSummary> SimParticleTrackerHitSummaryCollection;

} // namespace mu2e

#endif /* MCDataProducts_inc_SimParticleTrackerHitSummaryCollection_hh */
