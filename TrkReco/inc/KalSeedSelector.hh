//
// Base class for selecting a KalSeed.  Subclasses must be an art tool
//
#ifndef TrkReco_KalSeedSelector_hh
#define TrkReco_KalSeedSelector_hh

namespace mu2e {
  class KalSeed;
  class KalSeedSelector {
  public:
    // the following returns true if the kseed passes the selection
    virtual bool select(KalSeed const& kseed) const = 0;
    // the following returns true if 'test' is better than 'current', false otherwise
    virtual bool isBetter(KalSeed const& current, KalSeed const& test) const = 0;
    virtual ~KalSeedSelector() noexcept = default;
  };
}
#endif
